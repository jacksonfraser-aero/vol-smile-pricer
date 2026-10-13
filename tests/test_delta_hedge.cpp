#include "catch.hpp"
#include "delta_hedge.hpp"
#include <cmath>
#include <cstdint>
#include <random>
#include <vector>

namespace {
// Geometric Brownian motion path with `steps` increments over [0, T], so steps + 1 prices.
// Uses mt19937 (sequence fixed by the standard) with a hand-rolled Box-Muller transform, because
// std::normal_distribution's output is implementation-defined; this keeps tests deterministic
// across compilers and standard libraries.
std::vector<double> gbmPath(double S0, double mu, double sigma, double T, int steps, std::uint32_t seed) {
    std::mt19937 gen(seed);
    auto uniform = [&gen] { return (static_cast<double>(gen()) + 0.5) / 4294967296.0; }; // (0, 1)

    const double dt = T / steps;
    std::vector<double> path{S0};
    while (static_cast<int>(path.size()) < steps + 1) {
        double u1 = uniform(), u2 = uniform();
        double radius = std::sqrt(-2.0 * std::log(u1));
        for (double z : {radius * std::cos(2.0 * M_PI * u2), radius * std::sin(2.0 * M_PI * u2)}) {
            if (static_cast<int>(path.size()) == steps + 1) break;
            path.push_back(path.back() * std::exp((mu - 0.5 * sigma * sigma) * dt + sigma * std::sqrt(dt) * z));
        }
    }
    return path;
}
}

TEST_CASE("Hedge is self-financing at every rebalance", "[hedge][self-financing]") {
    const double S0 = 100.0, K = 100.0, T = 0.25, r = 0.05, sigma = 0.2;
    const int steps = 63;
    const double dt = T / steps;

    for (bool isCall : {true, false}) {
        DeltaHedger hedger(K, T, r, sigma, isCall);
        auto rows = hedger.simulate(gbmPath(S0, r, sigma, T, steps, 42));
        REQUIRE(rows.size() == static_cast<size_t>(steps + 1));

        for (size_t t = 1; t < rows.size(); ++t) {
            double expected = rows[t - 1].cashPosition * std::exp(r * dt)
                            - (rows[t].delta - rows[t - 1].delta) * rows[t].spot;
            INFO((isCall ? "call" : "put") << " t=" << t << " cash=" << rows[t].cashPosition
                 << " expected=" << expected << " diff=" << rows[t].cashPosition - expected);
            REQUIRE(rows[t].cashPosition == Approx(expected).epsilon(0).margin(1e-9));
        }
    }
}
