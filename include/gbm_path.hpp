#pragma once
#include <cmath>
#include <cstdint>
#include <random>
#include <vector>

// Geometric Brownian motion path with `steps` increments over [0, T], so steps + 1 prices.
// Uses mt19937 (sequence fixed by the standard) with a hand-rolled Box-Muller transform, because
// std::normal_distribution's output is implementation-defined; this keeps paths deterministic
// across compilers and standard libraries.
inline std::vector<double> gbmPath(double S0, double mu, double sigma, double T, int steps, std::uint32_t seed) {
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
