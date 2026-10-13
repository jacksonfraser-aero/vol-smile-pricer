#include "catch.hpp"
#include "binomial_tree.hpp"
#include <algorithm>
#include <cmath>

// Reference Black-Scholes prices for S=K=100, T=1, r=0.05, sigma=0.2 (Hull, standard ATM example).
static const double kBsAtmCall = 10.450583572185565;
static const double kBsAtmPut = 5.573526022256971;

TEST_CASE("Binomial European prices converge to Black-Scholes as steps increase", "[binomial][convergence]") {
    const double S = 100.0, K = 100.0, T = 1.0, r = 0.05, sigma = 0.2;

    for (bool isCall : {true, false}) {
        double reference = isCall ? kBsAtmCall : kBsAtmPut;
        auto error = [&](int steps) {
            return std::abs(BinomialTree(S, K, T, r, sigma, steps, isCall).price() - reference);
        };

        double e50 = error(50), e500 = error(500), e2000 = error(2000);
        INFO((isCall ? "call" : "put") << " errors: N=50 " << e50 << ", N=500 " << e500 << ", N=2000 " << e2000);

        REQUIRE(e500 < e50);
        REQUIRE(e2000 < e500);
        REQUIRE(e2000 < 1e-3);
    }
}

TEST_CASE("Put-call parity holds for binomial European prices", "[binomial][parity]") {
    const double S = 100.0, T = 0.75, r = 0.04, sigma = 0.3;
    const int steps = 300;

    for (double K : {70.0, 100.0, 130.0}) {
        double call = BinomialTree(S, K, T, r, sigma, steps, true).price();
        double put = BinomialTree(S, K, T, r, sigma, steps, false).price();
        INFO("K = " << K << ", C - P = " << call - put);
        REQUIRE(call - put == Approx(S - K * std::exp(-r * T)).margin(1e-9));
    }
}

TEST_CASE("Binomial prices respect no-arbitrage bounds", "[binomial][bounds]") {
    const double S = 100.0, r = 0.05;
    const int steps = 200;

    for (double K : {20.0, 80.0, 100.0, 120.0, 500.0}) {
        const double tol = 1e-12 * std::max(S, K); // floating-point slack, scaled to price magnitude
        for (double T : {0.01, 0.5, 5.0}) {
            for (double sigma : {0.05, 0.3, 1.0}) {
                double call = BinomialTree(S, K, T, r, sigma, steps, true).price();
                double put = BinomialTree(S, K, T, r, sigma, steps, false).price();
                double pvK = K * std::exp(-r * T);
                INFO("K=" << K << " T=" << T << " sigma=" << sigma << " C=" << call << " P=" << put);

                REQUIRE(call >= std::max(S - pvK, 0.0) - tol);
                REQUIRE(call <= S + tol);
                REQUIRE(put >= std::max(pvK - S, 0.0) - tol);
                REQUIRE(put <= pvK + tol);
            }
        }
    }
}
