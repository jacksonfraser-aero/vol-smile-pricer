#include "catch.hpp"
#include "binomial_tree.hpp"
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
