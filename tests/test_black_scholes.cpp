#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "black_scholes.hpp"
#include "binomial_tree.hpp"
#include <algorithm>

TEST_CASE("Binomial tree call roughly matches BS pricer", "[binomial]") {
    BinomialTree tree(42.0, 40.0, 0.5, 0.10, 0.20, 500, true);
    REQUIRE(tree.price() == Approx(4.76).epsilon(0.02));
}

TEST_CASE("Binomial tree put roughly matches BS pricer", "[binomial]") {
    BinomialTree tree(42.0, 40.0, 0.5, 0.10, 0.20, 500, false);
    REQUIRE(tree.price() == Approx(0.81).epsilon(0.03));
}

TEST_CASE("Call and put price match known textbook values", "[pricing]") {
    BlackScholes bs(42.0, 40.0, 0.5, 0.10, 0.20);

    REQUIRE(bs.callPrice() == Approx(4.76).epsilon(0.01));
    REQUIRE(bs.putPrice() == Approx(0.81).epsilon(0.01));
}

TEST_CASE("Put-call parity holds", "[pricing]") {
    BlackScholes bs(100.0, 100.0, 1.0, 0.05, 0.2);

    double call = bs.callPrice();
    double put = bs.putPrice();
    double lhs = call - put;
    double rhs = 100.0 - 100.0 * std::exp(-0.05 * 1.0);

    REQUIRE(lhs == Approx(rhs).epsilon(0.0001));
}

TEST_CASE("Deep ITM call approaches intrinsic value", "[pricing]") {
    BlackScholes bs(200.0, 50.0, 0.01, 0.05, 0.2);
    REQUIRE(bs.callPrice() == Approx(150.0).epsilon(0.05));
}

TEST_CASE("Greeks match known values for standard case", "[greeks]") {
    BlackScholes bs(42.0, 40.0, 0.5, 0.10, 0.20);

    REQUIRE(bs.delta(true) == Approx(0.7791).epsilon(0.01));
    REQUIRE(bs.delta(false) == Approx(-0.2209).epsilon(0.01));
    REQUIRE(bs.gamma() == Approx(0.0498).epsilon(0.02));
    REQUIRE(bs.vega() == Approx(8.813).epsilon(0.02));
}

TEST_CASE("Call delta approaches 1 deep ITM", "[greeks]") {
    BlackScholes bs(200.0, 50.0, 0.5, 0.05, 0.2);
    REQUIRE(bs.delta(true) == Approx(1.0).epsilon(0.01));
}

TEST_CASE("Put delta approaches -1 deep ITM", "[greeks]") {
    BlackScholes bs(20.0, 200.0, 0.5, 0.05, 0.2);
    REQUIRE(bs.delta(false) == Approx(-1.0).epsilon(0.01));
}
TEST_CASE("Put-call parity holds for ITM and OTM strikes", "[pricing][parity]") {
    const double S = 100.0, T = 0.75, r = 0.04, sigma = 0.3;

    for (double K : {70.0, 130.0}) {
        BlackScholes bs(S, K, T, r, sigma);
        INFO("K = " << K);
        REQUIRE(bs.callPrice() - bs.putPrice() == Approx(S - K * std::exp(-r * T)).margin(1e-10));
    }
}

TEST_CASE("Prices respect no-arbitrage bounds", "[pricing][bounds]") {
    const double S = 100.0, r = 0.05;

    for (double K : {20.0, 80.0, 100.0, 120.0, 500.0}) {
        const double tol = 1e-12 * std::max(S, K); // floating-point slack, scaled to price magnitude
        for (double T : {0.01, 0.5, 5.0}) {
            for (double sigma : {0.05, 0.3, 1.0}) {
                BlackScholes bs(S, K, T, r, sigma);
                double call = bs.callPrice(), put = bs.putPrice();
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

namespace {
struct Inputs { double S, K, T, r, sigma; };

double price(const Inputs& in, bool isCall) {
    BlackScholes bs(in.S, in.K, in.T, in.r, in.sigma);
    return isCall ? bs.callPrice() : bs.putPrice();
}

// Central difference of the BS price with respect to one input.
template <typename Bump>
double centralDiff(const Inputs& in, bool isCall, double h, Bump bump) {
    Inputs up = in, down = in;
    bump(up, h);
    bump(down, -h);
    return (price(up, isCall) - price(down, isCall)) / (2.0 * h);
}
}

TEST_CASE("Greeks match central finite differences of the BS price", "[greeks][finite-difference]") {
    const double S = 100.0, T = 0.5, r = 0.05, sigma = 0.25;
    // Bump sizes: relative for spot, absolute for the rest.
    // Tolerance is sized to the O(h^2) truncation error of a central difference at these bumps
    // (largest observed gap ~1.2e-6, on delta), not to floating-point noise.
    const double tol = 1e-5;
    const double hS = 1e-3 * S, hSigma = 1e-4, hR = 1e-4, hT = 1e-4;

    for (double K : {80.0, 100.0, 120.0}) {   // call ITM / ATM / OTM (put the reverse)
        const Inputs in{S, K, T, r, sigma};
        BlackScholes bs(S, K, T, r, sigma);

        for (bool isCall : {true, false}) {
            INFO("K=" << K << (isCall ? " call" : " put"));

            double fdDelta = centralDiff(in, isCall, hS, [](Inputs& x, double h) { x.S += h; });
            double fdVega = centralDiff(in, isCall, hSigma, [](Inputs& x, double h) { x.sigma += h; });
            double fdRho = centralDiff(in, isCall, hR, [](Inputs& x, double h) { x.r += h; });
            // theta is dV/dt in calendar time, i.e. -dV/dT
            double fdTheta = -centralDiff(in, isCall, hT, [](Inputs& x, double h) { x.T += h; });

            REQUIRE(bs.delta(isCall) == Approx(fdDelta).margin(tol));
            REQUIRE(bs.vega() == Approx(fdVega).margin(tol));
            REQUIRE(bs.rho(isCall) == Approx(fdRho).margin(tol));
            REQUIRE(bs.theta(isCall) == Approx(fdTheta).margin(tol));

            Inputs up = in, down = in;
            up.S += hS;
            down.S -= hS;
            double fdGamma = (price(up, isCall) - 2.0 * price(in, isCall) + price(down, isCall)) / (hS * hS);
            REQUIRE(bs.gamma() == Approx(fdGamma).margin(tol));
        }
    }
}
