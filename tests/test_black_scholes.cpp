#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "black_scholes.hpp"
#include "binomial_tree.hpp"

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