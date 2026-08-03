#include "delta_hedge.hpp"
#include "black_scholes.hpp"
#include <cmath>
#include <iostream>

DeltaHedger::DeltaHedger(double K, double T_total, double r, double sigma, bool isCall)
    : K_(K), T_total_(T_total), r_(r), sigma_(sigma), isCall_(isCall) {}

std::vector<HedgeStep> DeltaHedger::simulate(const std::vector<double>& prices, int stepsPerDay) {
    std::vector<HedgeStep> results;
    int n = prices.size();

    // sell one option at t=0, collect the premium
    BlackScholes bs0(prices[0], K_, T_total_, r_, sigma_);
    double premium = isCall_ ? bs0.callPrice() : bs0.putPrice();

    double cash = premium;   // cash from selling the option
    double heldShares = 0.0;
    double dt = (T_total_) / (n - 1);  // time per step, in years

    for (int i = 0; i < n; ++i) {
        double timeRemaining = T_total_ - i * dt;
        if (timeRemaining < 1e-6) timeRemaining = 1e-6; // avoid div-by-zero at expiry

        BlackScholes bs(prices[i], K_, timeRemaining, r_, sigma_);
        double delta = bs.delta(isCall_);

        // rebalance: buy/sell shares to match new delta
        double shareChange = delta - heldShares;
        cash -= shareChange * prices[i];   // pay for shares bought (or receive for sold)
        heldShares = delta;

        // accrue interest on cash at the risk-free rate for one step
        if (i > 0) {
            cash *= std::exp(r_ * dt);
        }

        double portfolioValue = cash + heldShares * prices[i];

        results.push_back({"", prices[i], delta, cash, portfolioValue});
    }

    // at expiry: settle the option payoff, and unwind the hedge
    double ST = prices.back();
    double payoff = isCall_ ? std::max(ST - K_, 0.0) : std::max(K_ - ST, 0.0);

    double finalCash = results.back().cashPosition;
    double finalShares = results.back().delta;

    // close out: sell remaining shares, pay the option payoff (we're short the option)
    double finalPortfolio = finalCash + finalShares * ST - payoff;

    pnl_ = finalPortfolio; // started with 0 capital beyond the premium received

    return results;
}

double DeltaHedger::finalPnL() const {
    return pnl_;
}