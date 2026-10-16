#pragma once
#include <vector>
#include <string>
#include <optional>

struct HedgeStep {
    std::string date;
    double spot;
    double delta;
    double cashPosition;
    double portfolioValue;
};

class DeltaHedger {
public:
    // sigma drives the hedge ratios; premiumSigma (defaults to sigma) prices the option sold at t=0.
    DeltaHedger(double K, double T_total, double r, double sigma, bool isCall,
                std::optional<double> premiumSigma = std::nullopt);
    // Rebalances once per input price: frequency is determined by the input path's spacing.
    std::vector<HedgeStep> simulate(const std::vector<double>& prices);
    double finalPnL() const;

private:
    double K_, T_total_, r_, sigma_, premiumSigma_;
    bool isCall_;
    double pnl_ = 0.0;
};

double smileImpliedVol(double strike);