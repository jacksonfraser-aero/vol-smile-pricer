#pragma once
#include <vector>
#include <string>

struct HedgeStep {
    std::string date;
    double spot;
    double delta;
    double cashPosition;
    double portfolioValue;
};

class DeltaHedger {
public:
    DeltaHedger(double K, double T_total, double r, double sigma, bool isCall);
    std::vector<HedgeStep> simulate(const std::vector<double>& prices, int stepsPerDay = 1);
    double finalPnL() const;

private:
    double K_, T_total_, r_, sigma_;
    bool isCall_;
    double pnl_ = 0.0;
};

double smileImpliedVol(double strike);