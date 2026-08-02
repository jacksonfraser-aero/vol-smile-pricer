#include "black_scholes.hpp"
#include <cmath>

BlackScholes::BlackScholes(double S, double K, double T, double r, double sigma)
    : S_(S), K_(K), T_(T), r_(r), sigma_(sigma) {}

double BlackScholes::d1() const {
    return (std::log(S_ / K_) + (r_ + 0.5 * sigma_ * sigma_) * T_) / (sigma_ * std::sqrt(T_));
}

double BlackScholes::d2() const {
    return d1() - sigma_ * std::sqrt(T_);
}

double BlackScholes::normCDF(double x) {
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

double BlackScholes::callPrice() const {
    return S_ * normCDF(d1()) - K_ * std::exp(-r_ * T_) * normCDF(d2());
}

double BlackScholes::putPrice() const {
    return K_ * std::exp(-r_ * T_) * normCDF(-d2()) - S_ * normCDF(-d1());
}

double BlackScholes::normPDF(double x) {
    return (1.0 / std::sqrt(2.0 * M_PI)) * std::exp(-0.5 * x * x);
}

double BlackScholes::delta(bool isCall) const {
    if (isCall) {
        return normCDF(d1());
    }
    return normCDF(d1()) - 1.0;
}

double BlackScholes::gamma() const {
    return normPDF(d1()) / (S_ * sigma_ * std::sqrt(T_));
}

double BlackScholes::vega() const {
    return S_ * normPDF(d1()) * std::sqrt(T_);
}

double BlackScholes::theta(bool isCall) const {
    double term1 = -(S_ * normPDF(d1()) * sigma_) / (2 * std::sqrt(T_));
    if (isCall) {
        double term2 = r_ * K_ * std::exp(-r_ * T_) * normCDF(d2());
        return term1 - term2;
    }
    double term2 = r_ * K_ * std::exp(-r_ * T_) * normCDF(-d2());
    return term1 + term2;
}

double BlackScholes::rho(bool isCall) const {
    if (isCall) {
        return K_ * T_ * std::exp(-r_ * T_) * normCDF(d2());
    }
    return -K_ * T_ * std::exp(-r_ * T_) * normCDF(-d2());
}