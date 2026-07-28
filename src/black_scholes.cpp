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