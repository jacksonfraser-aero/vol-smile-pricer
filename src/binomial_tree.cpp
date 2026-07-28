#include "binomial_tree.hpp"
#include <cmath>
#include <vector>
#include <algorithm>

BinomialTree::BinomialTree(double S, double K, double T, double r, double sigma, int steps, bool isCall)
    : S_(S), K_(K), T_(T), r_(r), sigma_(sigma), steps_(steps), isCall_(isCall) {}

double BinomialTree::price() const {
    double dt = T_ / steps_;
    double u = std::exp(sigma_ * std::sqrt(dt));
    double d = 1.0 / u;
    double p = (std::exp(r_ * dt) - d) / (u - d);
    double discount = std::exp(-r_ * dt);

    // terminal payoffs at maturity
    std::vector<double> values(steps_ + 1);
    for (int i = 0; i <= steps_; ++i) {
        double ST = S_ * std::pow(u, steps_ - i) * std::pow(d, i);
        values[i] = isCall_ ? std::max(ST - K_, 0.0) : std::max(K_ - ST, 0.0);
    }

    // step backward through the tree
    for (int step = steps_ - 1; step >= 0; --step) {
        for (int i = 0; i <= step; ++i) {
            values[i] = discount * (p * values[i] + (1 - p) * values[i + 1]);
        }
    }

    return values[0];
}