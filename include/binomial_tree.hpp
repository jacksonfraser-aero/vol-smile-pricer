#pragma once

class BinomialTree {
public:
    BinomialTree(double S, double K, double T, double r, double sigma, int steps, bool isCall);

    double price() const;

private:
    double S_, K_, T_, r_, sigma_;
    int steps_;
    bool isCall_;
};