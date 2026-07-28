#pragma once

class BlackScholes {
public:
    BlackScholes(double S, double K, double T, double r, double sigma);

    double callPrice() const;
    double putPrice() const;

private:
    double S_;      // spot price
    double K_;      // strike
    double T_;      // time to expiry (years)
    double r_;       // risk-free rate
    double sigma_;   // volatility

    double d1() const;
    double d2() const;
    static double normCDF(double x);
};