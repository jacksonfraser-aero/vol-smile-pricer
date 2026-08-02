#pragma once

class BlackScholes {
public:
    BlackScholes(double S, double K, double T, double r, double sigma);

    double callPrice() const;
    double putPrice() const;

    double delta(bool isCall) const;
    double gamma() const;
    double vega() const;
    double theta(bool isCall) const;
    double rho(bool isCall) const;

private:
    double S_, K_, T_, r_, sigma_;

    double d1() const;
    double d2() const;
    static double normCDF(double x);
    static double normPDF(double x);
};