#include "delta_hedge.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

std::vector<double> loadPrices(const std::string& path) {
    std::vector<double> prices;
    std::ifstream file(path);
    std::string line;
    std::getline(file, line);
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string dateStr, closeStr;
        std::getline(ss, dateStr, ',');
        std::getline(ss, closeStr, ',');
        prices.push_back(std::stod(closeStr));
    }
    return prices;
}

int main() {
    std::vector<double> prices = loadPrices("data/price_history.csv");
    std::cout << "Loaded " << prices.size() << " price points\n";

    double K = 350.0;
    double T = 26.0 / 365.0;
    double r = 0.05;
    double flatSigma = 0.27;
    double smileSigma = smileImpliedVol(K);

    std::cout << "Flat vol: " << flatSigma << ", smile-implied vol at K=" << K << ": " << smileSigma << "\n";

    DeltaHedger hedgerFlat(K, T, r, flatSigma, true);
    auto resultsFlat = hedgerFlat.simulate(prices);

    DeltaHedger hedgerSmile(K, T, r, smileSigma, true);
    auto resultsSmile = hedgerSmile.simulate(prices);

    std::ofstream out("data/hedge_comparison.csv");
    out << "step,spot,flat_delta,flat_portfolio,smile_delta,smile_portfolio\n";
    for (size_t i = 0; i < resultsFlat.size(); ++i) {
        out << i << "," << resultsFlat[i].spot << ","
            << resultsFlat[i].delta << "," << resultsFlat[i].portfolioValue << ","
            << resultsSmile[i].delta << "," << resultsSmile[i].portfolioValue << "\n";
    }

    std::cout << "Final P&L (flat vol):   " << hedgerFlat.finalPnL() << "\n";
    std::cout << "Final P&L (smile vol):  " << hedgerSmile.finalPnL() << "\n";
    std::cout << "Wrote data/hedge_comparison.csv\n";
    return 0;
}