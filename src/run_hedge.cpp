#include "delta_hedge.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

std::vector<double> loadPrices(const std::string& path) {
    std::vector<double> prices;
    std::ifstream file(path);
    std::string line;
    std::getline(file, line); // skip header

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

    double K = 310.0;       
    double T = 26.0 / 365.0; 
    double r = 0.05;
    double sigma = 0.27;    

    DeltaHedger hedger(K, T, r, sigma, true);
    auto results = hedger.simulate(prices);

    std::ofstream out("data/hedge_flatvol.csv");
    out << "step,spot,delta,cash,portfolio_value\n";
    for (size_t i = 0; i < results.size(); ++i) {
        out << i << "," << results[i].spot << "," << results[i].delta << ","
            << results[i].cashPosition << "," << results[i].portfolioValue << "\n";
    }

    std::cout << "Final P&L (flat vol hedge): " << hedger.finalPnL() << "\n";
    std::cout << "Wrote data/hedge_flatvol.csv\n";
    return 0;
}