#include "delta_hedge.hpp"
#include "gbm_path.hpp"
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// Monte Carlo version of run_hedge: flat-vol vs smile-vol hedging over many simulated paths.
// Both hedgers sell at the smile-implied vol, so the P&L difference reflects hedge ratios only.

double firstClose(const std::string& path) {
    std::ifstream file(path);
    std::string line, dateStr, closeStr;
    std::getline(file, line); // header
    std::getline(file, line);
    std::stringstream ss(line);
    std::getline(ss, dateStr, ',');
    std::getline(ss, closeStr, ',');
    return std::stod(closeStr);
}

struct Stats {
    double mean, sd;
};

Stats sampleStats(const std::vector<double>& xs) {
    double sum = 0.0;
    for (double x : xs) sum += x;
    double mean = sum / xs.size();
    double ss = 0.0;
    for (double x : xs) ss += (x - mean) * (x - mean);
    return {mean, std::sqrt(ss / (xs.size() - 1))};
}

int main() {
    const double S0 = firstClose("data/price_history.csv");
    const double K = 350.0;
    const double T = 26.0 / 365.0;
    const double r = 0.05;
    const int steps = 18;  // 19 prices, 19 rebalances, matching the historical run
    const int paths = 2000;
    const std::uint32_t seed = 20260708;
    const double flatSigma = 0.27;
    const double smileSigma = smileImpliedVol(K);

    std::cout << "S0 = " << S0 << ", K = " << K << ", flat vol = " << flatSigma
              << ", smile-implied vol = " << smileSigma << ", " << paths << " paths\n\n";

    struct Scenario {
        std::string name;
        double realisedVol;
    };
    const std::vector<Scenario> scenarios = {
        {"A: realised = smile-implied", smileSigma},
        {"B: realised = flat", flatSigma},
    };

    std::ofstream out("data/mc_hedge_comparison.csv");
    out << "scenario,realised_vol,flat_mean,flat_sd,smile_mean,smile_sd,diff_mean,diff_se,sd_ratio_smile_over_flat\n";

    std::cout << std::fixed << std::setprecision(4);
    for (const auto& sc : scenarios) {
        std::vector<double> flatPnL, smilePnL, diff;
        for (int k = 0; k < paths; ++k) {
            // Same path for both hedgers within a scenario.
            const auto path = gbmPath(S0, r, sc.realisedVol, T, steps, seed + k);

            DeltaHedger flat(K, T, r, flatSigma, true, smileSigma);
            DeltaHedger smile(K, T, r, smileSigma, true, smileSigma);
            flat.simulate(path);
            smile.simulate(path);

            flatPnL.push_back(flat.finalPnL());
            smilePnL.push_back(smile.finalPnL());
            diff.push_back(smile.finalPnL() - flat.finalPnL());
        }

        Stats f = sampleStats(flatPnL), s = sampleStats(smilePnL), d = sampleStats(diff);
        double diffSE = d.sd / std::sqrt(static_cast<double>(paths));
        double sdRatio = s.sd / f.sd;

        std::cout << sc.name << " (realised vol " << sc.realisedVol << ")\n"
                  << "  flat  P&L: mean " << f.mean << ", sd " << f.sd << "\n"
                  << "  smile P&L: mean " << s.mean << ", sd " << s.sd << "\n"
                  << "  smile - flat: mean " << d.mean << ", se " << diffSE << "\n"
                  << "  sd ratio (smile / flat): " << sdRatio << "\n\n";

        out << std::setprecision(6) << '"' << sc.name << "\"," << sc.realisedVol << ","
            << f.mean << "," << f.sd << "," << s.mean << "," << s.sd << ","
            << d.mean << "," << diffSE << "," << sdRatio << "\n";
        std::cout << std::setprecision(4);
    }

    std::cout << "Wrote data/mc_hedge_comparison.csv\n";
    return 0;
}
