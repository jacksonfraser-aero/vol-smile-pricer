#include "black_scholes.hpp"
#include "binomial_tree.hpp"
#include <iostream>
#include <fstream>
#include <vector>

int main() {
    double S = 42.0, K = 40.0, T = 0.5, r = 0.10, sigma = 0.20;

    BlackScholes bs(S, K, T, r, sigma);
    double bsCall = bs.callPrice();

    std::vector<int> stepCounts = {1, 2, 5, 10, 25, 50, 100, 200, 500, 1000};

    std::ofstream out("data/convergence.csv");
    out << "steps,binomial_price,bs_price,abs_diff\n";

    for (int steps : stepCounts) {
        BinomialTree tree(S, K, T, r, sigma, steps, true);
        double treePrice = tree.price();
        double diff = std::abs(treePrice - bsCall);

        out << steps << "," << treePrice << "," << bsCall << "," << diff << "\n";
        std::cout << "steps=" << steps << " tree=" << treePrice
                   << " bs=" << bsCall << " diff=" << diff << "\n";
    }

    out.close();
    std::cout << "Wrote data/convergence.csv\n";
    return 0;
}