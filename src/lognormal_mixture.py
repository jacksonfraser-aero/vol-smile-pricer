"""Lognormal-mixture model at a single maturity.

Reference: D. Brigo and F. Mercurio, "Lognormal-Mixture Dynamics and Calibration to Market
Volatility Smiles", International Journal of Theoretical and Applied Finance 5(4), 2002.

The terminal density of S_T is a mixture of lognormal densities, each with the common drift r
and its own volatility sigma_i, weighted by lambda_i > 0 with sum(lambda_i) = 1. European
option prices are therefore the same convex combination of Black-Scholes prices.
"""
from implied_vol import bs_call_price


def mixture_call_price(S, K, T, r, lambdas, sigmas):
    return sum(lam * bs_call_price(S, K, T, r, sigma) for lam, sigma in zip(lambdas, sigmas))
