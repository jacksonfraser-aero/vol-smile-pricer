"""Lognormal-mixture model at a single maturity.

Reference: D. Brigo and F. Mercurio, "Lognormal-Mixture Dynamics and Calibration to Market
Volatility Smiles", International Journal of Theoretical and Applied Finance 5(4), 2002.

The terminal density of S_T is a mixture of lognormal densities, each with the common drift r
and its own volatility sigma_i, weighted by lambda_i > 0 with sum(lambda_i) = 1. European
option prices are therefore the same convex combination of Black-Scholes prices.
"""
import numpy as np
from scipy.stats import norm

from implied_vol import bs_call_price


def _validate(lambdas, sigmas):
    if len(lambdas) != len(sigmas):
        raise ValueError(f"{len(lambdas)} lambdas but {len(sigmas)} sigmas")
    if any(lam <= 0 for lam in lambdas):
        raise ValueError(f"lambdas must be positive, got {list(lambdas)}")
    if abs(sum(lambdas) - 1.0) > 1e-12:
        raise ValueError(f"lambdas must sum to 1, got sum {sum(lambdas)!r}")
    if any(sigma <= 0 for sigma in sigmas):
        raise ValueError(f"sigmas must be positive, got {list(sigmas)}")


def mixture_call_price(S, K, T, r, lambdas, sigmas):
    _validate(lambdas, sigmas)
    return sum(lam * bs_call_price(S, K, T, r, sigma) for lam, sigma in zip(lambdas, sigmas))


def _bs_call_delta(S, K, T, r, sigma):
    d1 = (np.log(S / K) + (r + 0.5 * sigma**2) * T) / (sigma * np.sqrt(T))
    return norm.cdf(d1)


def mixture_call_delta(S, K, T, r, lambdas, sigmas):
    _validate(lambdas, sigmas)
    return sum(lam * _bs_call_delta(S, K, T, r, sigma) for lam, sigma in zip(lambdas, sigmas))


def _bs_put_price(S, K, T, r, sigma):
    d1 = (np.log(S / K) + (r + 0.5 * sigma**2) * T) / (sigma * np.sqrt(T))
    d2 = d1 - sigma * np.sqrt(T)
    return K * np.exp(-r * T) * norm.cdf(-d2) - S * norm.cdf(-d1)


def mixture_put_price(S, K, T, r, lambdas, sigmas):
    _validate(lambdas, sigmas)
    return sum(lam * _bs_put_price(S, K, T, r, sigma) for lam, sigma in zip(lambdas, sigmas))


def mixture_density(S, K, T, r, lambdas, sigmas):
    # Risk-neutral density of S_T at K: each component is lognormal with
    # ln S_T ~ N(ln S + (r - sigma_i^2 / 2) T, sigma_i^2 T).
    _validate(lambdas, sigmas)
    density = 0.0
    for lam, sigma in zip(lambdas, sigmas):
        mu = np.log(S) + (r - 0.5 * sigma**2) * T
        s = sigma * np.sqrt(T)
        density = density + lam * norm.pdf((np.log(K) - mu) / s) / (K * s)
    return density
