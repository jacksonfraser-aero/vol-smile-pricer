import numpy as np
from scipy.stats import norm

def bs_call_price(S, K, T, r, sigma):
    d1 = (np.log(S / K) + (r + 0.5 * sigma**2) * T) / (sigma * np.sqrt(T))
    d2 = d1 - sigma * np.sqrt(T)
    return S * norm.cdf(d1) - K * np.exp(-r * T) * norm.cdf(d2)

def bs_vega(S, K, T, r, sigma):
    d1 = (np.log(S / K) + (r + 0.5 * sigma**2) * T) / (sigma * np.sqrt(T))
    return S * norm.pdf(d1) * np.sqrt(T)

def implied_vol_newton(price, S, K, T, r, initial_guess=0.3, tol=1e-6, max_iter=100):
    sigma = initial_guess
    for _ in range(max_iter):
        model_price = bs_call_price(S, K, T, r, sigma)
        vega = bs_vega(S, K, T, r, sigma)

        if vega < 1e-8:  # avoid division blowup near-zero vega (deep ITM/OTM)
            return np.nan

        diff = model_price - price
        if abs(diff) < tol:
            return sigma

        sigma -= diff / vega

        if sigma <= 0:  # Newton-Raphson can overshoot into invalid territory
            sigma = 0.01

    return np.nan  # didn't converge