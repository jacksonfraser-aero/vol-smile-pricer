import numpy as np
from scipy.stats import norm

def bs_call_price(S, K, T, r, sigma):
    d1 = (np.log(S / K) + (r + 0.5 * sigma**2) * T) / (sigma * np.sqrt(T))
    d2 = d1 - sigma * np.sqrt(T)
    return S * norm.cdf(d1) - K * np.exp(-r * T) * norm.cdf(d2)

def bs_vega(S, K, T, r, sigma):
    d1 = (np.log(S / K) + (r + 0.5 * sigma**2) * T) / (sigma * np.sqrt(T))
    return S * norm.pdf(d1) * np.sqrt(T)

SIGMA_LO, SIGMA_HI = 1e-4, 5.0

def implied_vol_newton(price, S, K, T, r, initial_guess=0.3, tol=1e-8, max_iter=100):
    # Safeguarded (rtsafe-style) Newton: keep sigma bracketed, fall back to bisection when the
    # Newton step leaves the bracket or vega is tiny. tol is on the sigma step, not the price.
    lo, hi = SIGMA_LO, SIGMA_HI
    # call price is increasing in sigma, so no root in the bracket means no valid IV
    # (this also rejects prices below intrinsic or above S)
    if not (bs_call_price(S, K, T, r, lo) < price < bs_call_price(S, K, T, r, hi)):
        return np.nan

    sigma = min(max(initial_guess, lo), hi)
    for _ in range(max_iter):
        diff = bs_call_price(S, K, T, r, sigma) - price
        vega = bs_vega(S, K, T, r, sigma)

        # shrink the bracket: model price too high means the true sigma is lower
        if diff > 0:
            hi = sigma
        else:
            lo = sigma

        newton_ok = vega > 1e-8
        if newton_ok:
            next_sigma = sigma - diff / vega
            newton_ok = lo < next_sigma < hi
        if not newton_ok:
            next_sigma = 0.5 * (lo + hi)

        step = next_sigma - sigma
        sigma = next_sigma
        if abs(step) < tol:
            break
    else:
        return np.nan  # bracket still wide: didn't converge

    if bs_vega(S, K, T, r, sigma) < 1e-6:  # price insensitive to sigma: answer isn't meaningful
        return np.nan
    return sigma