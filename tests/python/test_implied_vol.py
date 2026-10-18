import itertools

import numpy as np

from implied_vol import bs_call_price, bs_vega

R = 0.05
STRIKES = [60, 80, 90, 100, 110, 120, 160]
EXPIRIES = [0.02, 0.1, 0.25, 0.5, 1.0, 2.0]
VOLS = [0.05, 0.1, 0.2, 0.4, 0.7, 1.0]


def test_call_price_within_no_arbitrage_bounds():
    S = 100.0
    violations = []
    for K, T, sigma in itertools.product(STRIKES, EXPIRIES, VOLS):
        C = bs_call_price(S, K, T, R, sigma)
        lower = max(S - K * np.exp(-R * T), 0.0)
        if not (lower - 1e-10 <= C <= S + 1e-10):
            violations.append((K, T, sigma, lower, C))
    assert violations == [], "\n".join(
        f"K={K} T={T} sigma={s}: lower={lo:.12g} C={c:.12g} S={S}" for K, T, s, lo, c in violations
    )


def test_vega_matches_finite_difference_in_sigma():
    S, h = 100.0, 1e-5
    mismatches = []
    for K, T, sigma in itertools.product(STRIKES, EXPIRIES, VOLS):
        fd = (bs_call_price(S, K, T, R, sigma + h) - bs_call_price(S, K, T, R, sigma - h)) / (2 * h)
        vega = bs_vega(S, K, T, R, sigma)
        if not np.isclose(vega, fd, rtol=1e-5, atol=1e-6):
            mismatches.append((K, T, sigma, vega, fd))
    assert mismatches == [], "\n".join(
        f"K={K} T={T} sigma={s}: vega={v:.10g} fd={f:.10g}" for K, T, s, v, f in mismatches
    )
