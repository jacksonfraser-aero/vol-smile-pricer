import itertools

import numpy as np

from implied_vol import bs_call_price

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
