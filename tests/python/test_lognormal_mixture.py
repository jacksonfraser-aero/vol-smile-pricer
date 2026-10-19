import itertools

from implied_vol import bs_call_price
from lognormal_mixture import mixture_call_price

R = 0.05
STRIKES = [60, 80, 90, 100, 110, 120, 160]
EXPIRIES = [0.02, 0.1, 0.25, 0.5, 1.0, 2.0]
VOLS = [0.05, 0.1, 0.2, 0.4, 0.7, 1.0]


def test_single_component_mixture_equals_black_scholes():
    S = 100.0
    mismatches = []
    for K, T, sigma in itertools.product(STRIKES, EXPIRIES, VOLS):
        mix = mixture_call_price(S, K, T, R, [1.0], [sigma])
        bs = bs_call_price(S, K, T, R, sigma)
        if mix != bs:
            mismatches.append((K, T, sigma, mix, bs))
    assert mismatches == [], "\n".join(
        f"K={K} T={T} sigma={s}: mixture={m!r} bs={b!r}" for K, T, s, m, b in mismatches
    )
