import itertools

import numpy as np
from scipy.integrate import quad

from implied_vol import bs_call_price
from lognormal_mixture import mixture_call_price, mixture_density, mixture_put_price

R = 0.05
STRIKES = [60, 80, 90, 100, 110, 120, 160]
EXPIRIES = [0.02, 0.1, 0.25, 0.5, 1.0, 2.0]
VOLS = [0.05, 0.1, 0.2, 0.4, 0.7, 1.0]
# (lambdas, sigmas): two- and three-component smiles, including a heavy low-weight tail
MIXTURES = [
    ([0.5, 0.5], [0.15, 0.35]),
    ([0.8, 0.2], [0.1, 0.6]),
    ([0.3, 0.5, 0.2], [0.1, 0.25, 0.6]),
    ([0.95, 0.05], [0.2, 1.5]),
]


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


def test_mixture_prices_satisfy_no_arbitrage_bounds_and_put_call_parity():
    S = 100.0
    strikes = range(40, 201, 5)
    violations = []
    for (lambdas, sigmas), K, T in itertools.product(MIXTURES, strikes, EXPIRIES):
        C = mixture_call_price(S, K, T, R, lambdas, sigmas)
        P = mixture_put_price(S, K, T, R, lambdas, sigmas)
        forward_intrinsic = S - K * np.exp(-R * T)
        call_ok = max(forward_intrinsic, 0.0) - 1e-10 <= C <= S + 1e-10
        put_ok = max(-forward_intrinsic, 0.0) - 1e-10 <= P <= K * np.exp(-R * T) + 1e-10
        parity_err = (C - P) - forward_intrinsic
        if not (call_ok and put_ok and abs(parity_err) <= 1e-10):
            violations.append((lambdas, sigmas, K, T, C, P, parity_err))
    assert violations == [], "\n".join(
        f"lambdas={l} sigmas={s} K={K} T={T}: C={c:.12g} P={p:.12g} parity_err={e:.3g}"
        for l, s, K, T, c, p, e in violations
    )


def test_density_is_non_negative_and_integrates_to_one():
    # Integrate over strikes 1e-6*S .. 1e4*S in x = ln K (dK = K dx); the component log-means
    # are passed as breakpoints so quad doesn't step over a narrow peak at short expiries.
    S = 100.0
    lo, hi = np.log(1e-6 * S), np.log(1e4 * S)
    grid = np.exp(np.linspace(lo, hi, 2001))
    failures = []
    for (lambdas, sigmas), T in itertools.product(MIXTURES, EXPIRIES):
        q = mixture_density(S, grid, T, R, lambdas, sigmas)
        log_means = [np.log(S) + (R - 0.5 * s**2) * T for s in sigmas]
        mass, _ = quad(
            lambda x: mixture_density(S, np.exp(x), T, R, lambdas, sigmas) * np.exp(x),
            lo, hi, points=log_means, limit=500, epsabs=1e-12, epsrel=1e-12,
        )
        if np.min(q) < 0 or abs(mass - 1.0) > 1e-6:
            failures.append((lambdas, sigmas, T, np.min(q), mass))
    assert failures == [], "\n".join(
        f"lambdas={l} sigmas={s} T={T}: min density={m:.3g} mass={a:.12g}"
        for l, s, T, m, a in failures
    )
