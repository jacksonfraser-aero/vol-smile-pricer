import itertools

import numpy as np
import pytest
from scipy.integrate import quad

from implied_vol import bs_call_price, implied_vol_newton
from lognormal_mixture import (
    mixture_call_delta,
    mixture_call_price,
    mixture_density,
    mixture_put_price,
    shifted_mixture_call_price,
    shifted_mixture_put_price,
)

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


def test_second_strike_derivative_of_call_price_recovers_density():
    # Breeden-Litzenberger: d2C/dK2 = exp(-rT) q(K). Central difference with h = 0.01 has
    # truncation error ~ h^2/12 * q'''' (well under 1e-6 here) and roundoff ~ eps*S/h^2 ~ 1e-10.
    S, h = 100.0, 1e-2
    strikes = range(40, 201, 5)
    mismatches = []
    for (lambdas, sigmas), K, T in itertools.product(MIXTURES, strikes, EXPIRIES):
        c_up = mixture_call_price(S, K + h, T, R, lambdas, sigmas)
        c_mid = mixture_call_price(S, K, T, R, lambdas, sigmas)
        c_dn = mixture_call_price(S, K - h, T, R, lambdas, sigmas)
        fd = np.exp(R * T) * (c_up - 2 * c_mid + c_dn) / h**2
        q = mixture_density(S, K, T, R, lambdas, sigmas)
        if not np.isclose(fd, q, rtol=1e-4, atol=1e-7):
            mismatches.append((lambdas, sigmas, K, T, fd, q))
    assert mismatches == [], "\n".join(
        f"lambdas={l} sigmas={s} K={K} T={T}: fd={f:.10g} density={d:.10g}"
        for l, s, K, T, f, d in mismatches
    )


def test_mixture_implied_vol_lies_between_component_vols():
    # The mixture price is a convex combination of BS prices and BS is increasing in sigma, so
    # the implied vol must lie in [min sigma_i, max sigma_i]. nan (solver declined) is tolerated,
    # as in test_implied_vol; a finite answer outside the range never is. Slack = solver tol.
    S, slack = 100.0, 1e-8
    strikes = range(40, 201, 5)
    in_range_count, nan_count, outside = 0, 0, []
    for (lambdas, sigmas), K, T in itertools.product(MIXTURES, strikes, EXPIRIES):
        price = mixture_call_price(S, K, T, R, lambdas, sigmas)
        iv = implied_vol_newton(price, S, K, T, R)
        if np.isnan(iv):
            nan_count += 1
        elif min(sigmas) - slack <= iv <= max(sigmas) + slack:
            in_range_count += 1
        else:
            outside.append((lambdas, sigmas, K, T, price, iv))
    print(f"\nmixture IVs: {in_range_count} in range, {nan_count} nan, {len(outside)} outside")
    assert in_range_count > 0
    assert outside == [], "\n".join(
        f"lambdas={l} sigmas={s} K={K} T={T}: price={p:.10g} iv={v:.10g}"
        for l, s, K, T, p, v in outside
    )


def test_mixture_delta_matches_finite_difference_in_spot():
    S, h = 100.0, 1e-4
    strikes = range(40, 201, 5)
    mismatches = []
    for (lambdas, sigmas), K, T in itertools.product(MIXTURES[:3], strikes, EXPIRIES):
        c_up = mixture_call_price(S + h, K, T, R, lambdas, sigmas)
        c_dn = mixture_call_price(S - h, K, T, R, lambdas, sigmas)
        fd = (c_up - c_dn) / (2 * h)
        delta = mixture_call_delta(S, K, T, R, lambdas, sigmas)
        if not np.isclose(delta, fd, rtol=1e-6, atol=1e-8):
            mismatches.append((lambdas, sigmas, K, T, delta, fd))
    assert mismatches == [], "\n".join(
        f"lambdas={l} sigmas={s} K={K} T={T}: delta={d:.10g} fd={f:.10g}"
        for l, s, K, T, d, f in mismatches
    )


INVALID_MIXTURES = {
    "length mismatch": ([0.5, 0.5], [0.2]),
    "zero lambda": ([0.0, 1.0], [0.2, 0.3]),
    "negative lambda": ([-0.2, 1.2], [0.2, 0.3]),
    "lambdas sum above 1": ([0.5, 0.5 + 1e-9], [0.2, 0.3]),
    "lambdas sum below 1": ([0.4, 0.5], [0.2, 0.3]),
    "zero sigma": ([0.5, 0.5], [0.0, 0.3]),
    "negative sigma": ([0.5, 0.5], [0.2, -0.3]),
}
PUBLIC_FUNCTIONS = [mixture_call_price, mixture_put_price, mixture_call_delta, mixture_density]


@pytest.mark.parametrize("case", INVALID_MIXTURES)
@pytest.mark.parametrize("fn", PUBLIC_FUNCTIONS, ids=lambda f: f.__name__)
def test_invalid_mixture_raises_value_error(fn, case):
    lambdas, sigmas = INVALID_MIXTURES[case]
    with pytest.raises(ValueError):
        fn(100.0, 100.0, 0.5, R, lambdas, sigmas)


def test_zero_shift_reproduces_unshifted_mixture_exactly():
    S = 100.0
    mismatches = []
    for (lambdas, sigmas), K, T in itertools.product(MIXTURES, STRIKES, EXPIRIES):
        shifted = shifted_mixture_call_price(S, K, T, R, 0.0, lambdas, sigmas)
        plain = mixture_call_price(S, K, T, R, lambdas, sigmas)
        if shifted != plain:
            mismatches.append((lambdas, sigmas, K, T, shifted, plain))
    assert mismatches == [], "\n".join(
        f"lambdas={l} sigmas={s} K={K} T={T}: shifted={a!r} unshifted={b!r}"
        for l, s, K, T, a, b in mismatches
    )


def test_shifted_mixture_satisfies_no_arbitrage_bounds_and_put_call_parity():
    # Non-negative shifts only: S_T >= alpha*S*exp(rT), so the asset stays non-negative and the
    # put bound P <= K*exp(-rT) must hold. (A negative shift lets S_T go below zero.)
    S = 100.0
    strikes = range(40, 201, 5)  # all above alpha*S*exp(rT) <= 0.3*100*exp(0.1) = 33.2
    violations = []
    for (lambdas, sigmas), alpha, K, T in itertools.product(MIXTURES, [0.1, 0.3], strikes, EXPIRIES):
        C = shifted_mixture_call_price(S, K, T, R, alpha, lambdas, sigmas)
        P = shifted_mixture_put_price(S, K, T, R, alpha, lambdas, sigmas)
        forward_intrinsic = S - K * np.exp(-R * T)
        call_ok = max(forward_intrinsic, 0.0) - 1e-10 <= C <= S + 1e-10
        put_ok = max(-forward_intrinsic, 0.0) - 1e-10 <= P <= K * np.exp(-R * T) + 1e-10
        parity_err = (C - P) - forward_intrinsic
        if not (call_ok and put_ok and abs(parity_err) <= 1e-10):
            violations.append((lambdas, sigmas, alpha, K, T, C, P, parity_err))
    assert violations == [], "\n".join(
        f"lambdas={l} sigmas={s} alpha={a} K={K} T={T}: C={c:.12g} P={p:.12g} parity_err={e:.3g}"
        for l, s, a, K, T, c, p, e in violations
    )


@pytest.mark.parametrize("fn", [shifted_mixture_call_price, shifted_mixture_put_price],
                         ids=lambda f: f.__name__)
def test_shifted_mixture_rejects_strike_at_or_below_shifted_floor(fn):
    S, T, alpha = 100.0, 1.0, 0.3
    floor = alpha * S * np.exp(R * T)
    for K in [floor, floor - 1.0, np.array([50.0, floor - 1.0])]:
        with pytest.raises(ValueError):
            fn(S, K, T, R, alpha, [0.5, 0.5], [0.2, 0.3])
