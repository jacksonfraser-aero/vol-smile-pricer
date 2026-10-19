"""Calibrate Brigo-Mercurio lognormal mixtures to the market smile and compare with the quadratic.

Fits N=2 (lambda, sigma1, sigma2), N=3 (5 free parameters) and shifted N=2 (alpha too, 4 free
parameters) mixtures by least squares on implied vols, compares them with the bid-ask IV band,
then checks every fitted smile for butterfly and call-spread arbitrage via the Breeden-Litzenberger
density exp(rT) * d2C/dK2 on a fine strike grid.

The shifted fit is reported twice: with alpha >= 0 (S_T >= alpha*S*exp(rT) >= 0, so the price
stays positive) and with alpha free. Free, the fit improves monotonically as alpha -> -inf, where
the shifted lognormal mixture tends to a normal (Bachelier) mixture; the alpha profile is saved to
data/shift_profile.csv.
"""
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from scipy.stats import norm

from implied_vol import bs_call_price, implied_vol_newton
from lognormal_mixture import mixture_density, shifted_mixture_call_price

# S, r and T exactly as src/smile_plot.py builds them (that script writes data/implied_vols.csv
# and isn't importable without re-running it, so the constants are repeated here).
RISK_FREE_RATE = 0.05
VALUATION_DATE = datetime(2026, 8, 2, 17, 51)

N_STARTS = 25
NAN_PENALTY = 1.0  # residual (in vol) when the mixture price has no usable implied vol
NEGATIVE_DENSITY_FLOOR = -1e-9  # below finite-difference roundoff (~eps * C / dK^2 ~ 1e-11)

chain = pd.read_csv("data/options_chain.csv")
chain_calls = chain[chain["type"] == "call"]
S = float(chain_calls["spot"].iloc[0])
expiry_date = datetime.strptime(chain_calls["expiry"].iloc[0], "%Y-%m-%d")
T = (expiry_date - VALUATION_DATE).days / 365.0
r = RISK_FREE_RATE

df = pd.read_csv("data/implied_vols.csv")
strikes = df["strike"].values
market_ivs = df["my_iv"].values
assert len(df) == 31 and np.allclose(df["spot"], S), "implied_vols.csv doesn't match the chain"
print(f"S = {S!r}, r = {r}, T = {T!r} ({(expiry_date - VALUATION_DATE).days}/365), rows = {len(df)}")


# --- parametrisations: unconstrained x -> (alpha, lambdas, sigmas) ---------------------------

# Shift bound: every strike used (down to the 0.5*S end of the arbitrage grid) must exceed
# alpha*S*exp(rT), so alpha < 0.5*exp(-rT).
ALPHA_MAX = 0.5 * np.exp(-r * T)


def unpack(x, n, shift):
    # first n-1 entries are weight logits (last component's logit pinned at 0), then n log-sigmas,
    # then, if shifted, y with alpha = ALPHA_MAX - exp(y) (shift="free", alpha < ALPHA_MAX) or
    # alpha = ALPHA_MAX * logistic(y) (shift="positive", 0 <= alpha < ALPHA_MAX). Clip logits so
    # no weight underflows to exactly 0 (validation requires lambda > 0).
    logits = np.append(np.clip(x[: n - 1], -30, 30), 0.0)
    weights = np.exp(logits - logits.max())
    lambdas = weights / weights.sum()  # for n=2 this is the logistic of x[0]
    sigmas = np.exp(x[n - 1 : 2 * n - 1])
    if shift == "free":
        alpha = ALPHA_MAX - np.exp(x[2 * n - 1])
    elif shift == "positive":
        alpha = ALPHA_MAX / (1 + np.exp(-x[2 * n - 1]))
    else:
        alpha = 0.0
    return alpha, lambdas, sigmas


def mixture_call(K, alpha, lambdas, sigmas):
    # alpha = 0 is exactly the unshifted mixture (tested in test_lognormal_mixture.py)
    return shifted_mixture_call_price(S, K, T, r, alpha, lambdas, sigmas)


def mixture_ivs(K, alpha, lambdas, sigmas, guesses):
    return np.array([
        implied_vol_newton(mixture_call(k, alpha, lambdas, sigmas), S, k, T, r, initial_guess=g)
        for k, g in zip(K, guesses)
    ])


def prob_negative(alpha, lambdas, sigmas):
    # S_T = floor + X_T with floor = alpha*S*exp(rT); for alpha < 0, S_T < 0 exactly when X_T < -floor
    floor = alpha * S * np.exp(r * T)
    if floor >= 0:
        return 0.0
    x0 = (1 - alpha) * S
    return sum(
        lam * norm.cdf((np.log(-floor / x0) - (r - 0.5 * sig**2) * T) / (sig * np.sqrt(T)))
        for lam, sig in zip(lambdas, sigmas)
    )


def residuals(x, n, shift):
    alpha, lambdas, sigmas = unpack(x, n, shift)
    res = mixture_ivs(strikes, alpha, lambdas, sigmas, market_ivs) - market_ivs
    return np.where(np.isnan(res), NAN_PENALTY, res)


def fit_mixture(n, rng, shift=None):
    fits = []
    for _ in range(N_STARTS):
        x0 = np.concatenate([rng.uniform(-3, 3, n - 1), np.log(rng.uniform(0.05, 1.0, n))])
        if shift == "free":
            x0 = np.append(x0, np.log(ALPHA_MAX - rng.uniform(-0.5, 0.45)))
        elif shift == "positive":
            x0 = np.append(x0, rng.uniform(-3, 3))
        fits.append(least_squares(residuals, x0, args=(n, shift), x_scale=1.0))
    best = min(fits, key=lambda f: f.cost)
    costs = np.array([f.cost for f in fits])
    n_at_best = np.sum(costs <= best.cost * (1 + 1e-6) + 1e-14)
    label = f"shifted ({shift}) N={n}" if shift else f"N={n}"
    print(f"{label}: {N_STARTS} starts, {n_at_best} reached the best cost {best.cost:.6e}")
    alpha, lambdas, sigmas = unpack(best.x, n, shift)
    order = np.argsort(sigmas)
    return alpha, lambdas[order], sigmas[order]


rng = np.random.default_rng(0)
quad_coeffs = np.polyfit(strikes, market_ivs, deg=2)  # same quadratic-in-strike fit as fit_smile.py
fits = {
    "quadratic": dict(n_params=3, iv=lambda K: np.polyval(quad_coeffs, K),
                      params=f"a={quad_coeffs[0]:.6e}, b={quad_coeffs[1]:.6e}, c={quad_coeffs[2]:.6f}"),
}
SHIFTED_POSITIVE = "shifted N=2, alpha>=0"
SHIFTED_FREE = "shifted N=2, free (->normal)"
# the free shift is fitted before the positive one so the random starts of the earlier fits
# are unchanged
for name, n, shift in [("mixture N=2", 2, None), ("mixture N=3", 3, None),
                       (SHIFTED_FREE, 2, "free"), (SHIFTED_POSITIVE, 2, "positive")]:
    alpha, lambdas, sigmas = fit_mixture(n, rng, shift)
    fits[name] = dict(
        n_params=2 * n - 1 + (shift is not None), alpha=alpha, lambdas=lambdas, sigmas=sigmas,
        iv=lambda K, a=alpha, l=lambdas, s=sigmas: mixture_ivs(
            np.atleast_1d(K), a, l, s, np.full(np.size(K), 0.3)),
        params=(f"alpha={alpha:.6g}, " if shift else "")
               + "lambdas=[" + ", ".join(f"{v:.6f}" for v in lambdas)
               + "], sigmas=[" + ", ".join(f"{v:.6f}" for v in sigmas) + "]",
    )

fits = {k: fits[k] for k in ["quadratic", "mixture N=2", "mixture N=3", SHIFTED_POSITIVE, SHIFTED_FREE]}

# --- alpha profile of the free shift: fix alpha, refit (lambda, sigmas) ----------------------

PROFILE_ALPHAS = [0.0, -1.0, -5.0, -20.0, fits[SHIFTED_FREE]["alpha"], -1000.0]


def profile_residuals(x, alpha):
    lam = 1 / (1 + np.exp(-np.clip(x[0], -30, 30)))
    res = mixture_ivs(strikes, alpha, np.array([lam, 1 - lam]), np.exp(x[1:]), market_ivs) - market_ivs
    return np.where(np.isnan(res), NAN_PENALTY, res)


# Two warm starts per alpha, each keeping a fitted solution's per-component dollar vol
# sigma_i * X0: the free-shift fit and the unshifted N=2 fit. Near alpha = 0 the free-shift start
# alone falls into a slightly worse local minimum, so keep whichever refit is better.
warm_starts = [fits[SHIFTED_FREE], fits["mixture N=2"]]
profile = []
for alpha in PROFILE_ALPHAS:
    refits = []
    for start in warm_starts:
        abs_vols = start["sigmas"] * (1 - start["alpha"]) * S
        x0 = np.r_[np.log(start["lambdas"][0] / start["lambdas"][1]), np.log(abs_vols / ((1 - alpha) * S))]
        refits.append(least_squares(profile_residuals, x0, args=(alpha,)))
    f = min(refits, key=lambda f: f.cost)
    lam = 1 / (1 + np.exp(-f.x[0]))
    profile.append(dict(alpha=alpha, rmse=np.sqrt(np.mean(f.fun**2)), lambda1=lam,
                        sigma1=np.exp(f.x[1]), sigma2=np.exp(f.x[2])))
profile = pd.DataFrame(profile)
profile.to_csv("data/shift_profile.csv", index=False)
print("\nFree-shift profile (alpha fixed, lambda and sigmas refitted):")
print(profile.to_string(index=False, float_format=lambda v: f"{v:.6g}"))
print("Saved data/shift_profile.csv")

# --- bid-ask band ------------------------------------------------------------------------

quotes = chain_calls.set_index("strike").loc[strikes, ["bid", "ask"]]
bids, asks = quotes["bid"].values, quotes["ask"].values
bid_ivs = np.array([implied_vol_newton(b, S, k, T, r) for b, k in zip(bids, strikes)])
ask_ivs = np.array([implied_vol_newton(a, S, k, T, r) for a, k in zip(asks, strikes)])
both = ~np.isnan(bid_ivs) & ~np.isnan(ask_ivs)
print(f"\nBid-ask IVs: {np.isnan(bid_ivs).sum()} bid and {np.isnan(ask_ivs).sum()} ask without an IV "
      f"(quote outside no-arbitrage bounds); median width {np.median(ask_ivs[both] - bid_ivs[both]):.4f} "
      f"vol points over the {both.sum()} strikes with both")


def inside_band(fitted_ivs):
    # Compared in price space so a bid below intrinsic (nan IV) still counts as a lower edge;
    # BS is increasing in sigma, so this is the IV band wherever both IVs exist.
    prices = bs_call_price(S, strikes, T, r, fitted_ivs)
    return (bids <= prices) & (prices <= asks)


# --- comparison table --------------------------------------------------------------------

print("\nFit quality (errors in vol points, 31 strikes)")
print(f"{'model':<30}{'params':>7}{'RMSE':>10}{'max |err|':>11}{'in bid-ask':>12}  fitted parameters")
for name, fit in fits.items():
    fit["fitted_ivs"] = fit["iv"](strikes)
    err = fit["fitted_ivs"] - market_ivs
    fit["in_band"] = inside_band(fit["fitted_ivs"])
    in_band = f"{fit['in_band'].sum()}/{len(strikes)}"
    print(f"{name:<30}{fit['n_params']:>7}{np.sqrt(np.mean(err**2)):>10.4f}"
          f"{np.max(np.abs(err)):>11.4f}{in_band:>12}  {fit['params']}")
free = fits[SHIFTED_FREE]
print(f"  '{SHIFTED_FREE}' is approaching the normal-mixture (Bachelier) limit: alpha is unbounded\n"
      f"  below and RMSE keeps falling as alpha -> -inf (profile above). P(S_T < 0) = "
      f"{prob_negative(free['alpha'], free['lambdas'], free['sigmas']):.3e}")

# --- arbitrage check: Breeden-Litzenberger density on a fine grid ---------------------------

grid = np.linspace(0.5 * S, 1.5 * S, 4001)
dK = grid[1] - grid[0]
inner = grid[1:-1]
print(f"\nArbitrage check: K in [{grid[0]:.2f}, {grid[-1]:.2f}], dK = {dK:.4f}, "
      f"negative means density < {NEGATIVE_DENSITY_FLOOR:g}")
for name, fit in fits.items():
    if name == "quadratic":
        grid_ivs = fit["iv"](grid)
        if np.any(grid_ivs <= 0):
            print(f"  {name}: fitted IV <= 0 for some strikes on the grid; BS price undefined there")
        calls = bs_call_price(S, grid, T, r, grid_ivs)
    else:
        calls = mixture_call(grid, fit["alpha"], fit["lambdas"], fit["sigmas"])
    density = np.exp(r * T) * (calls[2:] - 2 * calls[1:-1] + calls[:-2]) / dK**2
    fit["density"] = density
    mass = np.trapezoid(density, inner)

    negative = density < NEGATIVE_DENSITY_FLOOR
    if negative.any():
        # contiguous runs of negative density, reported as strike intervals
        edges = np.flatnonzero(np.diff(np.concatenate([[0], negative.astype(int), [0]])))
        runs = ", ".join(
            f"[{inner[a]:.2f}, {inner[b - 1]:.2f}] (min {density[a:b].min():.3e})"
            for a, b in zip(edges[::2], edges[1::2])
        )
        verdict = f"NEGATIVE on {runs}"
    else:
        verdict = "non-negative everywhere"
    print(f"  {name:<30} min density {density.min():.3e}; integral over grid {mass:.6f}; {verdict}")
    # Call-spread check: the implied CDF 1 + exp(rT) dC/dK must stay in [0, 1], i.e. calls can't
    # get dearer as the strike rises. A positive density can still break this (mass > 1).
    slope = np.exp(r * T) * np.diff(calls) / dK
    rising = np.flatnonzero(slope > 0)
    if rising.size:
        print(f"  {'':<30} CALL SPREAD ARBITRAGE: call price rises with strike for K >= "
              f"{grid[rising[0]]:.2f} (max implied CDF {1 + slope.max():.6f})")
    else:
        print(f"  {'':<30} call price non-increasing in strike (implied CDF in [0, 1])")
    if name != "quadratic":
        # density of S_T = floor + X_T is the mixture density of X_T at K - floor
        floor = fit["alpha"] * S * np.exp(r * T)
        exact = mixture_density((1 - fit["alpha"]) * S, inner - floor, T, r, fit["lambdas"], fit["sigmas"])
        print(f"  {'':<30} (max |FD - analytic mixture density| = {np.max(np.abs(density - exact)):.2e})")
        if floor < 0:
            p_negative = prob_negative(fit["alpha"], fit["lambdas"], fit["sigmas"])
            print(f"  {'':<30} negative shift: S_T >= {floor:.2f}, P(S_T < 0) = {p_negative:.3e}")

# --- outputs --------------------------------------------------------------------------------

pd.DataFrame({
    "strike": strikes,
    "market_iv": market_ivs,
    "bid_iv": bid_ivs,
    "ask_iv": ask_ivs,
    "quadratic_iv": fits["quadratic"]["fitted_ivs"],
    "mixture2_iv": fits["mixture N=2"]["fitted_ivs"],
    "mixture3_iv": fits["mixture N=3"]["fitted_ivs"],
    "shifted2_positive_iv": fits[SHIFTED_POSITIVE]["fitted_ivs"],
    "shifted2_free_iv": fits[SHIFTED_FREE]["fitted_ivs"],
    "quadratic_in_band": fits["quadratic"]["in_band"],
    "mixture2_in_band": fits["mixture N=2"]["in_band"],
    "mixture3_in_band": fits["mixture N=3"]["in_band"],
    "shifted2_positive_in_band": fits[SHIFTED_POSITIVE]["in_band"],
    "shifted2_free_in_band": fits[SHIFTED_FREE]["in_band"],
}).to_csv("data/mixture_fit.csv", index=False)
print("\nSaved data/mixture_fit.csv")

STYLES = {  # categorical slots 1-5 of the reference palette, plus line style as secondary encoding
    "quadratic": dict(color="#2a78d6", linestyle="-"),
    "mixture N=2": dict(color="#eb6834", linestyle="--"),
    "mixture N=3": dict(color="#1baf7a", linestyle=":"),
    SHIFTED_POSITIVE: dict(color="#eda100", linestyle="-."),
    SHIFTED_FREE: dict(color="#e87ba4", linestyle=(0, (5, 1, 1, 1, 1, 1))),
}
fig, (ax_iv, ax_q) = plt.subplots(2, 1, figsize=(9, 9), sharex=True)
curve_K = np.linspace(strikes.min(), strikes.max(), 300)
ax_iv.scatter(strikes, market_ivs, color="#444444", s=24, zorder=3, label="Market IV")
for name, fit in fits.items():
    ax_iv.plot(curve_K, fit["iv"](curve_K), linewidth=2, label=name, **STYLES[name])
ax_iv.set_ylabel("Implied volatility")
ax_iv.set_title(f"Smile fits: {df['underlying'].iloc[0]}, expiry {df['expiry'].iloc[0]}")
ax_iv.legend(frameon=False)

for name, fit in fits.items():
    ax_q.plot(inner, fit["density"], linewidth=2, label=name, **STYLES[name])
ax_q.axhline(0, color="#999999", linewidth=0.8)
ax_q.set_xlabel("Strike")
ax_q.set_ylabel(r"Implied density  $e^{rT}\,\partial^2 C / \partial K^2$")
ax_q.legend(frameon=False)
for ax in (ax_iv, ax_q):
    ax.axvline(S, color="#999999", linestyle="--", linewidth=0.8)
    ax.grid(alpha=0.25)
    ax.spines[["top", "right"]].set_visible(False)
plt.tight_layout()
plt.savefig("plots/mixture_vs_quadratic.png", dpi=150)
print("Saved plots/mixture_vs_quadratic.png")
