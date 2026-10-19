"""Calibrate Brigo-Mercurio lognormal mixtures to the market smile and compare with the quadratic.

Fits N=2 (lambda, sigma1, sigma2) and N=3 (5 free parameters) mixtures by least squares on
implied vols, then checks every fitted smile for butterfly arbitrage via the Breeden-Litzenberger
density exp(rT) * d2C/dK2 on a fine strike grid.
"""
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from implied_vol import bs_call_price, implied_vol_newton
from lognormal_mixture import mixture_call_price, mixture_density

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


# --- parametrisations: unconstrained x -> (lambdas, sigmas) ---------------------------------

def unpack(x, n):
    # first n-1 entries are weight logits (last component's logit pinned at 0), rest log-sigmas.
    # Clip logits so no weight underflows to exactly 0 (validation requires lambda > 0).
    logits = np.append(np.clip(x[: n - 1], -30, 30), 0.0)
    weights = np.exp(logits - logits.max())
    lambdas = weights / weights.sum()  # for n=2 this is the logistic of x[0]
    sigmas = np.exp(x[n - 1 :])
    return lambdas, sigmas


def mixture_ivs(K, lambdas, sigmas, guesses):
    return np.array([
        implied_vol_newton(mixture_call_price(S, k, T, r, lambdas, sigmas), S, k, T, r, initial_guess=g)
        for k, g in zip(K, guesses)
    ])


def residuals(x, n):
    lambdas, sigmas = unpack(x, n)
    res = mixture_ivs(strikes, lambdas, sigmas, market_ivs) - market_ivs
    return np.where(np.isnan(res), NAN_PENALTY, res)


def fit_mixture(n, rng):
    fits = []
    for _ in range(N_STARTS):
        x0 = np.concatenate([rng.uniform(-3, 3, n - 1), np.log(rng.uniform(0.05, 1.0, n))])
        fits.append(least_squares(residuals, x0, args=(n,), x_scale=1.0))
    best = min(fits, key=lambda f: f.cost)
    costs = np.array([f.cost for f in fits])
    n_at_best = np.sum(costs <= best.cost * (1 + 1e-6) + 1e-14)
    print(f"N={n}: {N_STARTS} starts, {n_at_best} reached the best cost {best.cost:.6e}")
    lambdas, sigmas = unpack(best.x, n)
    order = np.argsort(sigmas)
    return lambdas[order], sigmas[order]


rng = np.random.default_rng(0)
quad_coeffs = np.polyfit(strikes, market_ivs, deg=2)  # same quadratic-in-strike fit as fit_smile.py
fits = {
    "quadratic": dict(n_params=3, iv=lambda K: np.polyval(quad_coeffs, K),
                      params=f"a={quad_coeffs[0]:.6e}, b={quad_coeffs[1]:.6e}, c={quad_coeffs[2]:.6f}"),
}
for n in (2, 3):
    lambdas, sigmas = fit_mixture(n, rng)
    fits[f"mixture N={n}"] = dict(
        n_params=2 * n - 1, lambdas=lambdas, sigmas=sigmas,
        iv=lambda K, l=lambdas, s=sigmas: mixture_ivs(np.atleast_1d(K), l, s, np.full(np.size(K), 0.3)),
        params="lambdas=[" + ", ".join(f"{v:.6f}" for v in lambdas)
               + "], sigmas=[" + ", ".join(f"{v:.6f}" for v in sigmas) + "]",
    )

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
print(f"{'model':<14}{'params':>7}{'RMSE':>10}{'max |err|':>11}{'in bid-ask':>12}  fitted parameters")
for name, fit in fits.items():
    fit["fitted_ivs"] = fit["iv"](strikes)
    err = fit["fitted_ivs"] - market_ivs
    fit["in_band"] = inside_band(fit["fitted_ivs"])
    in_band = f"{fit['in_band'].sum()}/{len(strikes)}"
    print(f"{name:<14}{fit['n_params']:>7}{np.sqrt(np.mean(err**2)):>10.4f}"
          f"{np.max(np.abs(err)):>11.4f}{in_band:>12}  {fit['params']}")

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
        calls = mixture_call_price(S, grid, T, r, fit["lambdas"], fit["sigmas"])
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
    print(f"  {name:<14} min density {density.min():.3e}; integral over grid {mass:.6f}; {verdict}")
    # Call-spread check: the implied CDF 1 + exp(rT) dC/dK must stay in [0, 1], i.e. calls can't
    # get dearer as the strike rises. A positive density can still break this (mass > 1).
    slope = np.exp(r * T) * np.diff(calls) / dK
    rising = np.flatnonzero(slope > 0)
    if rising.size:
        print(f"  {'':<14} CALL SPREAD ARBITRAGE: call price rises with strike for K >= "
              f"{grid[rising[0]]:.2f} (max implied CDF {1 + slope.max():.6f})")
    else:
        print(f"  {'':<14} call price non-increasing in strike (implied CDF in [0, 1])")
    if name != "quadratic":
        exact = mixture_density(S, inner, T, r, fit["lambdas"], fit["sigmas"])
        print(f"  {'':<14} (max |FD - analytic mixture density| = {np.max(np.abs(density - exact)):.2e})")

# --- outputs --------------------------------------------------------------------------------

pd.DataFrame({
    "strike": strikes,
    "market_iv": market_ivs,
    "bid_iv": bid_ivs,
    "ask_iv": ask_ivs,
    "quadratic_iv": fits["quadratic"]["fitted_ivs"],
    "mixture2_iv": fits["mixture N=2"]["fitted_ivs"],
    "mixture3_iv": fits["mixture N=3"]["fitted_ivs"],
    "quadratic_in_band": fits["quadratic"]["in_band"],
    "mixture2_in_band": fits["mixture N=2"]["in_band"],
    "mixture3_in_band": fits["mixture N=3"]["in_band"],
}).to_csv("data/mixture_fit.csv", index=False)
print("\nSaved data/mixture_fit.csv")

STYLES = {  # categorical slots 1-3 of the reference palette, plus line style as secondary encoding
    "quadratic": dict(color="#2a78d6", linestyle="-"),
    "mixture N=2": dict(color="#eb6834", linestyle="--"),
    "mixture N=3": dict(color="#1baf7a", linestyle=":"),
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
