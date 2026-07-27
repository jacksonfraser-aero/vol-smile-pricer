# Project Plan — Vol Smile Pricer

## Motivation
My IB Extended Essay derived Black-Scholes but also stated its central
limitation plainly: the model does not adjust for a volatility which is
not constant. This project builds that flaw into something empirical
and tradeable — first the standard pricer, then proof the constant-vol
assumption breaks in real markets, then a quantified cost of ignoring it.

## Phases

**Phase 1 — Core pricer**
Black-Scholes closed-form (call/put) + binomial tree + analytical
Greeks (delta, gamma, theta, vega, rho). Baseline model — the one
Phase 2 onward shows is wrong in practice.

**Phase 2 — Implied volatility from real data**
Pull an options chain (SPX or a large-cap name) from a free source.
For every strike/expiry, use Newton-Raphson to invert the BS formula
and solve for market-implied volatility. Plot it — expect a smile or
skew, not a flat line.

**Phase 3 — Smile fitting**
Fit a parametric curve to the smile (quadratic-in-strike, or SVI).
This is the interpolation/extrapolation step real vol desks do
between quoted strikes.

**Phase 4 — Delta-hedging P&L**
Simulate hedging a short option position over a historical price
path two ways: (a) BS delta under a flat, wrong volatility assumption,
(b) delta computed from the fitted smile. Compare hedging P&L —
turns "the model is imperfect" into a dollar figure.

**Phase 5 (optional, if time) — Kelly project tie-in**
Feed smile-implied vol into the kelly-betting-simulator market-making
spread pricing, replacing its flat-vol assumption.

## Stack
- C++: pricing engine, Greeks, binomial tree, hedging simulation harness
- Python: options data pull, Newton-Raphson IV solver, smile fitting
- Mirrors research-in-Python / production-in-C++ split used at real
  quant desks.

## Day-by-day commit plan

| Date | Session | Commit |
|---|---|---|
| Mon 27 (eve) | C++ setup | Repo init: skeleton, README linking to IB essay derivation, project plan |
| Tue 28 (eve) | C++ | BS closed-form pricer class (call/put) + basic unit tests |
| Wed 29 (eve) | C++ | Binomial tree pricer, parameterized by step count |
| Thu 30 (am) | C++ | Convergence test: binomial → BS as steps increase, with output/plot |
| Fri 31 (full day) | C++ | Analytical Greeks — corrected vs. essay's original delta/gamma statements |
| Sun 2 (full day) | Python | Options chain data pull (SPX or large-cap) + Newton-Raphson IV solver; raw smile plot |
| Thu 6 (full day) | Python | Smile-fitting (quadratic-in-strike or SVI) + writeup tying back to essay's limitations section |
| Fri 7 (full day) | C++ | Delta-hedging simulation harness: flat-vol BS delta |
| Sat 8 (am) | C++/Python | Second hedging run using smile-implied delta; P&L comparison chart; README results section |
| Wed 12 (full day, if time) | Python | Optional: feed smile-implied vol into Kelly project (Phase 5) |

## Status
Phase 1 in progress.