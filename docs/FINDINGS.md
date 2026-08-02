# Empirical Test of the Constant-Volatility Assumption

## Motivation
My IB Extended Essay derived Black-Scholes and noted its central limitation:
the model assumes volatility is constant across strikes and expiries. This
is testable — if it were true, every option on the same underlying and
expiry should imply the same volatility when priced back through the model.

## Method
1. Pulled a live options chain for AAPL, expiry 2026-08-28 (~26 days out)
2. Filtered for liquidity (open interest > 5) and moneyness (strikes within
   70%-130% of spot) to exclude illiquid, unreliable quotes
3. For each remaining call strike, inverted the Black-Scholes formula using
   Newton-Raphson to solve for the volatility the market is actually
   pricing in, given the observed mid price
4. Plotted implied vol against strike

## Result
Spot was $308.91. Implied vol ranged from roughly 0.27 near the money
(around strikes 315-320) up to roughly 0.41 at the wings (strikes around
250 on the downside and 400 on the upside) — a spread of about 14 vol
points across a single underlying and expiry. The shape is a clean,
roughly symmetric smile with a clear minimum near spot.

This directly contradicts the flat-vol assumption: under true
Black-Scholes, this plot would be a horizontal line. Instead, every
strike implies a different volatility, which is exactly the flaw my
essay identified but didn't demonstrate empirically.

## Fit
Fit a quadratic-in-strike curve to the observed smile. [Fill in RMSE once
fit_smile.py runs.] The fit is used in Phase 4 to compute smile-implied
delta at strikes without direct quotes, as an input to the hedging
simulation.

## Why this matters
This is what Phase 4 tests directly: does hedging with the correct,
smile-implied volatility actually produce better P&L outcomes than
hedging with the flat assumption my essay's model relies on?