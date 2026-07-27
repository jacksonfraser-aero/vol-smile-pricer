# Vol Smile Pricer

## Motivation
Links back to the IB essay — one sentence stating the essay's own
admission that Black-Scholes doesn't adjust for non-constant volatility.
Link the essay/PDF if you have it in the repo or link out to it.

## What this project does
Short paragraph: builds the standard pricer, then empirically shows
where it breaks (the smile), fits a model for it, and quantifies
the P&L cost of ignoring it via delta-hedging simulation.

## Phases
1. BS + binomial tree + Greeks (baseline)
2. Implied vol extraction from real option chains (Newton-Raphson)
3. Smile fitting (quadratic / SVI)
4. Delta-hedging P&L: flat-vol vs smile-implied
5. (optional) feed into kelly-betting-simulator market-making pricing

## Status
Phase 1 in progress.

## Stack
C++ (pricing/hedging engine) + Python (data pull, smile fitting) —
mirrors research-in-Python / production-in-C++ split at real desks.
