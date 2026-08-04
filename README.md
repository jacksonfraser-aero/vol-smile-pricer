# Vol Smile Pricer

A pricing project that starts where my IB Extended Essay left off. The essay derived Black-Scholes but also pointed out its own weak spot — it assumes volatility is constant, which isn't true in real markets. This project builds that gap into something you can actually see and trade around, instead of just stating it as a limitation.

Plan was: build the standard pricer first, then pull real options data and show the constant-vol assumption breaking down (the "smile"), fit a model to it, and measure what it actually costs in P&L if you ignore it and hedge with the wrong assumption instead.

## Structure

1. Black-Scholes closed-form pricer, binomial tree, and the Greeks — the baseline model
2. Pull an options chain and back out implied vol per strike using Newton-Raphson
3. Fit a curve to the smile (quadratic-in-strike)
4. Simulate delta-hedging with flat vol vs. smile-implied vol and compare the P&L
5. (if time) feed the smile-implied vol into my other project, [kelly-betting-simulator](link)

## Stack

C++ for the pricer, binomial tree, Greeks, and the hedging simulation. Python for pulling data, solving for implied vol, and fitting the smile. Research in Python, engine in C++.

## Results

### Volatility smile
Pulled the AAPL options chain, expiry 2026-08-28 (~26 days out), spot at $308.91. Filtered to ~31 liquid call strikes (open interest > 5, moneyness within 70%-130% of spot). For each strike, inverted the Black-Scholes formula with Newton-Raphson to solve for market-implied volatility.

IV ranged from roughly 0.27 near the money (strikes around 315-320) up to ~0.41 at the wings (strikes 250 and 400) — a spread of about 14 vol points on a single underlying and expiry. Under true Black-Scholes this would be a flat line; instead every strike prices in a different volatility, directly contradicting the assumption my essay's derivation relies on.

Fit a quadratic-in-strike curve to the smile:

$$IV(K) = 2.267419 \times 10^{-5} \cdot K^2 - 1.468070 \times 10^{-2} \cdot K + 2.6491$$

RMSE of 0.0139 (1.4 vol points) against the ~31 market-implied points.

### Delta-hedging P&L: flat vol vs. smile-implied vol
Simulated hedging a short AAPL call, K=350, over a 19-trading-day historical price path leading into expiry. Compared two hedging strategies:

- **Flat vol (0.27):** the constant-volatility assumption my essay's model relies on
- **Smile-implied vol (0.288):** the volatility my fitted smile implies specifically at K=350

| | Final P&L |
|---|---|
| Flat vol hedge | $0.269 |
| Smile-implied vol hedge | $0.514 |

K=350 sits well away from the smile's minimum, so the flat 0.27 assumption underestimates the market-implied volatility there by roughly 1.8 vol points. Hedging with the smile-implied delta tracked the option's actual sensitivity more accurately across the path — the smile-vol hedge outperformed the flat-vol hedge by roughly $0.25 by expiry, on a position sized around a single contract. Scaled to institutional size, that gap is the real cost of the constant-vol assumption my essay identified.

See `plots/vol_smile_fitted.png` for the fitted smile and `plots/hedge_pnl_comparison.png` for the full hedging path.

## Status
Phases 1-4 complete. Phase 5 (Kelly project tie-in) optional, time permitting.