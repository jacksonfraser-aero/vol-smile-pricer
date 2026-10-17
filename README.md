# Vol Smile Pricer

![tests](https://github.com/jacksonfraser-aero/vol-smile-pricer/actions/workflows/tests.yml/badge.svg)

A pricing project that picks up where my IB EE left off. While deriving Black-Scholes in the EE, I emphasised its biggest limitation: it says volatility is fixed, but it's not in the real world. Here, I take that key flaw and turn it into something visible that can be actively traded, rather than just a known weakness.

What we were supposed to do: start with the standard pricer, then show real options data and why the constant-vol assumption breaks down (the "smile"), then fit a model to it and then test what it would actually cost to P&L if you ignore it and hedge using the wrong assumption.

## Structure

1. Black-Scholes closed-form pricer, binomial tree, and the Greeks - the baseline model
2. Pull an options chain and back out implied vol per strike using Newton-Raphson
3. Curve fitting to the smile (quadratic-in-strike)
4. Run a delta-hedge using flat vol or smile-implied vol and look at the P&L
5. (if time) plug the smile-implied vol into my other project: [kelly-betting-simulator](https://github.com/jacksonfraser-aero/kelly-betting-simulator)

## Stack

C++ for the pricer, binomial tree, Greeks, and the hedging simulation. Python for pulling data, solving for implied vol, and fitting the smile. Research in Python, engine in C++.

## Results

### Volatility smile

AAPL options chain, expiry 2026-08-28 (~26 days to expiry), spot $308.91. Filtered down to roughly 31 liquid call strikes (open interest > 5, within 70%-130% of spot). For each strike, inverted the Black-Scholes call using Newton-Raphson to find market-implied volatility.

IV was between about 0.27 at the near-the-money strikes (~315-320) up to about 0.41 at the wings (strikes 250 and 400) - a difference of about 14 vol points on a single underlying and expiry. True BS would make this a flat line; in reality, each strike prices in a different vol, directly contradicting the assumptions on which the derivation in my essay is based.

Fit a quadratic-in-strike curve to the smile:

$$IV(K) = 2.267419 \times 10^{-5} \cdot K^2 - 1.468070 \times 10^{-2} \cdot K + 2.6491$$

RMSE of 0.0139 (1.4 vol pts) on the ~31 market-implied points.

### P&L of the delta-hedge: flat vol and smile-implied vol

Simulated hedge of a short AAPL call, K=350, on a 19-trading-day historical price path to expiry. The two hedges sell the option at the same market price (the smile-implied vol at K=350). They only differ in the vol they use for the delta calculation:

- **Flat vol (0.27):** the constant volatility that my essay's model assumes
- **Smile-implied vol (0.288):** the implied volatility that my fitted smile produces at K=350

| | Final P&L |
|---|---|
| Flat vol hedge | $0.497 |
| Smile-implied vol hedge | $0.508 |

My initial take on this comparison had the smile hedge ending up $0.24 ahead, but it was an apples-to-oranges comparison: each hedger sold the option at its own vol, so the smile hedger received about $0.23 more premium at the start. That premium difference amounted to roughly 95% of the gap. At the same selling price, the gap on this path is $0.012.

One path won't tell you much, so I repeated this comparison over 2,000 simulated GBM paths (19 rebalances, same random paths for both hedges), once with the smile vol as the true realised vol, once with the flat vol as the truth:

| Realised vol | Smile less flat, average P&L (SE) | P&L std dev ratio (smile / flat) |
|---|---|---|
| 0.288 (smile is right) | -0.0071 (0.0054) | 0.961 |
| 0.270 (flat is right) | -0.0063 (0.0050) | 1.010 |

The hedge vol makes no measurable difference to the average P&L: both differences are within about 1.3 standard errors of zero. Hedging at the true vol does reduce the dispersion of P&L, as theory expects, but only by 1 to 4%, since with 19 rebalances the noise from discrete hedging (P&L standard deviation of ~1.0) is much larger than the effect of a 1.8 vol point difference. In practice, P&L depends on the premium relative to realised vol: when the option is written at 0.288 and the stock realises 0.27, both hedges make about $0.20 on average.

So at daily rebalancing, the smile matters for what you sell the option for, not how you hedge it.

Fitted smile at `plots/vol_smile_fitted.png`. Full hedging path at `plots/hedge_pnl_comparison.png`. Monte Carlo results at `data/mc_hedge_comparison.csv` (run `build/mc_hedge_comparison`).

## Testing

There is a Catch2 suite (`tests/`, run with `cmake --build build && ctest --test-dir build --output-on-failure`) on the C++ side that tests invariants and not hand-picked values:

- **Convergence rate:** binomial European prices converge to Black-Scholes at first order (error × N is bounded from N=50 to 4000)
- **Put-call parity:** $C - P = S - Ke^{-rT}$ for Black-Scholes and binomial prices
- **No-arbitrage bounds:** calls and puts lie within the model-free bounds for an entire grid of strikes, expiries and vols
- **Greeks vs finite differences:** delta, gamma, vega, theta and rho correspond to the central finite differences of the price at ITM, ATM and OTM strikes
- **Self-financing hedge:** at each rebalance, cash changes only by the interest earned plus the cost of the delta trade
- **P&L vs rebalancing frequency:** with realised vol equal to pricing vol, the P&L dispersion from hedging gets smaller with more frequent rebalancing

That self-financing test caught a real bug: the hedger was being rebalanced before interest was earned, so every trade earned interest for the interval before it actually happened. Moving the order of events brought the final P&L numbers down by around half a cent.

A separate test checks that premium vol and hedge vol can be set independently, which guards the fix for the premium contamination in my original flat vs smile comparison (explained in the hedging results above).

## Status

All phases 1-4 finished, plus a Monte Carlo check of the hedging comparison. Kelly project tie-in to be done if time available in phase 5.