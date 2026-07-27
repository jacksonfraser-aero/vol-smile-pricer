# Vol Smile Pricer

A pricing project that starts where my IB Extended Essay left off. The essay derived Black-Scholes but also pointed out its own weak spot — it assumes volatility is constant, which isn't true in real markets. This project builds that gap into something you can actually see and trade around, instead of just stating it as a limitation.

Plan is: build the standard pricer first, then pull real options data and show the constant-vol assumption breaking down (the "smile"), fit a model to it, and measure what it actually costs in P&L if you ignore it and hedge with the wrong assumption instead.

## Roughly how it's structured

1. Black-Scholes closed-form pricer, binomial tree, and the Greeks — the baseline model
2. Pull an options chain and back out implied vol per strike using Newton-Raphson — this is where the flat-vol assumption gets shown to be wrong
3. Fit a curve to the smile (quadratic or SVI)
4. Simulate delta-hedging with flat vol vs. smile-implied vol and compare the P&L
5. (if time) feed the smile-implied vol into my other project, [kelly-betting-simulator](link), which currently prices market-making spreads off a flat assumption

## Stack

C++ for the pricer, binomial tree, Greeks and the hedging simulation. Python for pulling data and fitting the smile. Basically research in Python, engine in C++.

## Status

Just getting started — pricer skeleton and test setup are in, nothing computes real numbers yet.