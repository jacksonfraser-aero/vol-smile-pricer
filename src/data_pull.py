import yfinance as yf
import pandas as pd

TICKER = "AAPL"  # SPY is the other obvious choice for an SPX proxy

def pull_chain(ticker: str, expiry_index: int = 5):
    tk = yf.Ticker(ticker)
    expiries = tk.options
    print(f"Available expiries: {expiries[:5]}...")

    expiry = expiries[expiry_index]  # skip 0-DTE/near-term noise, pick something a few weeks out
    chain = tk.option_chain(expiry)

    calls = chain.calls.copy()
    calls["type"] = "call"
    puts = chain.puts.copy()
    puts["type"] = "put"

    df = pd.concat([calls, puts], ignore_index=True)
    df["expiry"] = expiry
    df["underlying"] = ticker

    spot = tk.history(period="1d")["Close"].iloc[-1]
    df["spot"] = spot

    return df, expiry, spot

if __name__ == "__main__":
    df, expiry, spot = pull_chain(TICKER)
    df.to_csv("data/options_chain.csv", index=False)
    print(f"Spot: {spot:.2f}, Expiry: {expiry}")
    print(f"Pulled {len(df)} contracts, saved to data/options_chain.csv")