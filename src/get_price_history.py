import yfinance as yf

TICKER = "AAPL"

def fetch_history(ticker, period="2mo"):
    tk = yf.Ticker(ticker)
    hist = tk.history(period=period)
    hist = hist[["Close"]].reset_index()
    # keep only the last ~19 trading days, matching the ~26 calendar day option life
    hist = hist.tail(19)
    hist.to_csv("data/price_history.csv", index=False)
    print(f"Saved {len(hist)} rows to data/price_history.csv")

if __name__ == "__main__":
    fetch_history(TICKER)