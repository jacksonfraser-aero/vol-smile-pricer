import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from implied_vol import implied_vol_newton

RISK_FREE_RATE = 0.05  # rough approximation; refine later if time allows

df = pd.read_csv("data/options_chain.csv")
calls = df[df["type"] == "call"].copy()

spot = calls["spot"].iloc[0]
expiry_date = datetime.strptime(calls["expiry"].iloc[0], "%Y-%m-%d")
T = (expiry_date - datetime.now()).days / 365.0

calls = calls[(calls["bid"] > 0) & (calls["ask"] > 0)]
calls = calls[
    (calls["volume"] > 0) &
    (calls["openInterest"] > 5) &
    (calls["strike"] > spot * 0.7) &
    (calls["strike"] < spot * 1.3)
]
calls["mid"] = (calls["bid"] + calls["ask"]) / 2

calls["my_iv"] = calls.apply(
    lambda row: implied_vol_newton(row["mid"], spot, row["strike"], T, RISK_FREE_RATE),
    axis=1
)

calls = calls.dropna(subset=["my_iv"])
calls.to_csv("data/implied_vols.csv", index=False)

fig, ax = plt.subplots(figsize=(9, 5))
ax.scatter(calls["strike"], calls["my_iv"], label="My Newton-Raphson IV")
ax.axvline(spot, color="gray", linestyle="--", alpha=0.5, label=f"Spot (${spot:.2f})")
ax.set_xlabel("Strike")
ax.set_ylabel("Implied Volatility")
ax.set_title(f"Volatility Smile: {calls['underlying'].iloc[0]}, expiry {calls['expiry'].iloc[0]}")
ax.legend()
plt.tight_layout()
plt.savefig("plots/vol_smile_raw.png", dpi=150)
print("Saved plots/vol_smile_raw.png")