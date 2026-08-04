import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("data/hedge_comparison.csv")

fig, ax = plt.subplots(figsize=(9, 5))
ax.plot(df["step"], df["flat_portfolio"], marker="o", label="Flat vol hedge")
ax.plot(df["step"], df["smile_portfolio"], marker="o", label="Smile-implied vol hedge")
ax.axhline(0, color="gray", linestyle="--", alpha=0.5)
ax.set_xlabel("Trading day")
ax.set_ylabel("Portfolio value ($)")
ax.set_title("Delta-Hedging P&L: Flat Vol vs. Smile-Implied Vol")
ax.legend()
plt.tight_layout()
plt.savefig("plots/hedge_pnl_comparison.png", dpi=150)
print("Saved plots/hedge_pnl_comparison.png")