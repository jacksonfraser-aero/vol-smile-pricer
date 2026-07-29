import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("data/convergence.csv")

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(df["steps"], df["binomial_price"], marker="o", label="Binomial tree price")
ax.axhline(df["bs_price"].iloc[0], color="red", linestyle="--", label="Black-Scholes price")
ax.set_xscale("log")
ax.set_xlabel("Number of steps (log scale)")
ax.set_ylabel("Option price")
ax.set_title("Binomial Tree Convergence to Black-Scholes")
ax.legend()
plt.tight_layout()
plt.savefig("plots/convergence.png", dpi=150)
print("Saved plots/convergence.png")