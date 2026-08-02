import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv("data/implied_vols.csv")

strikes = df["strike"].values
ivs = df["my_iv"].values
spot = df["spot"].iloc[0]

# quadratic-in-strike fit: IV(K) = a*K^2 + b*K + c
coeffs = np.polyfit(strikes, ivs, deg=2)
a, b, c = coeffs

def fitted_iv(K):
    return a * K**2 + b * K + c

strike_range = np.linspace(strikes.min(), strikes.max(), 200)
fitted_curve = fitted_iv(strike_range)

fig, ax = plt.subplots(figsize=(9, 5))
ax.scatter(strikes, ivs, label="Market-implied IV (Newton-Raphson)", alpha=0.7)
ax.plot(strike_range, fitted_curve, color="red", label="Quadratic fit")
ax.axvline(spot, color="gray", linestyle="--", alpha=0.5, label=f"Spot (${spot:.2f})")
ax.set_xlabel("Strike")
ax.set_ylabel("Implied Volatility")
ax.set_title(f"Fitted Volatility Smile: {df['underlying'].iloc[0]}, expiry {df['expiry'].iloc[0]}")
ax.legend()
plt.tight_layout()
plt.savefig("plots/vol_smile_fitted.png", dpi=150)
print(f"Fit: IV(K) = {a:.6e}*K^2 + {b:.6e}*K + {c:.4f}")
print("Saved plots/vol_smile_fitted.png")

residuals = ivs - fitted_iv(strikes)
print(f"RMSE: {np.sqrt(np.mean(residuals**2)):.4f}")

fit_df = pd.DataFrame({"strike": strikes, "market_iv": ivs, "fitted_iv": fitted_iv(strikes)})
fit_df.to_csv("data/smile_fit.csv", index=False)