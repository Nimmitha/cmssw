import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import math
import os

# Read the CSV
# df = pd.read_csv("output/csvs/combined_scan_count.csv")
df = pd.read_csv("output/csvs/combined_scan_count_fit.csv")

# Gaussian model for fitting
def gaussian(x, norm, mean, sigma):
    return norm * np.exp(-0.5 * ((x - mean) / sigma)**2)

def gaussian_fit_vdm(tag, subdf):
    print(f"\n=== Fitting Gaussian for {tag} ===")
    subdf = subdf.sort_values(by="sep")

    xdata = subdf["sep"].values
    ydata = subdf["rate"].values
    yerr = subdf["rate_err"].values
    # ydata = subdf["rate_adjusted"].values
    # yerr = subdf["rate_adjusted_err"].values

    # Initial guess: norm, mean, sigma
    norm_guess = max(ydata)
    mean_guess = 0.0
    sigma_guess = 0.01

    try:
        popt, pcov = curve_fit(gaussian, xdata, ydata, sigma=yerr, p0=[norm_guess, mean_guess, sigma_guess], absolute_sigma=True)
        norm, mean, sigma = popt
        norm_err, mean_err, sigma_err = np.sqrt(np.diag(pcov))
    except RuntimeError:
        print(f"Fit failed for {tag}")
        return None, None

    height_at_0 = gaussian(0, norm, mean, sigma)

    # Plotting
    fig, (ax_top, ax_bottom) = plt.subplots(2, 1, figsize=(8, 8), gridspec_kw={'height_ratios': [3, 1]}, sharex=True)

    # Top plot: Fit and data
    ax_top.errorbar(xdata, ydata, yerr=yerr, fmt='o', label="Data", capsize=2)
    xfit = np.linspace(min(xdata) - 0.01, max(xdata) + 0.01, 1000)
    yfit = gaussian(xfit, *popt)
    ax_top.plot(xfit, yfit, 'r-', label="Gaussian Fit")
    ax_top.set_yscale("log")
    ax_top.set_ylabel("Rate")
    ax_top.legend()

    # Stat box
    textstr = '\n'.join((
        rf'$\mu = {mean:.6f} ± {mean_err:.6f}$',
        rf'$\sigma = {sigma:.4f} ± {sigma_err:.4f}$',
        rf'$Norm = {norm:.1f} ± {norm_err:.1f}$',
        rf'$Height(0) = {height_at_0:.2f}$'
    ))
    ax_top.text(0.60, 0.75, textstr, transform=ax_top.transAxes, fontsize=10,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    # Bottom plot: Pulls
    residuals = ydata - gaussian(xdata, *popt)
    pulls = residuals / yerr
    ax_bottom.axhline(0, color='gray', linestyle='--')
    ax_bottom.plot(xdata, pulls, 'ko')
    ax_bottom.set_ylabel("Pull")
    ax_bottom.set_xlabel("Separation")

    plt.tight_layout()
    os.makedirs("output/plots", exist_ok=True)
    # plt.savefig(f"output/plots/gaussian_fit_{tag}_adjusted.png")
    plt.savefig(f"output/plots/gaussian_fit_{tag}.png")
    plt.close()

    print(f"\n=== {tag} Fit Results ===")
    print(f"Mean         : {mean:.6f} ± {mean_err:.6f}")
    print(f"Sigma        : {sigma:.6f} ± {sigma_err:.6f}")
    print(f"Norm         : {norm:.1f} ± {norm_err:.1f}")
    print(f"Height@0     : {height_at_0:.2f} counts")

    return sigma, norm, sigma_err, norm_err

# Do fits for X and Y scans
sigmas = {}
heights = {}
sigma_errors = {}
height_errors = {}
for tag in ["X1", "Y1"]:
    subdf = df[df["type"] == tag]
    sigma, height, sigma_err, height_err = gaussian_fit_vdm(tag, subdf)
    sigmas[tag] = sigma
    heights[tag] = height
    sigma_errors[tag] = sigma_err
    height_errors[tag] = height_err

# === Visible cross-section calculation ===
### 2023
# n1 = 1.54e11    # protons in beam 1 (10 fill average)
# n2 = 1.56e11    # protons in beam 2 (10 fill average)
# n1_error = 0.03e11
# n2_error = 0.03e11
# ### 2024 
n1 = 1.52e11
n2 = 1.52e11
n1_error = 0.01e11
n2_error = 0.01e11

f_rev = 11245  # LHC revolution frequency in Hz
N_bunch = 2452    # number of colliding bunches (my 10 fills)

# sigma_vis = (2 * math.pi * heights["X1"] * heights["Y1"]) / (N_bunch * f_rev * n1 * n2)
sigma_vis = (2*math.pi * sigmas["X1"] * sigmas["Y1"]) * (heights["X1"] + heights["Y1"]) / (2 * N_bunch * f_rev * n1 * n2)
svis_err = 0
svis_err += sigma_errors["X1"]**2 / (sigmas["X1"]**2)
svis_err += sigma_errors["Y1"]**2 / (sigmas["Y1"]**2)
svis_err += (height_errors["X1"]**2 + height_errors["Y1"]**2) / ((heights["X1"] + heights["Y1"])**2)
svis_err += n1_error**2 / n1**2
svis_err += n2_error**2 / n2**2

sigma_vis_error = sigma_vis * math.sqrt(svis_err)

print(f"\n=== Visible Cross Section ===")
print(f"σ_vis = {sigma_vis*1e28 * 1000:.4f} nb")  # Convert cm² to nb
print(f"σ_vis error = {sigma_vis_error*1e28 * 1000:.4f} nb")  # Convert cm² to nb
