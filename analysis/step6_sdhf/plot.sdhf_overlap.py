#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.ndimage import gaussian_filter1d  # no longer used (Daniell smoothing below)

DANIELL_WIDTH_THZ = 0.5   # Daniell (boxcar) smoothing width

def daniell_smooth(y, w):
    """Daniell (moving-average) smoothing over w points (made odd)."""
    if w <= 1:
        return y.copy()
    if w % 2 == 0:
        w += 1
    k = np.ones(w) / w
    return np.convolve(y, k, mode="same")

def daniell_points(freq_hz, width_thz=DANIELL_WIDTH_THZ):
    df_thz = (freq_hz[1] - freq_hz[0]) / 1e12
    return int(round(width_thz / df_thz))
import os

press = np.array(["50bar", "100bar", "200bar", "250bar", "300bar", "350bar", "400bar"])
os.makedirs('plots/sdhf', exist_ok=True)

overlaps = np.zeros(len(press))

for i, p in enumerate(press):
    # Load fluid data
    dat_fluid = np.loadtxt(f"sdhf_consolidated/{p}_fluid.dat", skiprows=1)
    freq_hz = dat_fluid[:, 0]
    pssum_fluid = dat_fluid[:, 4]
    fTHz = freq_hz / 1e12

    # Load solid data
    dat_solid = np.loadtxt(f"sdhf_consolidated/{p}_solid.dat", skiprows=1)
    pssum_solid = -dat_solid[:, 4]  # Assume same freq

    # Smoothing
    win = daniell_points(freq_hz)   # 0.5 THz / 0.08 THz -> 7 points (0.56 THz)
    smoothed_pssum_fluid = daniell_smooth(pssum_fluid, win)
    smoothed_pssum_solid = daniell_smooth(pssum_solid, win)

    # Compute integrals
    int_product = integrate.trapezoid(smoothed_pssum_solid * smoothed_pssum_fluid, freq_hz)
    int_solid = integrate.trapezoid(smoothed_pssum_solid, freq_hz)
    int_fluid = integrate.trapezoid(smoothed_pssum_fluid, freq_hz)
    overlap = int_product / (int_solid * int_fluid) if int_solid * int_fluid != 0 else 0
    overlaps[i] = overlap

    # Plot overlap for this pressure
    fig, ax = plt.subplots(figsize=(10,6))
    ax.plot(fTHz, smoothed_pssum_solid, 'b-', label='Solid')
    ax.plot(fTHz, smoothed_pssum_fluid, 'r-', label='Liquid')
    min_ps = np.minimum(smoothed_pssum_solid, smoothed_pssum_fluid)
    ax.fill_between(fTHz, 0, min_ps, color='gray', alpha=0.3, label='Overlap Area')
    ax.set_xlabel(r'$\omega$ [THz]')
    ax.set_ylabel(r'$\mathrm{q}(\omega)\;\;\;\;[\mathrm{eV}\ \mathrm{ps}^{-1}\ \mathrm{Hz}^{-1}]$')
    ax.set_title(f"SDHF Overlap at {p}")
    ax.axhline(y=0, linestyle='--', color='k', alpha=0.5)
    ax.set_xlim([0., 15])
    ax.legend()
    plt.savefig(f"plots/sdhf/overlap_{p}.pdf")
    plt.close()

# Plot overlap vs pressure
press_num = np.array([int(p.replace('bar', '')) for p in press])
fig, ax = plt.subplots(figsize=(10,6))
ax.plot(press_num, overlaps, 'o-')
ax.set_xlabel('Pressure (bar)')
ax.set_ylabel('SDHF Overlap')
ax.set_title("SDHF Overlap vs Pressure")
plt.savefig("plots/sdhf/overlap_vs_pressure.pdf")
plt.close()

# Save overlaps
np.savetxt('sdhf_overlaps.dat', np.column_stack((press_num, overlaps)), header='pressure overlap', fmt='%d %.4e')

print("All plots saved in plots/sdhf/")
print("Overlaps saved in sdhf_overlaps.dat")
