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
flux_data = np.loadtxt('../step1_heat_flux/output.dat', skiprows=1)
lammps_flux = flux_data[:, 3]  # J_Wm2
lammps_flux_evps = flux_data[:, 2]  # ev/ps

os.makedirs('plots/sdhf', exist_ok=True)

integrals = np.zeros(len(press))
inplane_pcts = np.zeros(len(press))
outplane_pcts = np.zeros(len(press))

fig_comp, ax_comp = plt.subplots(figsize=(10,6))

for i, p in enumerate(press):
    # Assume data files are named like "outs/sdhf_{p}.txt" - adjust if needed
    dat = np.loadtxt(f"sdhf_consolidated/{p}_fluid.dat", skiprows=1)
    freq_hz = dat[:, 0]
    ps_x = dat[:, 1]
    ps_y = dat[:, 2]
    ps_z = dat[:, 3]
    pssum = dat[:, 4]

    fTHz = freq_hz / 1e12

    # Smoothing
    win = daniell_points(freq_hz)   # 0.5 THz / 0.08 THz -> 7 points (0.56 THz)
    smoothed_ps_x = daniell_smooth(ps_x, win)
    smoothed_ps_y = daniell_smooth(ps_y, win)
    smoothed_ps_z = daniell_smooth(ps_z, win)
    smoothed_pssum = daniell_smooth(pssum, win)

    # Plot 1: q(omega) vs omega for this pressure
    fig, ax = plt.subplots(figsize=(10,6))
    ax.plot(fTHz, pssum, 'k-', alpha=0.2, label='Raw')
    ax.plot(fTHz, smoothed_pssum, 'k-', label='Smoothed')
    ax.set_xlabel(r'$\omega$ [THz]')
    ax.set_ylabel(r'$\mathrm{q}(\omega)\;\;\;\;[\mathrm{eV}\ \mathrm{ps}^{-1}\ \mathrm{Hz}^{-1}]$')
    ax.set_title(f"SDHF for PFH at {p}")
    ax.axhline(y=0, linestyle='--', color='k', alpha=0.5)
    ax.set_xlim([0., 15])
    ax.set_ylim([-0.5e-14, 5e-14])
    ax.legend()
    plt.savefig(f"plots/sdhf/sdhf_{p}.pdf")
    plt.close()

    # Compute integral of q(omega) domega (in ev/ps)
    integral = integrate.trapezoid(smoothed_pssum, freq_hz)
    integrals[i] = integral

    # Cumulative integrals
    cum_integral_x = integrate.cumulative_trapezoid(smoothed_ps_x, freq_hz, initial=0)
    cum_integral_y = integrate.cumulative_trapezoid(smoothed_ps_y, freq_hz, initial=0)
    cum_integral_z = integrate.cumulative_trapezoid(smoothed_ps_z, freq_hz, initial=0)
    cum_integral_inplane = cum_integral_x + cum_integral_y
    cum_integral_total = cum_integral_inplane + cum_integral_z

    # Plot 2: Cumulative for this pressure
    fig, ax = plt.subplots(figsize=(10,6))
    ax.plot(fTHz, cum_integral_inplane, color='r', linestyle='-', label='in-plane')
    ax.plot(fTHz, cum_integral_z, color='b', linestyle='-', label='out-of-plane z')
    ax.plot(fTHz, cum_integral_total, color='k', linestyle='-', label='Total')
    ax.axhline(y=lammps_flux_evps[i], linestyle='--', color='k', alpha=0.5, label="LAMMPS Heat flux")
    ax.legend()
    ax.set_title(f"In-plane and out-of-plane cumulative heat flux for PFH at {p}")
    ax.set_xlabel(r'$\omega$ [THz]')
    ax.set_ylabel(r'$\int d\omega q(\omega)$ / [eV/ps]')
    ax.set_xlim([0., 40])
    plt.savefig(f"plots/sdhf/cumulative_{p}.pdf")
    plt.close()

    # Percentages at >30 THz (find index where fTHz > 30)
    idx_30 = np.searchsorted(fTHz, 30)
    total_at_30 = cum_integral_total[idx_30 - 1] if idx_30 > 0 else cum_integral_total[-1]
    inplane_at_30 = cum_integral_inplane[idx_30 - 1] if idx_30 > 0 else cum_integral_inplane[-1]
    outplane_at_30 = cum_integral_z[idx_30 - 1] if idx_30 > 0 else cum_integral_z[-1]
    inplane_pcts[i] = (inplane_at_30 / total_at_30) * 100 if total_at_30 != 0 else 0
    outplane_pcts[i] = (outplane_at_30 / total_at_30) * 100 if total_at_30 != 0 else 0

    # For comparative plot
    ax_comp.plot(fTHz, smoothed_pssum, label=p, alpha=0.7)

# Comparative plot of all SDHFs
ax_comp.set_xlabel(r'$\omega$ [THz]')
ax_comp.set_ylabel(r'$\mathrm{q}(\omega)\;\;\;\;[\mathrm{eV}\ \mathrm{ps}^{-1}\ \mathrm{Hz}^{-1}]$')
ax_comp.set_title("Comparative SDHF for all pressures")
ax_comp.axhline(y=0, linestyle='--', color='k', alpha=0.5)
ax_comp.set_xlim([0., 15])
ax_comp.set_ylim([-0.5e-14, 5e-14])
ax_comp.legend()
plt.savefig("plots/sdhf/comparative_sdhf.pdf")
plt.close()

# Plot integrals vs LAMMPS flux
fig, ax = plt.subplots(figsize=(10,6))
ax.plot(press, integrals, 'o-', label='Integral from SDHF')
ax.plot(press, lammps_flux_evps, 's-', label='LAMMPS flux')
ax.set_xlabel('Pressure')
ax.set_ylabel('Heat flux [eV/ps]')
ax.set_title("Integral of SDHF vs LAMMPS flux")
ax.legend()
plt.savefig("plots/sdhf/integral_vs_lammps.pdf")
plt.close()

# Plot in-plane and out-of-plane percentages
press_num = np.array([int(p.replace('bar', '')) for p in press])
fig, ax = plt.subplots(figsize=(10,6))
width = 20  # bar width
ax.bar(press_num - width/2, inplane_pcts, width=width, label='In-plane', color='r')
ax.bar(press_num + width/2, outplane_pcts, width=width, label='Out-of-plane', color='b')
ax.set_xlabel('Pressure (bar)')
ax.set_ylabel('Percentage (%)')
ax.set_title("In-plane vs Out-of-plane Percentages by Pressure")
ax.legend()
plt.savefig("plots/sdhf/percentages_comparison.pdf")
plt.close()

# Save percentages
#pct_data = np.column_stack((press, inplane_pcts, outplane_pcts))
#np.savetxt('percentages.dat', pct_data, header='pressure inplane_pct outplane_pct', fmt='%s %.2f %.2f')

print("All plots saved in plots/sdhf/")
print("Percentages saved in percentages.dat")
