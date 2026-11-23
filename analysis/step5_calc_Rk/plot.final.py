#!/opt/cray/pe/python/3.9.13.1/bin/python3
import numpy as np
import matplotlib.pyplot as plt

pressures = [50, 100, 200, 250, 300, 350, 400]
n_blocks = 5

# Prepare arrays for means and stds
rk_means = np.zeros(len(pressures))
rk_stds = np.zeros(len(pressures))
rk_hot_means = np.zeros(len(pressures))
rk_hot_stds = np.zeros(len(pressures))
rk_cool_means = np.zeros(len(pressures))
rk_cool_stds = np.zeros(len(pressures))
rk_avgs_from_file = np.zeros(len(pressures))  # From average file
rk_hot_from_file = np.zeros(len(pressures))
rk_cool_from_file = np.zeros(len(pressures))

dt_hot_means = np.zeros(len(pressures))
dt_hot_stds = np.zeros(len(pressures))
dt_cool_means = np.zeros(len(pressures))
dt_cool_stds = np.zeros(len(pressures))

therm_cond_means = np.zeros(len(pressures))
therm_cond_stds = np.zeros(len(pressures))

for idx, p in enumerate(pressures):
    # Load stats file
    stats_file = f"outs/{p}bar_kapitza_stats.dat"
    data = np.loadtxt(stats_file, skiprows=1)  # Skip header

    # Columns: Block(0), r2_hot(1), r2_cool(2), r2_pfh(3), therm_cond(4), dT_hot(5), dT_cool(6), Rk_hot(7), Rk_cool(8), Rk_av(9)
    rk_av_block = data[:, 9]
    rk_hot_block = data[:, 7]
    rk_cool_block = data[:, 8]
    dt_hot_block = data[:, 5]
    dt_cool_block = data[:, 6]
    therm_cond_block = data[:, 4]

    rk_means[idx] = np.mean(rk_av_block)
    rk_stds[idx] = np.std(rk_av_block)
    rk_hot_means[idx] = np.mean(rk_hot_block)
    rk_hot_stds[idx] = np.std(rk_hot_block)
    rk_cool_means[idx] = np.mean(rk_cool_block)
    rk_cool_stds[idx] = np.std(rk_cool_block)

    dt_hot_means[idx] = np.mean(dt_hot_block)
    dt_hot_stds[idx] = np.std(dt_hot_block)
    dt_cool_means[idx] = np.mean(dt_cool_block)
    dt_cool_stds[idx] = np.std(dt_cool_block)

    therm_cond_means[idx] = np.mean(therm_cond_block)
    therm_cond_stds[idx] = np.std(therm_cond_block)

    # Load average file
    avg_file = f"outs/{p}bar_kapitza_average.dat"
    avg_data = np.loadtxt(avg_file, skiprows=1)  # Skip header
    # Columns: r2_hot(0), r2_cool(1), r2_pfh(2), therm_cond(3), dT_hot(4), dT_cool(5), Rk_hot(6), Rk_cool(7), Rk_av(8)
    rk_avgs_from_file[idx] = avg_data[8]
    rk_hot_from_file[idx] = avg_data[6]
    rk_cool_from_file[idx] = avg_data[7]

# Plot (a) Rk vs pressure
fig, ax = plt.subplots(figsize=(8, 6))
# ax.errorbar(pressures, rk_hot_means, yerr=rk_hot_stds, fmt='o-', label='Rk_hot Block Mean ± Std')
# ax.errorbar(pressures, rk_cool_means, yerr=rk_cool_stds, fmt='s-', label='Rk_cool Block Mean ± Std')
ax.errorbar(pressures, rk_means, yerr=rk_stds, fmt='^-', label='Rk_av Block Mean ± Std')
# ax.plot(pressures, rk_hot_from_file, 'o--', label='Rk_hot From Average File')
# ax.plot(pressures, rk_cool_from_file, 's--', label='Rk_cool From Average File')
ax.plot(pressures, rk_avgs_from_file, '^--', label='Rk_av From Average File')
ax.set_xlabel('Pressure (bar)')
ax.set_ylabel('Rk (m$^2$K/W)')
ax.set_title('Kapitza Resistance vs Pressure')
ax.legend()
plt.savefig('plots/rk_vs_pressure.pdf')
plt.close()

# Plot (b) dT_hot and dT_cool vs pressure
fig, ax = plt.subplots(figsize=(8, 6))
ax.errorbar(pressures, dt_hot_means, yerr=dt_hot_stds, fmt='o-', label='dT_hot Mean ± Std')
ax.errorbar(pressures, dt_cool_means, yerr=dt_cool_stds, fmt='s-', label='dT_cool Mean ± Std')
ax.set_xlabel('Pressure (bar)')
ax.set_ylabel('ΔT (K)')
ax.set_title('Temperature Jumps vs Pressure')
ax.legend()
plt.savefig('plots/dt_vs_pressure.pdf')
plt.close()

# Plot (c) pfh_thermal_conductivity vs pressure
fig, ax = plt.subplots(figsize=(8, 6))
ax.errorbar(pressures, therm_cond_means, yerr=therm_cond_stds, fmt='o-')
ax.set_xlabel('Pressure (bar)')
ax.set_ylabel('Thermal Conductivity (W/mK)')
ax.set_title('PFH Thermal Conductivity vs Pressure')
plt.savefig('plots/therm_cond_vs_pressure.pdf')
plt.close()

print("Plots saved to plots/ directory.")
