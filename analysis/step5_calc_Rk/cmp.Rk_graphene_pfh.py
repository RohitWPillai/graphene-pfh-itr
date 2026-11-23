#!/opt/cray/pe/python/3.9.13.1/bin/python3
import numpy as np
import matplotlib.pyplot as plt

# pressures_pfh = [50, 200, 300, 350, 400]
pressures_pfh = [50, 100, 200, 250, 300, 350, 400]

# Prepare arrays for PFH
rk_means = np.zeros(len(pressures_pfh))
rk_stds = np.zeros(len(pressures_pfh))
rk_avgs_from_file = np.zeros(len(pressures_pfh))  # From average file

for idx, p in enumerate(pressures_pfh):
    # Load stats file
    stats_file = f"outs/{p}bar_kapitza_stats.dat"
    data = np.loadtxt(stats_file, skiprows=1)  # Skip header

    # Columns: Block(0), r2_hot(1), r2_cool(2), r2_pfh(3), therm_cond(4), dT_hot(5), dT_cool(6), Rk_hot(7), Rk_cool(8), Rk_av(9)
    rk_av_block = data[:, 9]

    rk_means[idx] = np.mean(rk_av_block)
    rk_stds[idx] = np.std(rk_av_block)

    # Load average file
    avg_file = f"outs/{p}bar_kapitza_average.dat"
    avg_data = np.loadtxt(avg_file, skiprows=1)  # Skip header
    # Columns: r2_hot(0), r2_cool(1), r2_pfh(2), therm_cond(3), dT_hot(4), dT_cool(5), Rk_hot(6), Rk_cool(7), Rk_av(8)
    rk_avgs_from_file[idx] = avg_data[8]


wrtArr = np.array(( pressures_pfh, rk_means, rk_stds ))
wrtArr = wrtArr.transpose()
np.savetxt("Rk_pfh.dat", wrtArr, header="Pressures [bar] Rk_mean [K/m^2W] Rk_std [K/m^2W]")

# Plot Rk vs pressure on semilogx
fig, ax = plt.subplots(figsize=(8, 6))
ax.errorbar(pressures_pfh, rk_means, yerr=rk_stds, fmt='o-', label='PFH Block Mean ± Std')
# ax.plot(pressures_pfh, rk_avgs_from_file, 'o--', label='PFH From Average File')
ax.set_xscale('log')
ax.set_xlabel('Pressure (bar)')
ax.set_ylabel('Rk (m$^2$K/W)')
ax.set_title('Kapitza Resistance vs Pressure (PFH)')
ax.legend(loc='upper left')
plt.savefig('plots/rk_vs_pressure_semilog.pdf')
plt.close()

print("Plot saved to plots/rk_vs_pressure_semilog.pdf")
