#!/opt/cray/pe/python/3.9.13.1/bin/python3
import numpy as np

pressures = [50, 100, 200, 250, 300, 350, 400]
num_blocks = 5
bulk_bins = 50  # Number of central bins to average over for bulk

# Function to compute bulk mean and std for a quantity (temp or dens)
def compute_bulk_mean_std(data_file, nbins_per_half):
    data = np.loadtxt(data_file)
    # Assuming data shape (2*nbins_per_half, 1 + num_blocks): col0 pos, col1-5 blocks
    right_data = data[0:nbins_per_half, :]
    left_data = data[nbins_per_half:, :]
    right_pos = right_data[:, 0]
    left_pos = left_data[:, 0]
    right_blocks = right_data[:, 1:].T  # shape (num_blocks, nbins_per_half)
    left_blocks = left_data[:, 1:].T  # shape (num_blocks, nbins_per_half)

    # Central indices
    central_start = (nbins_per_half - bulk_bins) // 2
    central_end = central_start + bulk_bins

    bulk_per_block = np.zeros(num_blocks)
    for b in range(num_blocks):
        right_bulk = np.mean(right_blocks[b, central_start:central_end])
        left_bulk = np.mean(left_blocks[b, central_start:central_end])
        bulk_per_block[b] = 0.5 * (right_bulk + left_bulk)

    mean_bulk = np.mean(bulk_per_block)
    std_bulk = np.std(bulk_per_block)
    return mean_bulk, std_bulk

# Collect data
data = np.zeros((len(pressures), 7))
for idx, p in enumerate(pressures):
    # T_bulk from pfh_temp.dat, nbins_per_half=250
    temp_file = f"outs/{p}bar_pfh_temp.dat"
    t_bulk_mean, t_bulk_std = compute_bulk_mean_std(temp_file, 250)

    # density_bulk from pfh_dens.dat, nbins_per_half=5000
    dens_file = f"outs/{p}bar_pfh_dens.dat"
    dens_bulk_mean, dens_bulk_std = compute_bulk_mean_std(dens_file, 5000)

    # thermal_conductivity from kapitza_stats.dat
    stats_file = f"outs/{p}bar_kapitza_stats.dat"
    stats_data = np.loadtxt(stats_file, skiprows=1)
    therm_cond_blocks = stats_data[:, 4]  # column 4: therm_cond
    therm_cond_mean = np.mean(therm_cond_blocks)
    therm_cond_std = np.std(therm_cond_blocks)

    data[idx] = [p, t_bulk_mean, t_bulk_std, dens_bulk_mean, dens_bulk_std, therm_cond_mean, therm_cond_std]

# Save to file
header = '# pressure T_bulk_mean T_bulk_std density_bulk_mean density_bulk_std therm_cond_mean therm_cond_std'
np.savetxt('bulk_stats.dat', data, header=header, fmt=['%d', '%.2f', '%.2f', '%.2f', '%.2f', '%.4e', '%.4e'])

print("Data written to bulk_stats.dat")
