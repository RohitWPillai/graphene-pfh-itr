#!/opt/cray/pe/python/3.9.13.1/bin/python3

import numpy as np

pressures = [50, 100, 200, 250, 300, 350, 400]
nsets = 80  # Assuming from bash script, adjust if needed
n_blocks = 5
block_size = nsets // n_blocks

assert nsets % n_blocks == 0, "nsets must be divisible by n_blocks"

# Prepare output arrays
means = np.zeros((len(pressures), 5))  # pressure, lc, lh, rc, rh
uncs = np.zeros((len(pressures), 4))  # lc, lh, rc, rh

for idx, p in enumerate(pressures):
    cold_lwall_pos = np.zeros((9, 3, nsets))
    hot_lwall_pos = np.zeros((9, 3, nsets))
    cold_rwall_pos = np.zeros((9, 3, nsets))
    hot_rwall_pos = np.zeros((9, 3, nsets))

    for i in range(nsets):
        dat1 = np.loadtxt(f"outs/{p}bar/out{i+1}.C.lwall")
        dat2 = np.loadtxt(f"outs/{p}bar/out{i+1}.H.lwall")
        dat3 = np.loadtxt(f"outs/{p}bar/out{i+1}.C.rwall")
        dat4 = np.loadtxt(f"outs/{p}bar/out{i+1}.H.rwall")
        cold_lwall_pos[:, :, i] = dat1
        hot_lwall_pos[:, :, i] = dat2
        cold_rwall_pos[:, :, i] = dat3
        hot_rwall_pos[:, :, i] = dat4

    # Extract the relevant positions
    lc_pos = cold_lwall_pos[-1, 1, :]
    lh_pos = hot_lwall_pos[-1, 1, :]
    rc_pos = cold_rwall_pos[-1, 1, :]
    rh_pos = hot_rwall_pos[-1, 1, :]

    # Compute block means
    lc_block_means = np.array([np.mean(lc_pos[j*block_size:(j+1)*block_size]) for j in range(n_blocks)])
    lh_block_means = np.array([np.mean(lh_pos[j*block_size:(j+1)*block_size]) for j in range(n_blocks)])
    rc_block_means = np.array([np.mean(rc_pos[j*block_size:(j+1)*block_size]) for j in range(n_blocks)])
    rh_block_means = np.array([np.mean(rh_pos[j*block_size:(j+1)*block_size]) for j in range(n_blocks)])

    # Overall means
    means[idx, 0] = p
    means[idx, 1] = np.mean(lc_block_means)
    means[idx, 2] = np.mean(lh_block_means)
    means[idx, 3] = np.mean(rh_block_means)
    means[idx, 4] = np.mean(rc_block_means)

    # Uncertainties (std across block means)
    uncs[idx, 0] = np.std(lc_block_means)
    uncs[idx, 1] = np.std(lh_block_means)
    uncs[idx, 2] = np.std(rh_block_means)
    uncs[idx, 3] = np.std(rc_block_means)

# Save means
np.savetxt('positions_means.dat', means, fmt=['%d', '%10.2f', '%10.2f', '%10.2f', '%10.2f'],
           header='pressure\tleft_cold\tleft_hot\tright_hot\tright_cold')

# Save uncertainties
np.savetxt('positions_uncertainties.dat', np.column_stack((pressures, uncs)),
           fmt=['%d', '%10.3f', '%10.3f', '%10.3f', '%10.3f'],
           header='pressure\tunc_left_cold\tunc_left_hot\tunc_right_hot\tunc_right_cold')
