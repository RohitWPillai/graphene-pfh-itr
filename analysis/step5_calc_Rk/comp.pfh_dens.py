#!/opt/cray/pe/python/3.9.13.1/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import sys

press=sys.argv[1]
start = 1
end = int(sys.argv[2]) # Note: In practice, set 'end' to a larger number divisible by 5 for meaningful block averaging
nsets = end - start + 1
nbins=5000
num_blocks=5
src=f"../step4_pfh_dens/outs/{press}bar"
#-----------------------------#
# PFH Densities
#-----------------------------#
lpfh_dens = np.zeros((nbins, nsets))
rpfh_dens = np.zeros((nbins, nsets))
for i in range( nsets ):
    left_pfh = np.loadtxt( f"{src}/out{i+1}.lpfh.dens" )
    right_pfh = np.loadtxt( f"{src}/out{i+1}.rpfh.dens" )
    lpfh_dens[:,i] = left_pfh[:,1]
    rpfh_dens[:,i] = right_pfh[:,1]
bin_pos = right_pfh[:,0]
def block_average_dens(data, num_blocks):
    dens = data # Densities in column 2, shape (9, nsets)
    nsets = dens.shape[1]
    block_size = nsets // num_blocks
    mean_dens_blocks = []
    for k in range(num_blocks):
        start = k * block_size
        end = (k + 1) * block_size if k < num_blocks - 1 else nsets
        block_dens = dens[:,start:end]
        mean_dens_k = np.mean(block_dens, axis=1)
        mean_dens_blocks.append(mean_dens_k)
    mean_dens_blocks = np.array(mean_dens_blocks) # (num_blocks, 9)
    return mean_dens_blocks
# left pfh (densities only)
dens_blocks_lpfh = block_average_dens( lpfh_dens, num_blocks )
dens_blocks_rpfh = block_average_dens( rpfh_dens, num_blocks )
# Write data to DAT files
with open(f'outs/{press}bar_pfh_dens.dat', 'w') as f:
    f.write('# Right half PFH\n')
    f.write('# Bin_Position Dens_Block1 Dens_Block2 Dens_Block3 Dens_Block4 Dens_Block5\n')
    for i in range(nbins):
        denss = " ".join(f"{dens_blocks_rpfh[k, i]:.4f}" for k in range(num_blocks))
        f.write(f"{bin_pos[i]:.4f} {denss}\n")
    f.write('# Left half PFH\n')
    f.write('# Bin_Position Dens_Block1 Dens_Block2 Dens_Block3 Dens_Block4 Dens_Block5\n')
    for i in range(nbins):
        denss = " ".join(f"{dens_blocks_lpfh[k, i]:.4f}" for k in range(num_blocks))
        f.write(f"{bin_pos[i]:.4f} {denss}\n")
print("Data written to pfh_dens.dat")
