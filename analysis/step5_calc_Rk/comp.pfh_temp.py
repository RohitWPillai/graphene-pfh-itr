#!/opt/cray/pe/python/3.9.13.1/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import sys

press=sys.argv[1]
start = 1
end = int(sys.argv[2]) # Note: In practice, set 'end' to a larger number divisible by 5 for meaningful block averaging
nsets = end - start + 1
nbins=250
num_blocks=5

src=f"../step3_pfh_temp/outs/{press}bar"
#-----------------------------#
# PFH temperature
#-----------------------------#
lpfh_temp = np.zeros((nbins, nsets))
rpfh_temp = np.zeros((nbins, nsets))
for i in range( nsets ):
    left_pfh = np.loadtxt( f"{src}/out{i+1}.lpfh.temp" )
    right_pfh = np.loadtxt( f"{src}/out{i+1}.rpfh.temp" )
    lpfh_temp[:,i] = left_pfh[:,1]
    rpfh_temp[:,i] = right_pfh[:,1]
bin_pos = right_pfh[:,0]
def block_average_temp(data, num_blocks):
    temp = data # Temperatures in column 2, shape (9, nsets)
    nsets = temp.shape[1]
    block_size = nsets // num_blocks
    mean_temp_blocks = []
    for k in range(num_blocks):
        start = k * block_size
        end = (k + 1) * block_size if k < num_blocks - 1 else nsets
        block_temp = temp[:,start:end]
        mean_temp_k = np.mean(block_temp, axis=1)
        mean_temp_blocks.append(mean_temp_k)
    mean_temp_blocks = np.array(mean_temp_blocks) # (num_blocks, 9)
    return mean_temp_blocks
# left pfh (temperatures only)
temp_blocks_lpfh = block_average_temp( lpfh_temp, num_blocks )
temp_blocks_rpfh = block_average_temp( rpfh_temp, num_blocks )
# Write data to DAT files
with open(f'outs/{press}bar_pfh_temp.dat', 'w') as f:
    f.write('# Right half PFH\n')
    f.write('# Bin_Position Temp_Block1 Temp_Block2 Temp_Block3 Temp_Block4 Temp_Block5\n')
    for i in range(nbins):
        temps = " ".join(f"{temp_blocks_rpfh[k, i]:.4f}" for k in range(num_blocks))
        f.write(f"{bin_pos[i]:.4f} {temps}\n")
    f.write('# Left half PFH\n')
    f.write('# Bin_Position Temp_Block1 Temp_Block2 Temp_Block3 Temp_Block4 Temp_Block5\n')
    for i in range(nbins):
        temps = " ".join(f"{temp_blocks_lpfh[k, i]:.4f}" for k in range(num_blocks))
        f.write(f"{bin_pos[i]:.4f} {temps}\n")
print("Data written to pfh_temp.dat")
