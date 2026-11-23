#!/opt/cray/pe/python/3.9.13.1/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import sys

press=sys.argv[1]
start = 1
end = int(sys.argv[2]) # Note: In practice, set 'end' to a larger number divisible by 5 for meaningful block averaging
nsets = end - start + 1

#-------------------------------------------#
# Graphene walls temperature
#-------------------------------------------#
nSheets=9
src=f"../step2_sheet_pos_and_temp/outs/{press}bar"
cool_lwall_pos_temp = np.zeros(( nSheets, 3, nsets ))
hot_lwall_pos_temp = np.zeros(( nSheets, 3, nsets ))
cool_rwall_pos_temp = np.zeros(( nSheets, 3, nsets ))
hot_rwall_pos_temp = np.zeros(( nSheets, 3, nsets ))
for i in range(start-1, end):

    dat1 = np.loadtxt(f"{src}/out{i+1}.C.lwall")
    dat2 = np.loadtxt(f"{src}/out{i+1}.H.lwall")
    dat3 = np.loadtxt(f"{src}/out{i+1}.C.rwall")
    dat4 = np.loadtxt(f"{src}/out{i+1}.H.rwall")

    cool_lwall_pos_temp[:,:,i] = dat1
    hot_lwall_pos_temp[:,:,i] = dat2
    cool_rwall_pos_temp[:,:,i] = dat3
    hot_rwall_pos_temp[:,:,i] = dat4

# Define block averaging functions
def block_average_pos(data, num_blocks):
    pos = data[:,1,:] # Positions in column 1, shape (9, nsets)
    nsets = pos.shape[1]
    block_size = nsets // num_blocks
    mean_pos_blocks = []
    for k in range(num_blocks):
        start = k * block_size
        end = (k + 1) * block_size if k < num_blocks - 1 else nsets
        block_pos = pos[:,start:end]
        mean_pos_k = np.mean(block_pos, axis=1)
        mean_pos_blocks.append(mean_pos_k)
    mean_pos_blocks = np.array(mean_pos_blocks) # (num_blocks, 9)
    overall_mean_pos = np.mean(mean_pos_blocks, axis=0) # (9,)
    error_pos = np.std(mean_pos_blocks, axis=0) # (9,)
    return overall_mean_pos, error_pos

def block_average_temp(data, num_blocks):
    temp = data[:,2,:] # Temperatures in column 2, shape (9, nsets)
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

# Set number of blocks
num_blocks = 5
# Compute block averages
#  right wall
mean_pos_cool_rwall, error_pos_cool_rwall = block_average_pos(cool_rwall_pos_temp, num_blocks)
temp_blocks_cool_rwall = block_average_temp(cool_rwall_pos_temp, num_blocks)
# Hot right wall
mean_pos_hot_rwall, error_pos_hot_rwall = block_average_pos(hot_rwall_pos_temp, num_blocks)
temp_blocks_hot_rwall = block_average_temp(hot_rwall_pos_temp, num_blocks)
#  left wall
mean_pos_cool_lwall, error_pos_cool_lwall = block_average_pos(cool_lwall_pos_temp, num_blocks)
temp_blocks_cool_lwall = block_average_temp(cool_lwall_pos_temp, num_blocks)
# Hot left wall
mean_pos_hot_lwall, error_pos_hot_lwall = block_average_pos(hot_lwall_pos_temp, num_blocks)
temp_blocks_hot_lwall = block_average_temp(hot_lwall_pos_temp, num_blocks)

# Write data to DAT files
with open(f'outs/{press}bar_cool_sheets.dat', 'w') as f:
    f.write('# Right half  sheets\n')
    f.write('# Sheet_Index Position Position_Error Temp_Block1 Temp_Block2 Temp_Block3 Temp_Block4 Temp_Block5\n')
    for i in range(nSheets):
        temps = " ".join(f"{temp_blocks_cool_rwall[k, i]:.4f}" for k in range(num_blocks))
        f.write(f"{i+1} {mean_pos_cool_rwall[i]:.4f} {error_pos_cool_rwall[i]:.4f} {temps}\n")
    f.write('# Left half  sheets\n')
    f.write('# Sheet_Index Position Position_Error Temp_Block1 Temp_Block2 Temp_Block3 Temp_Block4 Temp_Block5\n')
    for i in range(nSheets):
        temps = " ".join(f"{temp_blocks_cool_lwall[k, i]:.4f}" for k in range(num_blocks))
        f.write(f"{i+1} {mean_pos_cool_lwall[i]:.4f} {error_pos_cool_lwall[i]:.4f} {temps}\n")

with open(f'outs/{press}bar_hot_sheets.dat', 'w') as f:
    f.write('# Right half hot sheets\n')
    f.write('# Sheet_Index Position Position_Error Temp_Block1 Temp_Block2 Temp_Block3 Temp_Block4 Temp_Block5\n')
    for i in range(nSheets):
        temps = " ".join(f"{temp_blocks_hot_rwall[k, i]:.4f}" for k in range(num_blocks))
        f.write(f"{i+1} {mean_pos_hot_rwall[i]:.4f} {error_pos_hot_rwall[i]:.4f} {temps}\n")
    f.write('# Left half hot sheets\n')
    f.write('# Sheet_Index Position Position_Error Temp_Block1 Temp_Block2 Temp_Block3 Temp_Block4 Temp_Block5\n')
    for i in range(nSheets):
        temps = " ".join(f"{temp_blocks_hot_lwall[k, i]:.4f}" for k in range(num_blocks))
        f.write(f"{i+1} {mean_pos_hot_lwall[i]:.4f} {error_pos_hot_lwall[i]:.4f} {temps}\n")

print("Data written to cool_sheets.dat and hot_sheets.dat")
