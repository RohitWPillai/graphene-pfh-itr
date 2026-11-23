#!/opt/cray/pe/python/3.9.13.1/bin/python3

import numpy as np
import matplotlib.pyplot as plt

pressures_pfh = [50, 100, 200, 250, 300, 350, 400]
wall_pos = np.loadtxt("../step2_sheet_pos_and_temp/positions_means.dat")

nBins=10000
plotUpto=int(nBins/2)

i=0
fig, ax = plt.subplots(figsize=(8, 6))
for p in pressures_pfh:

    dat = np.loadtxt(f"outs/{p}bar_pfh_dens.dat")
    refPos = wall_pos[i,3]

    ax.plot( dat[0:plotUpto:10,0]-refPos, np.mean( dat[0:plotUpto:10,1:], axis=1), label=f'{p} bar', linestyle='-', marker='none' )

ax.set_xlim([0.,35])
ax.legend()
ax.set_xlabel(r"Distance from graphene wall [Ang]")
ax.set_ylabel(r"$\rho$ [kg/m^3]")

fig.savefig("plots/pfh_dens_vs_press.pdf")
