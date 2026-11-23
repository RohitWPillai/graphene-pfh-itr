#!/opt/cray/pe/python/3.9.13.1/bin/python3
import numpy as np
import scipy
from scipy.stats import linregress
import matplotlib.pyplot as plt
import os
import sys
import subprocess

pressures = ['50', '100', '200', '250', '300', '350', '400']
fitFrom = 500
eV_to_J = 1.60218e-19  # eV to Joule
ps_to_s = 1e-12  # ps to s
ang2_to_m2 = 1e-20  # angstrom^2 to m^2

os.makedirs('plots', exist_ok=True)

with open('output.dat', 'w') as f:
    f.write('# pressure\tR2\theat_flux_evps\tJ_Wm2\n')
    for press in pressures:

        subCmd   = f"grep 'BOUNDS' /scratch-nvme/<project>/<project>/<user>/PFH/{press}bar/rep_1/dump1.lpfh -m 1 -A 3 |"
        command1 =  subCmd + "head -2 | tail -1 | awk '{print $2}'"
        command2 =  subCmd + "head -3 | tail -1 | awk '{print $2}'"
        result1 = subprocess.run( command1, shell=True, capture_output=True, text=True )
        result2 = subprocess.run( command2, shell=True, capture_output=True, text=True )

        # print( result1 )
        # print( result2 )

        # print( result1.stdout.strip )
        # print( result2.stdout.strip )

        xhi = float( result1.stdout.strip() )
        yhi = float( result2.stdout.strip() )

        Area_ang2 = (2*xhi)*(2*yhi)     # in Ang^2
        convert = eV_to_J / (ps_to_s * Area_ang2 * ang2_to_m2)  # conversion from eV/ps to J/m^2 s

        dat1 = np.loadtxt(f"qInOut/out.{press}_bar")
        t1 = np.arange(0, dat1.shape[0]) * 1e4 * 1e-3  # in ps
        hO1 = dat1[:,1]  # in ev
        hI1 = dat1[:,2]  # in ev
        avg = 0.5 * (-hO1 + hI1)
        slope, intercept, r_value, p_value, std_err = linregress(t1[fitFrom:], avg[fitFrom:])
        power_eVps = slope * 0.5  # in ev/ps
        power = power_eVps * convert
        r_squared = r_value * r_value
        f.write(f"{press}\t{r_squared:.5f}\t{power_eVps:.5e}\t{power:.5e}\n")

        fig, ax = plt.subplots(figsize=(10,8))
        ax.plot(t1*1e-3, -hO1, 'b-', label=r'$E_{In}$', linewidth=2)
        ax.plot(t1*1e-3, hI1, 'c-', label=r'-$E_{Out}$', linewidth=2)
        ax.plot(t1*1e-3, avg, 'g-', label=r'Average', linewidth=2)
        linhIO = slope * t1 + intercept
        ax.plot(t1*1e-3, linhIO, 'k-', label=r'Linear Fit', linewidth=2)
        ax.axvline(x=t1[fitFrom]*1e-3, color='gray', linestyle='--', label='Transient upto here')
        ax.legend(fontsize=16)
        ax.set_title(f"Press={press} bar, Slope=Q={power_eVps:.3f} eV/ps", fontsize=16)
        ax.set_xlabel(r"Time [ns]", fontsize=14)
        ax.set_ylabel(r"Energy In/Out [eV]", fontsize=14)
        plt.savefig(f"plots/qInOut_{press}bar.pdf")
        plt.close()
