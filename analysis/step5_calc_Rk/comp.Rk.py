#!/opt/cray/pe/python/3.9.13.1/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import sys
import scipy
from scipy.stats import linregress

# Declare variables
j = float(sys.argv[1]) # in units of W/m^2
press = sys.argv[2] # in bar

#------------------
# wall temperatures
#------------------
nSheets = 9
fitRegion = 100
skipPlot = 20 # skip points for plotting purposes
skipSheets = 3 # skip the sheets subject to thermalization

# load all relevant files
cool_sheets = np.loadtxt( f"outs/{press}bar_cool_sheets.dat" )
hot_sheets = np.loadtxt( f"outs/{press}bar_hot_sheets.dat" )
pfh_temp = np.loadtxt( f"outs/{press}bar_pfh_temp.dat" )
pfh_dens = np.loadtxt( f"outs/{press}bar_pfh_dens.dat" )
nBins = int(pfh_temp.shape[0] / 2)
fitStart = int((nBins - fitRegion)/2)
fitEnd = int((nBins + fitRegion)/2)
nBins2 = int(pfh_dens.shape[0] / 2)

# Compute averages over blocks
hot_sheets_pos = hot_sheets[0:nSheets,1] # right-side domain positions (mean)
hot_sheets_R_temp_av = np.mean(hot_sheets[0:nSheets,3:8], axis=1)
hot_sheets_L_temp_av = np.mean(hot_sheets[nSheets:,3:8], axis=1)
hot_sheets_temp_av = 0.5 * (hot_sheets_R_temp_av + hot_sheets_L_temp_av)
cool_sheets_pos = cool_sheets[0:nSheets,1] # right-side domain positions (mean)
cool_sheets_R_temp_av = np.mean(cool_sheets[0:nSheets,3:8], axis=1)
cool_sheets_L_temp_av = np.mean(cool_sheets[nSheets:,3:8], axis=1)
cool_sheets_temp_av = 0.5 * (cool_sheets_R_temp_av + cool_sheets_L_temp_av)
pfh_pos = pfh_temp[0:nBins,0]
pfh_temp_av = 0.5 * (np.mean(pfh_temp[0:nBins,1:6], axis=1) + np.mean(pfh_temp[nBins:,1:6], axis=1))
pfh_bins = pfh_dens[0:nBins2,0]
pfh_dens_av = 0.5 * (np.mean(pfh_dens[0:nBins2,1:6], axis=1) + np.mean(pfh_dens[nBins2:,1:6], axis=1))

#------ Linear Fitting ------#
# hot sheets
slope, intercept, r_value, p_value, std_err = linregress( hot_sheets_pos[skipSheets:], hot_sheets_temp_av[skipSheets:] )
hot_fit = slope * hot_sheets_pos[skipSheets:] + intercept
r2_fit_hot = r_value * r_value
# cool sheets
slope, intercept, r_value, p_value, std_err = linregress( cool_sheets_pos[skipSheets:], cool_sheets_temp_av[skipSheets:] )
cool_fit = slope * cool_sheets_pos[skipSheets:] + intercept
r2_fit_cool = r_value * r_value
# pfh temp profile
slope, intercept, r_value, p_value, std_err = linregress( pfh_pos[fitStart:fitEnd], pfh_temp_av[fitStart:fitEnd] )
pfh_temp_fit = slope * pfh_pos + intercept
r2_fit_pfh = r_value * r_value
dTdx = 1e10 * abs((pfh_temp_fit[-1]-pfh_temp_fit[0]) / (pfh_pos[-1]-pfh_pos[0])) # in units of K/m (1e10 is for conversion from Ang to m)
pfh_therm_cond = j / dTdx # W/m^2 / (K/m) = W/mK
# Computing Rk - calculate at both hot and cool interface
# first step determine dT_hot and dT_cool
pfh_temp_hot = slope * hot_sheets_pos[-1] + intercept # extrapolate the linear temperature profile to the hot graphene sheet
pfh_temp_cool = slope * cool_sheets_pos[-1] + intercept # extrapolate the linear temperature profile to the cool graphene sheet
dT_hot = hot_sheets_temp_av[-1] - pfh_temp_hot # in units of K
dT_cool = pfh_temp_cool - cool_sheets_temp_av[-1] # in units of K
Rk_hot = dT_hot / j
Rk_cool = dT_cool / j
Rk_av = 0.5 * (Rk_hot + Rk_cool)

#------ Plotting ------#
fig, ax = plt.subplots(figsize=(12,6))
ax2 = ax.twinx()
# Plotting graphene sheets temperature
ax.plot( hot_sheets_pos, hot_sheets_temp_av, 'D--', label='Hot Graphene', c='k', markerfacecolor='none' )
ax.plot( cool_sheets_pos, cool_sheets_temp_av, 's--', label='Cool Graphene', c='k', markerfacecolor='none' )
ax.plot( hot_sheets_pos[skipSheets:], hot_fit, c='k', alpha=0.5, linestyle='--' )
ax.plot( cool_sheets_pos[skipSheets:], cool_fit, c='k', alpha=0.5, linestyle='--' )
ax.plot( pfh_pos[skipPlot:-skipPlot], pfh_temp_av[skipPlot:-skipPlot], 'r.', label='PFH' )
ax.plot( pfh_pos, pfh_temp_fit, c='coral', linestyle='--', linewidth=2.5 )
# plot pfh density
ax2.plot( pfh_bins, pfh_dens_av, 'b-', label='Density' )
ax2.tick_params(axis='y', color='blue', labelcolor='blue')
ax2.set_ylabel('Density [kg/m$^3$]', color='b')
# Highlight the background from x=2 to x=5
ax.axvspan( hot_sheets_pos[0]-1.5, hot_sheets_pos[2]+1.5, color='red', alpha=0.1 )
ax.axvspan( cool_sheets_pos[0]-1.5, cool_sheets_pos[2]+1.5, color='blue', alpha=0.1 )
ax.set_xlabel('Position (Å)')
ax.set_ylabel('Temperature (K)')
ax.set_title(f"Graphene PFH temperature profile at {press} bar (Averaged over blocks)")
ax.legend()
plt.savefig(f"plots/{press}bar_average.pdf")

# plt.show()
# Write average stats into a file
with open(f'outs/{press}bar_kapitza_average.dat', 'w') as f:
    f.write('# r2_hot r2_cool r2_pfh therm_cond dT_hot dT_cool Rk_hot Rk_cool Rk_av\n')
    f.write(f"{r2_fit_hot:.4f} {r2_fit_cool:.4f} {r2_fit_pfh:.4f} {pfh_therm_cond:.4e} {dT_hot:.4f} {dT_cool:.4f} {Rk_hot:.4e} {Rk_cool:.4e} {Rk_av:.4e}\n")
