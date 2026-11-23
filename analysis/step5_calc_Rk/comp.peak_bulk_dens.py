#!/opt/cray/pe/python/3.9.13.1/bin/python3
import numpy as np
import scipy
import scipy.signal

pressures = [50, 100, 200, 250, 300, 350, 400]
bulk_bins = 50 # Number of central bins to average over for bulk
search_range = 4     # Angstroms to search for the first peak from the interface
offset = 0.0  # Offset from interface to start searching
nbins_per_half = 5000

# Collect data
data = np.zeros((len(pressures), 5))
for idx, p in enumerate(pressures):
    # Load density data
    dens_file = f"outs/{p}bar_pfh_dens.dat"
    dens_data = np.loadtxt(dens_file)
    dens = 0.5 * ( dens_data[0:nbins_per_half, :] + dens_data[nbins_per_half:, :] )
    pos = dens_data[0:nbins_per_half, 0]

    # Average density profiles
    avg_dens = np.mean( dens[:, 1:], axis=1 )

    # Bulk density (average of left and right central)
    central_start = (nbins_per_half - bulk_bins) // 2
    central_end = central_start + bulk_bins
    bulk_dens = np.mean( avg_dens[central_start:central_end] )

    # Load hot sheets to get hot interface positions
    hot_sheets_file = f"outs/{p}bar_hot_sheets.dat"
    hot_sheets = np.loadtxt(hot_sheets_file)
    right_hot_interface = hot_sheets[8, 1]  # Right inner hot sheet

    # Right peak (positions increase from hot)
    start_idx   = np.searchsorted( pos, right_hot_interface + offset )
    end_idx     = np.searchsorted( pos, right_hot_interface + search_range)
    print( "pressure = %d"%(p) )
    print( right_hot_interface )
    print( pos[start_idx], pos[end_idx] )
    sub_dens    = avg_dens[start_idx:end_idx]

    print( start_idx, end_idx ) #, sub_dens )
    if len(sub_dens) > 0:

        peaks = scipy.signal.find_peaks(sub_dens)[0]
        if len(peaks) > 0:

            sub_sub_dens = sub_dens[peaks]
            arg = np.argmax( sub_sub_dens )     
            # peak_pos = pos[peak_idx]
            # peak_dens_arg = sub_sub_dens[arg]
            # print( "Index corresponding to max density =%d"%(arg) )
            # print( "Position corresponding to max density =%14.5f"%(peak_pos) )
            # print( "Max density in the numerous peaks is =%14.5f"%(peak_dens) )
            # print( "Max density in the numerous peaks using arg =%14.5f"%(peak_dens_arg) )
            # print( "--------------------------------------------------------" )
            # print( "                                                        " )

            rel_idx = peaks[arg]  # Select the first detected peak (leftmost in the search range)
            peak_idx = start_idx + rel_idx
            peak_dens = np.max( sub_sub_dens )
            peak_dist = pos[peak_idx] - right_hot_interface
        else:
            peak_dens = np.nan
            peak_dist = np.nan
    else:
        peak_dens = np.nan
        peak_dist = np.nan

    data[idx] = [p, bulk_dens, peak_dens, peak_dist, peak_dens/bulk_dens]

# Save to file
header = '# pressure bulk_dens peak_dens peak_dist reduced_dens'
np.savetxt('density_stats.dat', data, header=header, fmt=['%d', '%.2f', '%.2f', '%.2f', '%.3f'])
print("Data written to density_stats.dat")
