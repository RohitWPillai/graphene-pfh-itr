# Data files: columns and units

The files in this folder are kept as produced. Two of them carry wrong units in their header lines; the
correct units are in the table below.

| File | Columns and units | Header correction |
|---|---|---|
| `Rk_pfh.dat` | pressure (bar), R_K mean (m² K/W), R_K standard deviation over 5 blocks (m² K/W) | the header says `[K/m^2W]`; the unit is m² K/W |
| `heat_flux_J.dat` | pressure (bar), R² of the linear fit, heat flow per path (eV/ps), J (W/m²) | |
| `density_stats.dat` | pressure (bar), bulk density (kg/m³), first-peak density at the hot interface (kg/m³), peak distance from the graphene sheet (Å), ratio of the two densities | |
| `heat_flux_qInOut/out.<P>_bar` | timestep (1 fs), cumulative energy exchanged by the hot thermostat (eV; negative = added), by the cold thermostat (eV); every 10 ps | |
| `per_pressure/<P>bar_hot_sheets.dat`, `_cool_sheets.dat` | sheet index, mean position (Å), position error (Å), temperature in blocks 1 to 5 (K) | |
| `per_pressure/<P>bar_pfh_temp.dat` | bin position (Å), temperature in blocks 1 to 5 (K) | |
| `per_pressure/<P>bar_pfh_dens.dat` | bin position (Å), density in blocks 1 to 5 (kg/m³) | |
| `per_pressure/<P>bar_kapitza_stats.dat` | block, R² of the hot, cold and PFH fits, PFH thermal conductivity (W/m K), temperature jump hot and cold (K), R_K hot, cold and mean (m² K/W) | |
| `per_pressure/<P>bar_kapitza_average.dat` | the same quantities from the block-averaged profiles | |
| `spectral_heat_flux/<P>bar_fluid.dat`, `_solid.dat` | frequency (Hz), spectral heat flux in x, y, z and total (eV ps⁻¹ Hz⁻¹); integrated over frequency it gives the heat flow per path in eV/ps; the solid side has the opposite sign | the header says `[eV/s]`; the unit is eV ps⁻¹ Hz⁻¹ |
| `spectral_heat_flux/sdhf_overlaps.dat` | pressure (bar), overlap ∫SF / (∫S ∫F) (1/Hz) | |
| `sheet_positions/positions_means.dat`, `positions_uncertainties.dat` | pressure (bar), mean z of the left-cold, left-hot, right-hot and right-cold contact sheets and their uncertainties (Å) | |
| `summary.json` | the main results per pressure, with units; written by `tools/verify_results.py --write` | |

Row layouts (which rows belong to which slab or domain) are described in `AGENTS.md`.
