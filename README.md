# Heat flow across graphene / n-perfluorohexane interfaces (NEMD, 50 to 400 bar)

Non-equilibrium molecular dynamics (NEMD, LAMMPS) of steady heat flow across interfaces between multilayer
graphene and liquid n-perfluorohexane (PFH, C6F14), a flexible six-carbon chain molecule, at seven pressures from
50 to 400 bar. The repository contains the simulation inputs kept from the runs (changes are listed at the end),
the analysis code, the reduced data, the figures, and a script that recomputes the numbers in this README from the
data.

The results are the interfacial thermal (Kapitza) resistance R_K, the layering of PFH against graphene, and the
spectral heat flux on both sides of the interface. R_K is the temperature jump at the interface divided by the
heat flux across it; the spectral heat flux is that flux resolved by vibrational frequency. These results are
preliminary and have not been peer reviewed.

Authors: Thejas Hulikal Chakrapani and Rohit Pillai, School of Engineering, University of Edinburgh.
Simulations run on ARCHER2, the UK national supercomputer, in 2025.

## Check the numbers yourself

```
python3 tools/verify_results.py
```

This takes under a minute and needs Python 3.9 or later with numpy, scipy and matplotlib. It re-runs the
shipped R_K, density and spectral-overlap scripts on a copy of `data/`, recomputes R_K and the heat flux with
independent code, and checks the two results tables and the numbers quoted in the results text and the known
issues against the data. It prints PASS or FAIL for each check.

## Simulated system

- **System.** A periodic box contains two PFH slabs, each 400 molecules, separated by stacks of multilayer graphene.
  Each stack is 18 sheets thick, made of two 9-sheet walls back to back. The three outer sheets of each wall
  (six sheets per stack, in the middle of the stack) are thermostatted: one stack at 350 K, the other at 250 K.
  Heat flows from the hot stack, through six free graphene sheets, across a graphene/PFH interface, through
  the PFH, across a second interface and through six more free sheets into the cold stack. Each slab has one
  hot and one cold interface, so there are four interfaces in all. Everything outside the thermostatted sheets
  is integrated without a thermostat (NVE).
- **Size.** The simulated system has 44,800 atoms: 16,000 in PFH and 28,800 carbon atoms in 36 graphene sheets.
  It is made by placing two copies of the 22,400-atom configuration built in step 0 (stages are listed under
  Simulation stages; `simulation/initial_config/combined_system.dat`, box 49.2 × 42.6 × 122.8 Å) end to end. The lateral box
  relaxes during the first NPT stage; the heat-flux cross-section is 2150 Å² (2152 Å² at 100 bar). The gap between
  the graphene sheets that face each other across a slab is 65.5 Å at 50 bar and 61.9 Å at 400 bar.
- **Force fields** (all in `in.forcefield`):
  - PFH: OPLS-AA for perfluoroalkanes (Watkins and Jorgensen, J. Phys. Chem. A 105, 4118, 2001), with partial
    charges and PPPM electrostatics; built with moltemplate and packmol.
  - Carbon-carbon within a graphene sheet: optimised Tersoff (Lindsay and Broido, Phys. Rev. B 81, 205441, 2010).
  - Carbon-carbon between sheets: Lennard-Jones, ε = 2.39 meV, σ = 3.412 Å (Girifalco et al., Phys. Rev. B 62,
    13104, 2000).
  - Graphene-PFH: Lennard-Jones (`lj/smooth/linear`, cut-off 15 Å) with Lorentz-Berthelot mixing.
- **Pressure.** Set in the first NPT stage by a barostat on the box length normal to the walls; the box is then
  fixed at the mean length and the rest of the run is at constant volume.

## Simulation stages

The driver (`simulation/submit_workflow`, ARCHER2 SLURM) runs these stages for each pressure, with a 1 fs timestep
from step 2 on:

| Stage | What it does | Length |
|---|---|---|
| step0_init | build graphene, pack PFH, combine (one 22,400-atom cell) | |
| step1_npt | relax laterally at 0 bar, compress normal to the walls to the target pressure, average the box length, fix it | 5.2 ns (3.2 ns at 100 bar) |
| step2_nvt | equilibrate at 300 K | 2 ns |
| step3_nve | duplicate the cell into the two-slab system, equilibrate | 2 ns |
| step4_therm | 300 K start, then thermostat the hot and cold sheets until the heat flow is steady | 1 ns + 7.5 ns |
| step5_dump | production, positions and velocities every 5 fs | 4 × 20 dump files × 250 ps = 20 ns |

The input files kept from each run are in `simulation/as_run/pressure_<P>bar/`. They are the same at all
pressures except 100 bar, which was run in an earlier batch: its starting configuration was built separately
(in LAMMPS `real` rather than `metal` units at step 0, with no explicit packmol seed), its box-averaging stage
in step 1 is 1 ns rather than 3 ns, a few `recenter` fixes are in different places, and its dumps were written
to a folder without `rep_1/`. The force-field coefficients are the same. Whether it used the same velocity seeds
as the other pressures is not recorded in the files.

## Analysis

- **Heat flux J.** The energy added at the hot stack and the energy removed at the cold stack grow linearly once
  the flow is steady. The heat flow is the slope of their average, fitted after the first 5 ns of thermostatting; half of it crosses each
  interface (the heat leaves the hot stack in both directions). J is that half divided by the cross-section
  (`analysis/step1_heat_flux`).
- **Temperature jump.** The PFH temperature profile, averaged over the two mirror-image slabs, is fitted with a
  straight line over its central 100 of 250 bins and extrapolated to the position of the graphene sheet in
  contact with the PFH. The jump is the difference between that sheet's temperature and the extrapolated PFH
  temperature (`analysis/step5_calc_Rk/comp.Rk.block.py`).
- **R_K** = jump / J, computed separately at the hot and the cold interface and then averaged. The interface
  sheets are at about 341 K (hot) and 259 K (cold), so the averaged R_K mixes two temperatures.
- **Errors.** The 20 ns of production are split into 5 blocks of 4 ns. The error bar is the standard deviation of
  the 5 block values (population, ddof = 0), as in the figures. The standard error of the mean is half of it.
- **Spectral heat flux** (`analysis/step6_sdhf`). The heat flux across the two hot interfaces (the only
  ones whose forces are written to the dumps) is resolved by frequency on the fluid and the solid side, with a
  frequency spacing of 0.08 THz. It is smoothed with a 7-point moving-average (Daniell) window, 0.56 THz wide,
  the nearest to the 0.5 THz chosen. Smoothing changes the integrated flux by less than 0.7 %.

[docs/analysis-choices.md](docs/analysis-choices.md) lists every analysis choice, its value and where it is set.

## Results

R_K, from `data/Rk_pfh.dat` and `data/per_pressure/*_kapitza_stats.dat`, in 10⁻⁸ m² K/W. Each pressure is the
target of the preparation stage (see Limitations).

| Pressure (bar) | R_K (mean ± sd of 5 blocks) | R_K hot interface | R_K cold interface |
|---|---|---|---|
| 50 | 4.29 ± 0.21 | 3.90 | 4.69 |
| 100 | 4.34 ± 0.28 | 4.05 | 4.63 |
| 200 | 4.02 ± 0.21 | 3.54 | 4.50 |
| 250 | 4.14 ± 0.09 | 3.81 | 4.46 |
| 300 | 4.13 ± 0.21 | 3.57 | 4.69 |
| 350 | 3.86 ± 0.26 | 3.33 | 4.39 |
| 400 | 4.01 ± 0.38 | 3.46 | 4.57 |

PFH density, from `data/density_stats.dat`. The bulk density is the mean over the central 50 of 5000 bins, about
0.6 Å of liquid. The first peak is the density maximum of the PFH layer next to the graphene sheet, measured at
the hot interface:

| Pressure (bar) | Bulk density (kg/m³) | First-peak density (kg/m³) | Peak distance from the graphene sheet (Å) |
|---|---|---|---|
| 50 | 1647 | 2584 | 3.18 |
| 100 | 1676 | 2647 | 3.19 |
| 200 | 1711 | 2704 | 3.16 |
| 250 | 1732 | 2768 | 3.15 |
| 300 | 1756 | 2798 | 3.14 |
| 350 | 1776 | 2836 | 3.11 |
| 400 | 1791 | 2861 | 3.12 |

From 50 to 400 bar the bulk density rises by 9 % and the first density peak by 11 %, while R_K falls by 6.5 %,
from 4.29 to 4.01 × 10⁻⁸ m² K/W. The fall is small against the scatter: a weighted straight-line fit through the
seven pressures gives a fall of 0.10 × 10⁻⁸ m² K/W per 100 bar, which is 1.2 times its uncertainty when the block
standard deviations are used as error bars and 2.5 times when the standard errors are used. A pressure dependence is
therefore suggested, not established. R_K at the hot interface is 13 to 24 % lower than at the cold one,
depending on pressure, and lower in 34 of the 35 time blocks. The steady heat flux was 4.2 to 4.7 × 10⁸ W/m²
(`data/heat_flux_J.dat`). Integrated over frequency, the spectral heat flux equals the heat flow measured at
the thermostats to within 1.6 % at every pressure, on both sides of the interface.

Key figures: `figures/rk_vs_pressure.pdf`, `figures/pfh_dens_vs_press.pdf`,
`figures/spectral_heat_flux/overlap_vs_pressure.pdf`, `figures/temperature_profiles/`.

## Limitations

- **One run per pressure.** The error bars come from five consecutive 4 ns blocks of a single run, so they may
  understate the uncertainty. They do not include the uncertainty of J or of the
  fitting choices.
- **Large driving.** The thermostats are 100 K apart and J is about 4.5 × 10⁸ W/m². No smaller-driving run tests
  whether R_K is independent of the driving. The hot and cold interfaces give different values, which could come
  from their different temperatures or from the driving.
- **Choices not varied.** The thermostat damping and thickness, the system size (21.5 nm² cross-section, a
  6.2 to 6.5 nm liquid gap), the fitting window and the plane where the jump is measured were each tried at one
  value only.
- **Pressure.** Each pressure is the target of the preparation stage (step 1). The pressure the liquid exerts on
  the walls was computed during equilibration (steps 2 and 3), but those outputs are not in this repository,
  and it was not computed during the heat-flow run.
- **Temperature.** Local temperatures are computed from the atoms' kinetic energy without subtracting a local
  flow velocity. The PFH has no net flow, so the error should be small, but this is not shown.
- **Force field.** The PFH model is taken from the literature; this repository does not check its density or
  thermal conductivity against experiment.

## Known issues in the shipped code

The code is kept as run, so that it matches the results. Only the first two issues change a number reported
here, R_K by 0.2 % and 0.02 %:

- `analysis/step1_heat_flux/compute_heat_flux.py` builds its time axis from row numbers rather than the recorded
  timesteps; the energy files repeat a timestep at each restart. Using the recorded timesteps raises the fitted
  heat flow by 0.2 %, which would lower R_K by about the same.
- The PFH temperature and density codes (`analysis/step3_pfh_temp`, `analysis/step4_pfh_dens`) use a carbon mass
  of 12.017; the simulations use 12.011 (the graphene sheet code uses 12.011). Correcting it lowers the PFH
  density by about 0.01 % and the PFH temperatures by 0.015 %, which raises R_K at the hot interface by 0.3 %,
  lowers it at the cold interface by 0.2 % and raises the averaged R_K by 0.02 %.
- The header of `data/Rk_pfh.dat` gives the unit as `[K/m^2W]`; it is m² K/W.
- `analysis/step7_vdos/calc_vdos.sh` passes a 1 fs frame interval; production frames are 5 fs apart, so its
  frequency axis would be five times too large. No vibrational density of states (VDOS) results are included.
- In the drivers, each `if [ $? -ne 0 ]` follows an `echo`, so it cannot catch a failed LAMMPS step.
- Rebuilding the starting configuration needs a packmol binary, the PFH molecule file `c6f14.xyz` and
  moltemplate's OPLS-AA library `oplsaa2024.lt` (imported by `pfh.lt`, distributed with moltemplate), none of
  which is included.
- In `comp.peak_bulk_dens.py` the comment says "first detected peak"; the code takes the highest peak within
  4 Å of the sheet.
- The spectral overlap in `data/spectral_heat_flux/sdhf_overlaps.dat` is ∫SF / (∫S ∫F) over frequency, where S
  and F are the solid-side and fluid-side spectral heat fluxes; it has units of 1/Hz, not a dimensionless fraction.
- `step5_dump` iteration 1 reads `in.groups`; iterations 2 to 4 restart with the groups already defined. The
  iterations otherwise differ only in their dump file ranges.

## Reproducing the results

- **From the data here:** R_K, its block statistics, the density statistics and the heat flux, with
  `tools/verify_results.py`.
- **Needs the trajectories (about 20 ns per pressure, dumped every 5 fs; not included because of their size):**
  the per-block temperature and density profiles (analysis steps 2 to 4), the spectral heat flux (step 6) and
  the VDOS (step 7, code only, no results here).
- **Rerun the simulations:** LAMMPS stable 29 Aug 2024 with the MOLECULE, KSPACE and MANYBODY packages. The
  SLURM scripts are for ARCHER2; set the paths and the account in them for your machine (user paths and the
  project code are replaced by `<user>`, `<project>` and `/path/to/`), including the dump destination `dest` in
  `simulation/template/step5_dump/in.varDefs`, which LAMMPS reads. moltemplate and packmol are needed only to
  rebuild the starting configuration, which is included: to skip the rebuild, copy
  `simulation/initial_config/combined_system.dat` to `simulation/step0_init/prep/combine/`, where the driver
  looks for it. The driver's `pressures` list holds 50, 200, 300 and 400 bar; 250 and 350 bar were run with the
  same inputs, and 100 bar with those in `simulation/as_run/pressure_100bar/`.

## Repository map

```
simulation/
  submit_workflow(_sequential)   drivers: run step0 to step5 for each pressure on ARCHER2 (SLURM)
  template/                      clean inputs copied into each pressure folder
  as_run/pressure_<P>bar/        the input files kept from each run
  step0_init/                    configuration build (graphene generation, PFH packing, combination)
  initial_config/combined_system.dat
                                 starting configuration for 50 and 200 to 400 bar (100 bar used its own,
                                 as_run/pressure_100bar/step1_npt/in.config)
analysis/                        code only (C++, Python, SLURM), steps 1 to 7 in order
data/
  README.md                      columns and units of every data file, with two header corrections
  Rk_pfh.dat, density_stats.dat  summary tables
  heat_flux_J.dat                steady heat flux per pressure (W/m²), from step 1
  heat_flux_qInOut/              thermostat energy in and out (eV), every 10 ps
  per_pressure/                  block-averaged sheet temperatures, PFH temperature and density profiles,
                                 Kapitza statistics (the inputs and outputs of step 5)
  spectral_heat_flux/            spectral heat flux per pressure (fluid and solid side) and overlaps
  sheet_positions/               mean graphene sheet positions and uncertainties
  summary.json                   the main results per pressure in one machine-readable file, written by
                                 tools/verify_results.py
figures/                         PDFs from steps 1, 5 and 6
tools/verify_results.py          the check described above
docs/analysis-choices.md         every analysis and simulation choice, with its value and where it is set
AGENTS.md                        orientation for AI coding agents working in this repository
THIRD-PARTY.md, LICENSES/        the files taken from other software and their licences
```

## Changes from the files as run

The files are as run, with these exceptions, none of which changes a simulation setting or a reported number:

1. User names, home paths, the ARCHER2 project code and an example budget code are replaced by `<user>`,
   `<project>` and `/path/to/` in SLURM headers and comments, executable paths and dump destinations (one of
   which, `dest`, is a LAMMPS variable).
2. `simulation/submit_workflow` had been damaged in transfer: comment markers were missing from 19 lines and its
   indentation had become non-breaking spaces, so bash could not parse it. The markers and spaces are restored;
   `simulation/submit_workflow_sequential` shows the same lines as comments.
3. Material not needed for this study is removed: a plotted comparison in
   `analysis/step5_calc_Rk/cmp.Rk_graphene_pfh.py` and its figure, comments that pointed to files outside this
   repository, a commented-out reference density without a source, an unused input file from an earlier
   workflow (`runFiles/npt/in.fix`), and an earlier, unused version of the spectral heat flux code
   (`analysis/sdhf/`).
4. The `pressures` list in `simulation/submit_workflow_sequential` now names this release's pressures.
5. moltemplate's OPLS-AA parameter library (`oplsaa.lt`, five copies) is removed; it is distributed with
   moltemplate and is needed only to rebuild the starting configuration.
6. Each copy of the Tersoff potential file carries a notice of its modification on its second line, as its
   GPL-2.0 licence requires (see `THIRD-PARTY.md`).
7. `.DS_Store` files are removed.

## Citation and licence

Cite the Zenodo record: all versions, https://doi.org/10.5281/zenodo.24458687; this version (0.1.0),
https://doi.org/10.5281/zenodo.24458688. `CITATION.cff` gives the same. Code (`analysis/`, `simulation/`,
`tools/`) is under the MIT licence (`LICENSE`); data, figures and text are under CC BY 4.0 (`LICENSE-data`). The
Tersoff potential files come from LAMMPS and stay under its GPL-2.0 licence (`THIRD-PARTY.md`).
