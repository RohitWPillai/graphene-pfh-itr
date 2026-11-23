# Simulation and analysis choices

Every setting below was chosen by a person, and each one can move the result. The table gives the value used,
the file where it is set, and what it controls. Paths are relative to the repository root; `template/` means
`simulation/template/`. The values are the same at every pressure unless a row says otherwise.

## Simulation

| Choice | Value | Set in | Controls |
|---|---|---|---|
| Graphene sheets per wall | 9 (36 in the two-slab system) | `template/step0_init/graf_pfh.prep_workflow.sh` (`SHEETS`) | wall thickness; how far the thermostats sit from the interface |
| Initial sheet spacing | 3.35 Å | same file (`sheet_sep`) | starting graphite stacking |
| PFH slab thickness at build | 60 Å, 400 molecules per slab | same file (`slab_width`); `calc.PFH_slab_len.py` | liquid thickness; length of the linear part of the temperature profile |
| PFH force field | OPLS-AA for perfluoroalkanes (Watkins and Jorgensen 2001) | `template/*/in.forcefield` | liquid structure, density, vibrational spectrum |
| Intra-sheet carbon | optimised Tersoff (Lindsay and Broido 2010) | `in.forcefield`, `in.SiC.tersoff` (C C C entry, B = 430, h = -0.930) | graphene phonons |
| Inter-sheet carbon | LJ, ε = 2.39 meV, σ = 3.412 Å | `in.forcefield` (`eps_CC`, `sig_CC`) | coupling between graphene sheets |
| Graphene-PFH coupling | LJ `lj/smooth/linear`, Lorentz-Berthelot mixing (geometric ε, arithmetic σ) | `in.forcefield` | the interface itself; the largest single lever on R_K |
| Mixing within PFH | geometric ε and σ (OPLS convention) | `in.forcefield` | PFH liquid |
| Cut-offs | LJ 15 Å, real-space Coulomb 13 Å, PPPM accuracy 1e-6 | `in.forcefield` | accuracy and cost |
| Special bonds | LJ and Coulomb 0, 0, 0.5 | `in.forcefield` | OPLS 1-4 scaling |
| Pressure control | NPT in x and y at 0 bar, then in z at the target pressure; box fixed at the mean z length | `template/step1_npt/in.run` | the pressure each run represents |
| Box-length averaging | 3 ns (1 ns at 100 bar) | `step1_npt/in.varDefs` (`Npnpt`) | precision of the fixed box length, hence of the pressure |
| Equilibration | 2 ns NVT at 300 K, then 2 ns NVE | `step2_nvt`, `step3_nve` | starting state for the heat-flux run |
| Two-slab construction | two copies of the equilibrated cell, end to end | `template/step3_nve/in.run` (`read_data ... shift`) | mirror symmetry: two slabs, four interfaces |
| Thermostatted sheets | 3 outer sheets of each wall (types 4 to 6 and 13 to 15) | `template/step4_therm/in.groups` (`therm`) | where heat enters and leaves |
| Thermostat | Nosé-Hoover (`fix nvt`), damping 0.5 ps | `step4_therm/in.varDefs` (`damp`), `in.run` (`qIn`, `qOut`) | how heat is injected |
| Hot and cold temperatures | 350 K and 250 K | `step4_therm/in.varDefs` (`tempHot`, `tempCold`) | heat flux (about 4.5 × 10⁸ W/m²) and the temperatures of the two interfaces |
| Time to steady state | 7.5 ns of thermostatting before production | `step4_therm/in.varDefs` (`Nitherm`) | whether the profiles are steady |
| Timestep | 1 fs (step 1 ramps from 0.01 fs) | `in.run` of each stage | integration accuracy |
| Production | 80 dump files of 250 ps (20 ns), dumps every 5 fs | `template/step5_dump/in.varDefs` (`Nrunblock`, `Ndump`), driver (`maxDump`, 4 iterations) | statistics; spectral resolution |
| Velocity seeds | 9677132, 6773221, 8700004, 1632452, the same for all pressures in the driver (not recorded for 100 bar) | `simulation/submit_workflow` | independent starts (not independent between pressures) |

## Analysis

| Choice | Value | Set in | Controls |
|---|---|---|---|
| Transient skipped in the heat-flux fit | first 500 samples of 10 ps (5 ns) of the 27.5 ns record | `analysis/step1_heat_flux/compute_heat_flux.py` (`fitFrom`) | J |
| Heat per interface | half the thermostat power | same file (`slope * 0.5`) | J (the hot stack feeds two interfaces) |
| Cross-section | box area from the first production dump (2150 Å²; 2152 Å² at 100 bar) | same file (`BOUNDS`) | J |
| PFH temperature bins | 250 per slab | `analysis/step3_pfh_temp/compute.temp.slurm` (`nbins`) | profile resolution |
| PFH density bins | 5000 per slab | `analysis/step4_pfh_dens/compute.dens.slurm` (`nbins`) | layering resolution |
| Linear-fit window for the PFH temperature | central 100 of 250 bins | `analysis/step5_calc_Rk/comp.Rk.block.py` (`fitRegion`) | slope and extrapolation, hence the jump |
| Where the jump is measured | at the position of the graphene sheet in contact with the PFH (not at a dividing surface in the liquid) | same file (`hot_sheets_pos[-1]`, `cool_sheets_pos[-1]`) | the value of R_K: it includes the gap between the sheet and the first PFH layer |
| Hot and cold interfaces | computed separately, then averaged | same file (`Rk_av`) | R_K reported; the hot one is 13 to 24 % lower |
| Mirror averaging | the two slabs are averaged before fitting | same file | noise |
| Blocks for the error | 5 blocks of 16 dump sets (4 ns each) | `run.sh` (`nSets=80`), `comp.Rk.block.py` (`nBlocks`) | error bar |
| Error bar | population standard deviation of the 5 blocks (ddof = 0) | `plot.final.py`, `cmp.Rk_graphene_pfh.py` (`np.std`) | size of the error bar; the standard error is half of it |
| Bulk density | mean of the central 50 of 5000 bins | `comp.peak_bulk_dens.py` (`bulk_bins`) | bulk density |
| First density peak | highest peak within 4 Å of the hot-wall sheet | `comp.peak_bulk_dens.py` (`search_range`) | peak density and distance |
| Interfaces in the spectral heat flux | the two hot-wall interfaces only (forces tallied for `lhotwall`, `rhotwall`) | `template/step5_dump/in.run.iter1` (`force/tally`) | which interfaces the spectra describe |
| Spectral heat flux smoothing | Daniell (moving-average) window of 0.5 THz, 7 points of 0.08 THz | `analysis/step6_sdhf/plot.sdhf_*.py` | the plotted spectra; changes the integral by less than 0.7 % |
