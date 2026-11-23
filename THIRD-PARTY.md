# Third-party files

## Tersoff potential file (`SiC.tersoff`, `in.SiC.tersoff`)

72 identical copies, one in each simulation stage folder (list them with `git ls-files | grep SiC.tersoff`).

- **Source:** the `potentials/SiC.tersoff` file distributed with LAMMPS (https://www.lammps.org), as its header
  records: contributor Aidan Thompson; values from J. Tersoff, Phys. Rev. B 39, 5566 (1989) and its erratum,
  Phys. Rev. B 41, 3248 (1990); file from Rutuparna Narulkar.
- **Modified:** the C-C-C entry carries the optimised Tersoff parameters of L. Lindsay and D. A. Broido,
  Phys. Rev. B 81, 205441 (2010) (B = 430, h = -0.930), in place of the original carbon values. The Si entries
  are as distributed and are not used by these simulations.
- **Licence:** LAMMPS is distributed under the GNU General Public License, version 2 (GPL-2.0). These files keep
  that licence; its text is in `LICENSES/GPL-2.0.txt`. The repository's MIT and CC BY 4.0 licences do not apply
  to them. Each copy carries a notice of the change on its second line (the date of the change, in 2025, was not
  recorded).

## moltemplate input text (`pfh.lt`, `system.lt`, `README_setup.sh`)

The moltemplate input files and setup notes under `simulation/*/step0_init/**/PFH/` follow the examples
distributed with moltemplate (https://moltemplate.org), and some of their comment text may be adapted from those
examples. moltemplate is under the MIT licence, copyright (c) 2013 Andrew Jewett, University of California
Santa Barbara; its licence text is in `LICENSES/moltemplate-MIT.txt`.

## Not included

moltemplate's OPLS-AA library `oplsaa2024.lt`, which `pfh.lt` imports to rebuild the starting configuration, is
distributed with moltemplate (https://moltemplate.org) and is not included here. The force-field coefficients
used in the simulations are written out in each stage's `in.forcefield`.
