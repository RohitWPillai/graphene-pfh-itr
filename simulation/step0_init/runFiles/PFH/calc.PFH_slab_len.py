#!/usr/bin/python3
import sys
import math

# Check for command-line arguments; if not provided, use default values.
if len(sys.argv) < 4:
    lx = 2 * 24.93      # Default lateral dimension in x (Å)
    ly = 2 * 21.56      # Default lateral dimension in y (Å)
    slab_width = 80.0   # Default pfh slab thickness in Å
else:
    lx = float(sys.argv[1])
    ly = float(sys.argv[2])
    slab_width = float(sys.argv[3])

# Constants
mPFH = (12.011*6 + 18.998*14)*1e-3 # kg/mol (mass of one mole of pfh)
nAv = 6.022e23      # 1/mol (Avogadro's number)
rho = 1.76e3        # kg/m^3 (density; pfh ~1e3 kg/m^3)

# Compute the volume of the pfh slab in Å^3
V_box = lx * ly * slab_width

# The mass of pfh required (in kg) equals density * volume in m^3.
# Note: 1 m^3 = 1e30 Å^3.
mass_required = rho * (V_box / 1e30)

# Number of moles required = mass_required / mPFH
moles_required = mass_required / mPFH

# Number of pfh molecules required = moles_required * nAv
N = moles_required * nAv

# Round N to the nearest 100.
N_rounded = int(round(N / 100.0) * 100)

# Output only the rounded number.
print(N_rounded)
