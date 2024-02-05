# Program name: example.py

import sys
import math
sys.path.insert(0, '../main/')
from generator import SimulationBox

box = SimulationBox(50, 50, 50)

for i in range(4):
    box.generateChain(50, 1.53, math.pi*(109.1/180))

box.structure("test_structure")
box.settings("test_settings")

pairwise_coeffs = [0.1984, 3.62390, 10.5]
bond_coeffs = [350, 1.53]
angle_coeffs = [50, 109.1]
torsion_coeffs = [1, -3.0, 0, 4, 0]

box.DreidingUA(pairwise_coeffs, bond_coeffs, angle_coeffs, torsion_coeffs)

box.equilibrate(10000,
                0.5,
                500)

box.quench(50000,
           0.5,
           500,
           500,)

box.quench(50000,
           0.5,
           500,
           100)

box.quench(50000,
           0.5,
           100,
           100)

box.run(lammps_path="~/Research/research-code/test_software/lammps/src/lmp_serial")
