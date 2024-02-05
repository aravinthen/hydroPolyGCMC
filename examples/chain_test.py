# Program name: example.py

import sys
import math
sys.path.insert(0, '../main/')
from generator import SimulationBox

box = SimulationBox(5, 5, 5)
 
for i in range(4):
    box.generateChain(50, 1.53, math.pi*(109.1/180), 2)

box.structure("test_structure")

