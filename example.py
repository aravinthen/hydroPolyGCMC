# Program name: example.py
# 

import sys
import math
from generator import SimulationBox

box = SimulationBox(30, 30, 30)

for i in range(10):
    box.generateChain(50, 0.5, math.pi*(i/10))

box.structure("test")
