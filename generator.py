# Program Name: generator.py
# Author: Aravinthen Rajkumar
# Description: Generates a polymer configuration

import numpy as np
import math
import random
import matplotlib.pyplot as plt

class SimulationBox:
    def __init__(self, x, y, z):

        # dimensions of box volume
        self.x = x
        self.y = y
        self.z = z

        # number of walks
        self.num_walks = 0
        
        # Contains details on the random walks generated.
        # This is necessary to convert between function scope and global scope.
        # When the beads are inputted into a lammps structure file, they need a
        # global number.
        self.walk_details = {}
        
        self.beads = [] # Contains a list of all required information for beads.
                        # 0. the chain number
                        # 1. the bead number
                        # 2, 3, 4. the positions

    def generateChain(self, index, bondlength, dihedral):
        # this function generates a random walk within a volume
        # periodic boundary conditions are employed

        def bond_vec(bondlength, dihedral, bv):
            # calculates a random bond vector

            # use the dihedral angle to obtain the position angle 
            gamma = math.pi - dihedral

            # randomly generate theta
            theta = random.uniform(0,2*math.pi)

            # calculate bond vector components
            x = bondlength*math.sin(gamma)*math.cos(theta)
            y = bondlength*math.sin(gamma)*math.sin(theta)
            z = bondlength*math.cos(gamma)

            vec = np.array([x,y,z])
            
            # calculate rotation based on whether a previous bond vector
            # was generated. you can't just leave the bond vector as it is,
            # it must be generated with respect to the previous position to
            # avoid an upward bias
            bv = bv*np.linalg.norm(bv)
            # taking the dot product with the z-axis, so just take zth
            # component
            # using v_ to distinguish from bond vector
            v_gamma = math.acos(bv[2]) # rotation over y axis
            v_theta = math.acos(bv[1]) # rotation over z axis
            
            rotmat_y = np.array([[math.cos(v_gamma), 0, math.sin(v_gamma)],
                                 [0, 1, 0],
                                 [-math.sin(v_gamma), 0, math.cos(v_gamma)]])
            
            rotmat_z = np.array([[math.cos(v_theta), -math.sin(v_theta), 0],
                                 [math.sin(v_theta), math.cos(v_theta), 0],
                                 [0, 0, 1]])
                       
            new_vec = np.dot(np.dot(rotmat_z, rotmat_y), vec)
            
            return new_vec

        walk = [] # all beads within this walk will be stored here.
                          # they will later be checked and then read into the
                          # global bead list. Keep function scope separate to
                          # class scope please!
        
        # generate initial position
        posn = [random.uniform(0, self.x),
                random.uniform(0, self.y),
                random.uniform(0, self.z)]

        # generate the bond vector outside to use for the rotation matrix
        seed_bv = np.random.rand(3)
        seed_bv = seed_bv/np.linalg.norm(seed_bv)
        
        bv = bond_vec(bondlength, dihedral, seed_bv)

        # Indices:
        # 0. the chain number
        # 1. the bead number
        # 2, 3, 4. the positions
        walk.append([self.num_walks,
                     0,
                     posn[0],
                     posn[1],
                     posn[2]])
        
        for bead in range(1, index):
            posn = posn + bv

            # account for periodicity
            posn[0] = posn[0] % self.x
            posn[1] = posn[1] % self.y
            posn[2] = posn[2] % self.z

            walk.append([self.num_walks,
                         bead,
                         posn[0],
                         posn[1],
                         posn[2]])

            prev_bv = bv            
            bv = bond_vec(bondlength, dihedral, bv)

        for i in walk:
            self.beads.append(i)

        # add global walk parameters to walk_details
        self.walk_details[self.num_walks] = index
        
        # increment number of walks.
        self.num_walks += 1            


box = SimulationBox(30, 30, 30)
box.generateChain(5, 0.1, math.pi/2)

x = []
y = []
z = []

print("testing angles")
for j in range(1, len(box.beads)):
    b1 = np.array([box.beads[j][2] - box.beads[j-1][2],
                   box.beads[j][3] - box.beads[j-1][3],
                   box.beads[j][4] - box.beads[j-1][4]])

    b2 = np.array([box.beads[j-1][2] - box.beads[j-2][2],
                   box.beads[j-1][3] - box.beads[j-2][3],
                   box.beads[j-1][4] - box.beads[j-2][4]])
    
    print(np.linalg.norm(b1),
          np.linalg.norm(b2))
    
    print(math.acos(np.dot(np.array(b1),np.array(b2))/(np.linalg.norm(b1)*np.linalg.norm(b2))))
    print(" ")

for i in box.beads:
    x.append(i[2])
    y.append(i[3])
    z.append(i[4])

fig = plt.figure()
ax = fig.add_subplot(projection='3d')

plt.plot(x, y, z, "-")
plt.show()
