# Program Name: generator.py
# Author: Aravinthen Rajkumar
# Description: Generates a polymer configuration

import numpy as np
import math
import random
import matplotlib.pyplot as plt
from datetime import date
today = date.today()

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
            # avoid an upward bias.

            # essentially, we are transforming our coordinate system such
            # that this is the new z-axis
            n = bv/np.linalg.norm(bv)
            k = np.array([0,0,1]) # old x-axis
            
            # we do so using quaternion transformations
            theta = math.acos(np.dot(n, k))
            b = np.cross(k, n)
            b = b/np.linalg.norm(b)

            # quaternion components
            q0 = math.cos(theta/2)
            q1 = math.sin(theta/2)*b[0]
            q2 = math.sin(theta/2)*b[1]
            q3 = math.sin(theta/2)*b[2]

            Q = np.array([[ q0*q0 + q1*q1 - q2*q2 - q3*q3 , 2*(q1*q2 - q0*q3), 2*(q1*q3 + q0*q2)],
                          [ 2*(q1*q2 + q0*q3), q0*q0 - q1*q1 + q2*q2 - q3*q3, 2*(q2*q3 - q0*q1)],
                          [ 2*(q1*q3 - q0*q2),  2*(q3*q2 + q0*q1), q0*q0 - q1*q1 - q2*q2 + q3*q3]])

            # this should work fine!
            new_vec = np.dot(Q, vec)
            
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

    def structure(self, filename):
        f = open(filename, "w")
        
        f.write(f"\
#-----------------------------------------------------------------------------------\n\
# HYDROPOLY - STRUCTURAL DATA FILE                                                  \n\
#-----------------------------------------------------------------------------------\n\
# FILENAME:  {filename}                                                        \n\
# DATE: {today}                                                                     \n\
#-----------------------------------------------------------------------------------\n\
                                                                                    \n\
{len(self.beads)} atoms                                                             \n\
{len(self.beads)-self.num_walks} bonds                                              \n\
                                                                                    \n\
1 atom types                                                                        \n\
1 bond types                                                                        \n\
                                                                                    \n\
0.0000 {self.x} xlo xhi                                                             \n\
0.0000 {self.y} ylo yhi                                                             \n\
0.0000 {self.z} zlo zhi                                                             \n\
                                                                                    \n\
                                                                                    \n\
Masses                                                                              \n\
                                                                                    \n\
1 1.0                                                                               \n\
                                                                                    \n\
Atoms                                                                               \n\
                                                                                    \n\
")

        for bead in range(len(self.beads)):
            struc = self.beads[bead][0]+1
            x, y, z = self.beads[bead][-3], self.beads[bead][-2], self.beads[bead][-1]
            f.write(f"\t{bead+1}\t{struc}\t 1 \t{x:.5f}\t\t{y:.5f}\t\t{z:.5f}\t\n")

        f.write("\n\
Bonds   \n\
   \n")
        
        global_num = 0
        bond_num = 0
        for walk in self.walk_details:
            for bond in range(self.walk_details[walk]-1):
                bond_num+=1
                global_num+=1
                f.write(f"\t{global_num}\t 1 \t{bond_num}\t{bond_num+1}\n")
            bond_num+=1
                    
            
        f.close()






