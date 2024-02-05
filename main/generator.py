# Program Name: generator.py
# Author: Aravinthen Rajkumar
# Description: Generates a polymer configuration

import numpy as np
import math
import random
import matplotlib.pyplot as plt
from datetime import date
import os
today = date.today()

class SimulationBox:
    def __init__(self, x, y, z, cellnums):

        # dimensions of box volume
        self.x = x
        self.y = y
        self.z = z
        self.z = cellnums
        self.box_size = [x,y,z]

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

        # checks and figures
        self.structure_ready = False
        self.structure_name = None
        self.settings_file = None

        # potential info
        self.pairwise = None
        self.bond = None
        self.angle = None
        self.torsion = None

        self.Cells = []
        index_list = ((i,j,k)
                      for i in range(self.numsx)
                      for j in range(self.numsy)
                      for k in range(self.numsz))
        
        for index in index_list:
            cell_positions = np.array([round(self.x*index[0], 16),
                                       round(self.y*index[1], 16),
                                       round(self.z*index[2], 16)])
            self.Cells.append(self.Cell(index[0], index[1], index[2], cell_positions, self.cellside))

        

    def generateChain(self, index, bondlength, dihedral, cutoff):
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
        self.beads.append([self.num_walks,
                           0,
                           posn[0],
                           posn[1],
                           posn[2]])
        
        for bead in range(1, index):
            trial_bv = bond_vec(bondlength, dihedral, bv)
            trial_posn = posn + trial_bv

            # account for periodicity
            trial_posn[0] = trial_posn[0] % self.x
            trial_posn[1] = trial_posn[1] % self.y
            trial_posn[2] = trial_posn[2] % self.z

            attempts = 0
            issue = 0
            valid = False
            while valid==False:
                for ex_beads in self.beads:
                    delta = trial_posn - [ex_beads[-3], ex_beads[-2], ex_beads[-1]]
                    for i in range(3):
                        if delta[i] > 0.5*self.box_size[i]:
                            delta[i] -= self.box_size[i]
                        elif (delta[i] <= -0.5*self.box_size[i]):
                            delta[i] += self.box_size[i]
                            
                    distance = np.linalg.norm(delta)

                    if distance < cutoff:
                        issue=+1
                    
                if issue == 0:
                    valid = True
                else:
                    trial_bv = bond_vec(bondlength, dihedral, bv)
                    trial_posn = posn + trial_bv
                    # account for periodicity
                    trial_posn[0] = trial_posn[0] % self.x
                    trial_posn[1] = trial_posn[1] % self.y
                    trial_posn[2] = trial_posn[2] % self.z
                    issue = 0
                    attempts+=1

            bv = trial_bv
            posn = trial_posn
            self.beads.append([self.num_walks,
                               bead,
                               posn[0],
                               posn[1],
                               posn[2]])

        # add global walk parameters to walk_details
        self.walk_details[self.num_walks] = index
        
        # increment number of walks.
        self.num_walks += 1            

    def structure(self,
                  filename):
        
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
{len(self.beads)-2*self.num_walks} angles                                           \n\
{len(self.beads)-3*self.num_walks} dihedrals                                        \n\
                                                                                    \n\
1 atom types                                                                        \n\
1 bond types                                                                        \n\
1 angle types                                                                       \n\
1 dihedral types                                                                    \n\
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
        
        # Bond writing
        global_num = 0
        bond_num = 0
        for walk in self.walk_details:
            for bond in range(self.walk_details[walk]-1):
                bond_num+=1
                global_num+=1
                f.write(f"\t{global_num}\t 1 \t{bond_num}\t{bond_num+1}\n")
            bond_num+=1

        f.write("\n\
Angles   \n\
   \n")
        a_global_num = 0
        angle_num = 0
        for walk in self.walk_details:
            for angle in range(self.walk_details[walk]-2):
                angle_num+=1
                a_global_num+=1
                f.write(f"\t{a_global_num}\t 1 \t{angle_num}\t{angle_num+1}\t{angle_num+2}\n")
            angle_num+=2

        f.write("\n\
Dihedrals   \n\
   \n")

        d_global_num = 0
        dihedral_num = 0
        for walk in self.walk_details:
            for angle in range(self.walk_details[walk]-3):
                dihedral_num+=1
                d_global_num+=1
                f.write(f"\t{d_global_num}\t 1 \t{dihedral_num}\t{dihedral_num+1}\t{dihedral_num+2}\t{dihedral_num+3}\n")
            dihedral_num+=3
        
        f.close()
        
        self.structure_ready = True
        self.structure_name = filename

    def settings(self,
                 filename,
                 nlist=[10,10000],
                 nskin=1.0,):
        """
        Initializes the settings file.
        """
        if self.structure_ready == False:
            raise EnvironmentError("You must create a structure file before settings can be defined.")

        self.settings_file = filename
        f = open(self.settings_file, "w")
        f.write(f"\
#-----------------------------------------------------------------------------------\n\
# HYDROPOLY - SIMULATION FILE                                                       \n\
#-----------------------------------------------------------------------------------\n\
# FILENAME:  {filename}                                                             \n\
# DATE: {today}                                                                     \n\
#-----------------------------------------------------------------------------------\n\
units           real                                                                \n\
boundary	p p p                                                               \n\
atom_style	molecular                                                           \n\
log 		log.{self.settings_file}.txt                                        \n\
read_data	{self.structure_name}                                               \n\
                                                                                    \n\
neighbor	0.4 bin                                                             \n\
neigh_modify	every {nlist[0]} one {nlist[1]}                                     \n\
")
        f.close()


        
    def DreidingUA(self,
                   pairwise,
                   bond,
                   angle,
                   torsion):
        """
        Initializes the settings file.
        """
        
        self.pairwise = pairwise
        self.bond = bond
        self.angle = angle
        self.torsion = torsion
        
        eps, sig, cutoff = pairwise[0], pairwise[1], pairwise[2]
        Kb, r0 = bond[0], bond[1]
        Ka, theta0 = angle[0], angle[1]
        A1, A2, A3, A4, A5 = torsion[0], torsion[1], torsion[2], torsion[3], torsion[4]
        
        f = open(self.settings_file, "a")
        f.write(f"\n\
# Potential - DREIDING_UA                                                           \n\
pair_style	lj/cut {cutoff}                                                     \n\
pair_coeff	1 1 {eps} {sig} {cutoff}                                            \n\
bond_style      harmonic                                                            \n\
bond_coeff	1 {Kb} {r0}                                                         \n\
angle_style     harmonic                                                            \n\
angle_coeff	1 {Kb} {theta0}                                                     \n\
dihedral_style	multi/harmonic                                                      \n\
dihedral_coeff	1 {A1} {A2} {A3} {A4} {A5}                                          \n\
")        
        f.close()

    def equilibrate(self,
                    runtime,
                    timestep,
                    temp,
                    temp_f=None,
                    limit=0.05,
                    thermo= 1000,
                    damp= 10.0,
                    dump= 1000,
                    seed=random.randrange(0,99999)):

        if temp_f == None:
            temp_f = temp
        
        f = open(self.settings_file, "a")
        f.write(f"\n\
# Equilibration Stage                                                                         \n\
velocity 	all create {temp} {seed}                                                      \n\
fix		1 all nve/limit {limit}                                                       \n\
fix		2 all langevin {temp} {temp_f} {damp} {seed}                               \n\
thermo_style	custom step temp                                                              \n\
thermo          {thermo}                                                                      \n\
timestep	1                                                                             \n\
run		{runtime}                                                                     \n\
unfix 1                                                                                       \n\
unfix 2                                                                                       \n\
")        
        f.close()

    def quench(self,
               runtime,
               timestep,
               qtemp,
               qtemp_f,
               limit=0.5,
               thermo=1000,
               tdamp=50.0,
               pdamp=1000.0,
               dump=1000,
               drag=2.0,
               seed=random.randrange(0,99999)):
        
        f = open(self.settings_file, "a")
        f.write(f"\n\
dump            1 all cfg {dump} dump.{self.settings_file}_*.cfg mass type xs ys zs fx fy fz  \n\
velocity 	all create {qtemp} 1231                                                       \n\
fix		3 all npt temp {qtemp} {qtemp_f} {tdamp} iso 0.0 0.0 {pdamp} drag {drag}      \n\
fix		4 all momentum 1 linear 1 1 1                                                 \n\
thermo_style	custom step temp  lx ly lz press pxx pyy pzz                                  \n\
thermo          {thermo}                                                                      \n\
timestep	1                                                                             \n\
run		{runtime}                                                                     \n\
unfix 3                                                                                       \n\
unfix 4                                                                                       \n\
undump          1                                                                             \n\
")        
        f.close()
        
    def run(self, lammps_path, mpi=0):
        os.system(f"{lammps_path} < {self.settings_file}")

# Program name: example.py
box = SimulationBox(100, 100, 100)
 
for i in range(5):
    print(f"WALK {i}")
    box.generateChain(2000, 1.53, math.pi*(109.1/180), 0.5)

box.structure("test_structure")
