#!/usr/bin/env python3
# -*- coding: iso-8859-15 -*-

import numpy as np
from md import *
import sys
from os import system
clear = lambda: system('clear')


Lx, Ly, Lz = 700, 700, 700
grid = Grid(n=20, L=[Lx,Ly,Lz], neighbors_dist=1)

# Particles
ncb = 7
num_atoms = ncb**3
positions = [np.array([x, y, z])
             for x in np.linspace(50, Lx-50, ncb)
             for y in np.linspace(50, Ly-50, ncb)
             for z in np.linspace(50, Lz-50, ncb)
             ]

atoms = [Atom(
              grid=grid,
              pos=positions[id],
              vel=np.random.uniform(-1.0E3, 1.0E3, 3),
              rad=5,
              mass=1,
              color = [175,175,175],
              id=id,
         )
         for id in range(num_atoms)
]

# Set CM momentum to zero
CMP = np.sum(np.array([atom.calc_P() for atom in atoms]), axis=0)
CMP_peratom = CMP/num_atoms
for atom in atoms:
    atom.vel = atom.vel - CMP_peratom

# Set kinetic energy to specific val
Ekin = 4.3E3

# Loop
frame = 0
run = True
while run:
    frame += 1

    # Put atom in grid
    grid.clear()
    for atom in atoms:
        grid.insert(atom)

    # Calculate neighbors for each atom
    for atom in atoms:
        atom.calc_neighbors()

    # Mechanism
    for atom1 in atoms:
        atom1.step1(dt=0.001)
        for atom2 in atom1.neighbors:
            atom1.SoftS(atom2, e=1E3)
    for atom in atoms:
        atom.step2(dt=0.0001)

#    # Measurements
#    total_KE = sum([atom.calc_KE() for atom in atoms])
#    CMP = np.sum(np.array([atom.calc_P() for atom in atoms]), axis=0)
#
    # Output
    print('{} {} {}'.format(frame, total_KE, CMP))

