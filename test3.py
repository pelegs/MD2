#!/usr/bin/env python3
# -*- coding: iso-8859-15 -*-

import numpy as np
from md import *
import sys
from sys import stderr


Lx, Ly, Lz = 1000, 1000, 1000
grid = Grid(n=5, L=[Lx,Ly,Lz], neighbors_dist=1)

# Particles
num_atoms = len(grid.cells)
positions = [cell.center for cell in grid.cells]

atoms = [Atom(
              grid=grid,
              pos=positions[id],
              vel=np.zeros(3),
              rad=10,
              mass=1,
              color = [175,175,175],
              id=id,
         )
         for id in range(num_atoms)
]

# Set kinetic energy to specific val
Ekin = 1.7E4
Ekin_perAtom = Ekin/num_atoms
for atom in atoms:
    atom.vel = set_norm(
        np.random.uniform(-1, 1, 3),
        np.sqrt(2*Ekin_perAtom*atom.m_)
    )

# Set CM momentum to zero
CMP = np.sum(np.array([atom.calc_P() for atom in atoms]), axis=0)
CMP_perAtom = CMP/num_atoms
for atom in atoms:
    atom.vel = atom.vel - CMP_perAtom

# Loop
frame = 0
num_frames = 5000
while frame < num_frames:
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
        atom1.step1(dt=0.005)
        for atom2 in atom1.neighbors:
            atom1.SoftS(atom2, e=1E3)
    for atom in atoms:
        atom.step2(dt=0.005)

#    # Measurements
#    total_KE = sum([atom.calc_KE() for atom in atoms])
#    CMP = np.sum(np.array([atom.calc_P() for atom in atoms]), axis=0)
#
    # Output
#    print(frame, total_KE, CMP)
    stderr.write('\rFrame: {:04d}'.format(frame))
    print(num_atoms)
    print('Frame: {}'.format(frame))
    for atom in atoms:
        print('He {}'.format(' '.join(map(str, atom.pos))))
stderr.write('\n')
