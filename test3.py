#!/usr/bin/env python3
# -*- coding: iso-8859-15 -*-

import numpy as np
from md import *
import sys
from sys import stderr


Lx, Ly, Lz = 1000, 1000, 1000
grid = Grid(n=10, L=[Lx,Ly,Lz], neighbors_dist=1)

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
Ekin = 1.7E5
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

# Set time step
dt = 0.005

# NDarray that stores all pos and vel data
num_steps = 350
data = np.zeros(shape=(num_steps, num_atoms, 6))

# Loop
step = 0
while step < num_steps:
    # Put atom in grid
    grid.clear()
    for atom in atoms:
        grid.insert(atom)

    # Calculate neighbors for each atom
    for atom in atoms:
        atom.calc_neighbors()

    # Mechanism
    for atom1 in atoms:
        atom1.step1(dt)
        for atom2 in atom1.neighbors:
            atom1.SoftS(atom2, e=1E3)
    for atom in atoms:
        atom.step2(dt)

    # Measurements
    for i, atom in enumerate(atoms):
        data[step, i, :3] = atom.pos
        data[step, i, 3:] = atom.vel

    # Update step count
    stderr.write('\rStep: {:05d}/{:05d}'.format(
        step, num_steps
    ))
    step += 1

# Output
np.save('test3.data', data)
stderr.write('\n')
