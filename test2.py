#!/usr/bin/env python3
# -*- coding: iso-8859-15 -*-

import numpy as np
from md import *
import sys
import pygame
from pygame.locals import *
from pygame.color import THECOLORS as COLORS
from os import system
clear = lambda: system('clear')


def draw_cell(surface, cell, fill=[0,0,0]):
    rect = (cell.ul_corner[0], cell.ul_corner[1],
            cell.L[0], cell.L[1])
    pygame.draw.rect(surface, fill, rect)
    pygame.draw.rect(surface, [25, 25, 25], rect, width=3)


def draw_atom(surface, atom):
    pygame.draw.circle(
        surface, atom.color, atom.pos[:2].astype(int), atom.rad
    )


def draw_vec(surface, pos, vec, color=[255,0,0]):
    pygame.draw.line(
        surface, color, pos.astype(int), (pos+vec).astype(int), width=1
    )


width, height = 700, 700
grid = Grid(n=20, L=[width,height,100], neighbors_dist=1)

pygame.init()

fps = 60
fpsClock = pygame.time.Clock()

screen = pygame.display.set_mode((width, height))
cells_bg = pygame.Surface((width, height))

# Grid
for index1d, cell in enumerate(grid.cells):
    draw_cell(cells_bg, cell, [0,0,0])

# Particles
n_sqrt = 15
num_atoms = n_sqrt**2
positions = [np.array([x, y, 50])
             for x in np.linspace(50, width-50, n_sqrt)
             for y in np.linspace(50, height-50, n_sqrt)
             ]

atoms = [Atom(
              grid=grid,
#              pos=np.random.uniform(0, width, 3),
              pos=positions[id],
              vel=np.random.uniform(-150, 150, 3),
              rad=5,
              mass=1,
#              color = np.random.randint(100,255,3),
              color = [175,175,175],
              id=id,
         )
         for id in range(num_atoms)
]

atoms[0].color = [255,0,0]

for atom in atoms:
    atom.pos[2] = 50
    atom.vel[2] = 0.0

# Game loop
run = True
while run:
    for event in pygame.event.get():
        if event.type == QUIT:
            run = False
            break

    # Test neighbors
#    x, y = pygame.mouse.get_pos()
#    atoms[0].pos = np.array([x, y, 50])

    # Put atom in grid
    grid.clear()
    for atom in atoms:
        grid.insert(atom)

    # Calculate neighbors for each atom
    for atom in atoms:
        atom.calc_neighbors()

#    for atom in atoms:
#        if atom in atoms[0].neighbors:
#            atom.color = [0,200,0]
#        elif atom is not atoms[0]:
#            atom.color = [175,175,175]

    # Mechanism
    for atom1 in atoms:
        atom1.step1(dt=0.001)
        for atom2 in atom1.neighbors:
            atom1.SoftS(atom2, e=1E3)
    for atom in atoms:
        atom.step2(dt=0.0001)
    total_KE = sum([atom.calc_KE() for atom in atoms])
#    print('\r{:0.3f}'.format(total_KE))

    # Reset screen
    screen.fill([0,0,0])
    screen.blit(cells_bg, [0,0])

    # Draw
    for atom in atoms:
        draw_atom(screen, atom)

    # Update
    pygame.display.flip()
    fpsClock.tick(fps)

pygame.quit()
