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
grid = Grid(n=10, L=[width,height,width], neighbors_dist=1)

pygame.init()

fps = 60
fpsClock = pygame.time.Clock()

screen = pygame.display.set_mode((width, height))
cells_bg = pygame.Surface((width, height))

# Grid
for index1d, cell in enumerate(grid.cells):
    draw_cell(cells_bg, cell, [0,0,0])

# Particles
num_atoms = 75
atoms = [Atom(
              grid=grid,
              pos=np.random.uniform(0, width, 3),
              vel=np.random.uniform(-1, 1, 3),
              rad=5,
              mass=1,
              color = np.random.randint(100,255,3),
              id=id,
         )
         for id in range(num_atoms)
]

for atom in atoms:
    atom.pos[2] = width
    atom.vel[2] = 0.0

# Game loop
run = True
while run:
    for event in pygame.event.get():
        if event.type == QUIT:
            run = False
            break

    # Put atom in grid
    grid.clear()
    for atom in atoms:
        grid.insert(atom)

    # Calculate neighbors for each atom
    for atom in atoms:
        atom.calc_neighbors()

    # Mechanism
    for atom1 in atoms:
        atom1.step1(dt=0.0001)
        for atom2 in atom1.neighbors:
            atom1.SoftS(atom2, e=1E3)
    for atom in atoms:
        atom.step2(dt=0.0001)
    total_KE = sum([atom.calc_KE() for atom in atoms])
    print('\r{:0.3f}'.format(total_KE))

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
