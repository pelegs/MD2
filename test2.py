#!/usr/bin/env python3
# -*- coding: iso-8859-15 -*-

import numpy as np
from md import *
import sys
import pygame
from pygame.locals import *


def draw_cell(surface, cell, fill=[0,0,0]):
    rect = (cell.ul_corner[0], cell.ul_corner[1],
            cell.L[0], cell.L[1])
    pygame.draw.rect(surface, fill, rect)
    pygame.draw.rect(surface, [25, 25, 25], rect, width=3)


def draw_atom(surface, atom, color):
    pygame.draw.circle(surface, color, atom.pos[:2].astype(int), atom.rad)


def draw_vec(surface, pos, vec, color=[255,0,0]):
    pygame.draw.line(
        surface, color, pos.astype(int), (pos+vec).astype(int), width=1
    )



width, height = 100, 100
grid = Grid(n=2, L=[width,height,200], neighbors_dist=2)

pygame.init()

fps = 60
fpsClock = pygame.time.Clock()

screen = pygame.display.set_mode((width, height))
cells_bg = pygame.Surface((width, height))

# Grid
for index1d, cell in enumerate(grid.cells):
    draw_cell(cells_bg, cell, [0,0,0])

# Particles
num_atoms = 10
atoms = [Atom(
              pos=np.random.uniform(0, width, 3),
              vel=np.random.uniform(-1000, 1000, 3),
              rad=5,
              mass=1,
         )
         for _ in range(num_atoms)
]

forces = []

# Game loop
run = True
while run:
    for event in pygame.event.get():
        if event.type == QUIT:
            run = False
            break
#            sys.exit()

    # Mechanism
    for atom1 in atoms:
        for atom2 in [x for x in atoms if x is not atom1]:
            atom1.SoftS(atom2, e=1E6)
        atom1.step(grid, dt=0.001)

    # Reset screen
    screen.fill([0,0,0])
    screen.blit(cells_bg, [0,0])

    # Draw
    for atom in atoms:
        draw_atom(screen, atom, [255,0,0])

    # Update
    pygame.display.flip()
    fpsClock.tick(fps)

pygame.quit()
