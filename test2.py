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



width, height = 800, 800
grid = Grid(n=15, L=[width,height,800], neighbors_dist=2)

pygame.init()

fps = 60
fpsClock = pygame.time.Clock()

screen = pygame.display.set_mode((width, height))
cells_bg = pygame.Surface((width, height))

# Grid
for index1d, cell in enumerate(grid.cells):
    draw_cell(cells_bg, cell, [0,0,0])

# Particles
atom1 = Atom(
    pos=np.array([400,400,400]),
    vel=np.array([-250,0,0]),
    rad=20,
    mass=1E3,
)
atom2 = Atom(
    pos=np.array([300,400,400]),
    vel=np.array([0,0,0]),
    rad=7,
    mass=1,
)

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
    E=1E5
    ff = atom1.SoftS(atom2, e=E)
    forces.append((dist(atom1.pos, atom2.pos), atom1.F[0]))
    F1 = set_norm(atom1.F, 30)
    atom2.SoftS(atom1, e=E)
    F2 = set_norm(atom2.F, 30)
    atom1.step(dt=0.001)
    atom2.step(dt=0.001)

    # Reset screen
    screen.fill([0,0,0])
    screen.blit(cells_bg, [0,0])

    # Draw
    draw_atom(screen, atom1, [255,0,0])
    draw_vec(
        screen,
        atom1.pos[:2],
        F1[:2],
        [255,255,255]
    )
    draw_atom(screen, atom2, [0,255,0])
    draw_vec(
        screen,
        atom2.pos[:2],
        F2[:2],
        [255,255,255]
    )

    # Update
    pygame.display.flip()
    fpsClock.tick(fps)

pygame.quit()

# Data?
with open('forces.data', 'w') as file:
    for p in forces:
        file.write('{}\n'.format(' '.join(map(str, p))))
