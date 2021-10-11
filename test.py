#!/usr/bin/env python3
# -*- coding: iso-8859-15 -*-

import numpy as np
from md import Grid
import sys
import pygame
from pygame.locals import *


def draw_cell(surface, grid, index, fill=False):
    ulc = grid.get_coordinates(index)[0]
    dims = grid.L/grid.n
    rect = (ulc[0], ulc[1], dims[0], dims[1])
    if fill is not False:
        pygame.draw.rect(surface, fill, rect)
    pygame.draw.rect(surface, [255, 255, 255], rect, width=2)


width, height = np.random.randint(200,800,2)
grid = Grid(n=np.random.randint(5,15), L=[width,height,0])

pygame.init()

fps = 60
fpsClock = pygame.time.Clock()

screen = pygame.display.set_mode((width, height))

for idx in grid.indices:
    draw_cell(screen, grid, idx)
idx_fill = (0,0)

# Game loop
while True:
    for event in pygame.event.get():
        if event.type == QUIT:
            pygame.quit()
            sys.exit()

    # Mechanism
    prev_idx_fill = idx_fill
    cursor_pos = pygame.mouse.get_pos()
    idx_fill = tuple(grid.get_indices(cursor_pos))

    # Draw
    if idx_fill != prev_idx_fill:
        draw_cell(screen, grid, prev_idx_fill, [0,0,0])
    else:
        draw_cell(screen, grid, idx_fill, [0,255,0])
    pygame.draw.circle(screen, [255,0,0], cursor_pos, 3)

    # Update
    pygame.display.flip()
    fpsClock.tick(fps)
