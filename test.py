#!/usr/bin/env python3
# -*- coding: iso-8859-15 -*-

import numpy as np
from md import Grid, get_neighboring_indices
import sys
import pygame
from pygame.locals import *


def draw_cell(surface, grid, index, fill=[255,0,0]):
    ulc = grid.get_coordinates(index)[0]
    dims = grid.L/grid.n
    rect = (ulc[0], ulc[1], dims[0], dims[1])
    pygame.draw.rect(surface, fill, rect)
    pygame.draw.rect(surface, [255, 255, 255], rect, width=2)


width, height = 800, 800
grid = Grid(n=15, L=[width,height,0])

pygame.init()

fps = 60
fpsClock = pygame.time.Clock()

screen = pygame.display.set_mode((width, height))


# Game loop
while True:
    for event in pygame.event.get():
        if event.type == QUIT:
            pygame.quit()
            sys.exit()

    # Mechanism
    cursor_pos = pygame.mouse.get_pos()
    cursor_idx = tuple(grid.get_indices(cursor_pos))
    neighbor_idx = get_neighboring_indices(cursor_idx, grid.n, 1, dim=2)

    # Draw
    for idx in grid.indices:
        if idx in neighbor_idx:
            if tuple(idx) == tuple(cursor_idx):
                draw_cell(screen, grid, idx, [0,255,0])
            else:
                draw_cell(screen, grid, idx, [0,150,0])
        else:
            draw_cell(screen, grid, idx, [0,0,0])
    pygame.draw.circle(screen, [255,0,0], cursor_pos, 3)

    # Update
    pygame.display.flip()
    fpsClock.tick(fps)
