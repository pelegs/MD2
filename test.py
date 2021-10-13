#!/usr/bin/env python3
# -*- coding: iso-8859-15 -*-

import numpy as np
from md import Grid, get_neighboring_indices
import sys
import pygame
from pygame.locals import *


def draw_cell(surface, cell, fill=[0,0,0]):
    rect = (cell.ul_corner[0], cell.ul_corner[1],
            cell.L[0], cell.L[1])
    pygame.draw.rect(surface, fill, rect)
    pygame.draw.rect(surface, [255, 255, 255], rect, width=2)


width, height = 800, 800
grid = Grid(n=15, L=[width,height,800], neighbors_dist=2)

pygame.init()

fps = 60
fpsClock = pygame.time.Clock()

screen = pygame.display.set_mode((width, height))

for index1d, cell in enumerate(grid.cells):
    draw_cell(screen, cell, [0,0,0])
cursor_idx_1d = 0
current_cell = grid.cells[0]

# Game loop
while True:
    for event in pygame.event.get():
        if event.type == QUIT:
            pygame.quit()
            sys.exit()

    # Mechanism
    cursor_idx_1d_prev = cursor_idx_1d
    prev_cell = grid.cells[cursor_idx_1d_prev]
    cursor_pos = pygame.mouse.get_pos()
    cursor_pos_3d = (cursor_pos[0], cursor_pos[1], 400)
    cursor_idx_3d = tuple(grid.get_index_3d_from_pos(cursor_pos_3d))
    cursor_idx_1d = grid.get_index_1d_from_3d(cursor_idx_3d)
    current_cell = grid.cells[cursor_idx_1d]

    # Draw
    if cursor_idx_1d != cursor_idx_1d_prev:
        draw_cell(screen, prev_cell, [0,0,0])
        prev_cell = current_cell
        for neighbor in grid.cells[cursor_idx_1d_prev].neighbors:
            draw_cell(screen, neighbor, [0,0,0])
    for neighbor in grid.cells[cursor_idx_1d].neighbors:
        draw_cell(screen, neighbor, [0,150,0])
    draw_cell(screen, current_cell, [0,255,0])
    pygame.draw.circle(screen, [255,0,0], cursor_pos, 3)

    # Update
    pygame.display.flip()
    fpsClock.tick(fps)
