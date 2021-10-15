#!/usr/bin/env python3
# -*- coding: iso-8859-15 -*-

import numpy as np
from numpy.linalg import norm
from itertools import product


# ----------------- General maths functions  ----------------- #

def wrap(x, L):
    """
    Wraps x using the periodic boundry [0,L).
    Example: given x=3 with L=2, the return
    value is w=1. Given x=-1 with L=4, the
    return value is w=3.
    """
    return x-L*np.floor(x/L)

def get_neighboring_indices(indices, n=10, step=1):
    i, j, k = indices
    x_ind = [x%n for x in np.arange(i-step, i+step+1, 1)]
    y_ind = [y%n for y in np.arange(j-step, j+step+1, 1)]
    z_ind = [z%n for z in np.arange(k-step, k+step+1, 1)]
    return list(product(x_ind, y_ind, z_ind))


# ----------------- Vector related functions ----------------- #

def norm_sqr(v):
    return np.sum([x**2 for x in v])

def dist_sqr(v1, v2):
    return norm_sqr(v1-v2)

def dist(v1, v2):
    return norm(v1-v2)

def normalize(v):
    n = norm(v)
    if n != 0:
        return v/n
    else:
        return np.zeros(3)

def set_norm(v, n):
    return n*normalize(v)

def look_at(v1, v2):
    return normalize(v2-v1)


# ----------------- Atom class  ----------------- #

class Atom:
    def __init__(self,
                 pos=np.zeros(3), vel=np.zeros(3),
                 mass=1, rad=1, id=-1, color=[255,255,255]):
        self.pos = pos
        self.vel = vel
        self.mass = mass
        self.m_ = mass ** -1.0
        self.rad = rad
        self.id = id
        self.color = color
        self.a = np.zeros(3)
        self.a_next = np.zeros(3)
        self.reset_force()

    def reset_force(self):
        self.F = np.zeros(3)

    def sum_forces(self):
        return np.sum(self.F, axis=0)

    def add_force(self, F):
        self.F = np.vstack((self.F, F))

    def SoftS(self, atom2, e=1.0, force=None):
        if force is not None:
            self.add_force(-force)
            return
        else:
            d2 = dist_sqr(self.pos, atom2.pos)
            d6 = d2**3
            d13 = d2**6.5
            s6 = (self.rad+atom2.rad)**6
            s12 = s6**2
            dir = look_at(self.pos, atom2.pos)
            SF = 24 * e * s6 * (d6-2*s6) / d13
            self.add_force(SF*dir)
        return SF

    def step(self, grid, dt=0.0001):
        F = self.sum_forces()
        self.a = F*self.m_
        self.reset_force()
        self.vel = self.vel + self.a*dt
        self.pos = self.pos + self.vel*dt

        # Boundry conditions
        for i in range(3):
            if not (0 <= self.pos[i] <= grid.L[i]):
                self.pos[i] = wrap(self.pos[i], grid.L[i])

    def step1(self, dt=0.001):
        self.pos = self.pos + self.vel*dt + 0.5*self.a*dt**2

    def step2(self, dt=0.001):
        self.a = self.a_next
        self.a_next = self.F*self.m_
        self.reset_force()
        self.vel = self.vel + 0.5*(self.a+self.a_next)*dt

    def calc_P(self):
        self.P = self.m * self.vel

    def calc_KE(self):
        self.KE = 0.5 * self.m * norm_sqr(self.vel)


# ----------------- Cell class  ----------------- #

class Cell:
    def __init__(self, index1d, index3d, ul_corner, L=np.ones(3)):
        self.index1d = index1d
        self.index3d = index3d
        self.ul_corner = ul_corner
        self.L = L
        self.neighbors = []
        self.particles = []

    def clear(self):
        self.particles.clear()

    def insert(self, obj):
        self.particles.append(obj)

    def remove(self, obj):
        if obj in self.particles:
            self.particles.remove(obj)

    def set_neighbors(self, neighbor_cells):
        self.neighbors = neighbor_cells

    def __str__(self):
        return 'idx1d: {}, idx3d: {}, corners: {}'.format(
            self.index1d, self.index3d, self.corners
        )


# ----------------- Grid class  ----------------- #

class Grid:
    def __init__(self, n=10, L=np.ones(3), neighbors_dist=1):
        self.n = n
        self.L = np.array(L)
        self.cell_dim = self.L/n
        self.bins = np.linspace((0,0,0), self.L, self.n+1).T
        self.indices = [
            x for x in product(range(n), range(n), range(n))
        ]
        self.neighbors_dist = neighbors_dist

        # Create cells
        self.cells = []
        for index1d, index3d in enumerate(self.indices):
            corners = self.get_cell_ul_corner(index3d)
            new_cell = Cell(index1d, index3d, corners, self.L/self.n)
            self.cells.append(new_cell)

        # Set neighbors for each cell
        for cell in self.cells:
            neighbors = [self.cells[self.get_index_1d_from_3d(index3d)]
                         for index3d in get_neighboring_indices(
                            cell.index3d, self.n, step=self.neighbors_dist
                        )]
            cell.set_neighbors(neighbors)

    def reset(self):
        for cell in self.cells:
            cell.reset()

    def insert(self, obj):
        # Check validity of index
        for i in range(3):
            if not (0 <= index[i] < self.n):
                raise ValueError('Index {} out of range.'.format(i))
                return

        # Add obj to cell
        index3d = self.get_index_3d_from_pos(obj.pos)
        index1d = self.get_index_1d_from_3d(index3d)
        self.cells[index1d].insert(obj)

    def get_index_3d_from_pos(self, pos):
        """
        Returns the cell index for the position `pos`.
        NOTE: periodic boundry conditions.
        """
        pos_wrapped = [wrap(x, L) for x, L in zip(pos, self.L)]
        index = [np.digitize(x, b) for x, b in zip(pos_wrapped, self.bins)]
        # subtract 1 from each index because np starts counting from 1
        index = [n-1 for n in index]
        return index

    def get_index_3d_from_1d(self, idx1d):
        return self.indices[idx1d]

    def get_index_1d_from_3d(self, idx3d):
        i, j, k = idx3d
        return i*self.n**2 + j*self.n + k - 1

    def get_cell_ul_corner(self, index3d):
        i, j, k = index3d
        return np.array([self.bins[0,i], self.bins[1,j], self.bins[2,k]])


# ----------------- Main  ----------------- #

if __name__ == '__main__':
    grid = Grid(n=10, L=[800,600,400])
    for cell in grid.cells[0].neighbors:
        print(cell.index3d)
