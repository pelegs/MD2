#!/usr/bin/env python3
# -*- coding: iso-8859-15 -*-

import numpy as np
from numpy.linalg import norm


# -------- General maths functions  ------- #

def wrap(x, L):
    """
    Wraps x using the periodic boundry [0,L).
    Example: given x=3 with L=2, the return
    value is w=1. Given x=-1 with L=4, the
    return value is w=3.
    """
    return x-L*np.floor(x/L)


# -------- Vector related functions ------- #

def norm_sqr(v):
    return np.sum([x**2 for x in v])

def dist_sqrt(v1, v2):
    return norm_sqr(v1-v2)

def dist(v1, v2):
    return norm(v1-v2)


# -------- Atom class  ------- #

class Atom:
    def __init__(self,
                 pos=np.zeros(3), vel=np.zeros(3),
                 mass=1, rad=1):
        self.pos = pos
        self.vel = vel
        self.mass = mass
        self.m_ = mass ** -1.0
        self.rad = rad
        self.a = np.zeros(3)
        self.a_next = np.zeros(3)
        self.reset_force()

    def reset_force(self):
        self.F = np.zeros(3)

    def add_force(self, F):
        self.F = self.F + F

    def LJ(self, atom2, eps=1.0, force=None):
        if force:
            self.add_force(-force)
            return
        else:
            d = dist(self.pos, atom2.pos)
            F = -24*eps*self.rad**6*(d**6-2*self.rad**6)/d**13
            self.add_force(F)

    def step1(dt=0.001):
        self.pos = self.pos + self.vel*dt + 0.5*self.a*dt**2

    def step2(dt=0.001):
        self.a = self.a_next
        self.a_next = self.F*self.m_
        self.reset_force()
        self.vel = self.vel + 0.5*(self.a+self.a_next)*dt

    def calc_P(self):
        self.P = self.m * self.vel

    def calc_KE(self):
        self.KE = 0.5 * self.m * norm_sqr(self.vel)


# -------- Grid class  ------- #

class Grid:
    def __init__(self, n=10, L=np.ones(3)):
        self.n = n
        self.L = np.array(L)
        self.cell_dim = self.L/n
        self.bins = np.linspace((0,0,0), self.L, self.n+1).T
        self.reset()
        self.indices = [(x,y) for x in range(n) for y in range(n)]

    def reset(self):
        self.cells = [[[[]
                      for _ in range(self.n)]
                      for _ in range(self.n)]
                      for _ in range(self.n)
        ]

    def insert(obj, index=np.zeros(3)):
        # Check validity of index
        for i in range(3):
            if not (0 <= index[i] < self.n):
                raise ValueError('Index {} out of range.'.format(i))
                return

        # Add obj to cell
        ix, iy, iz = self.get_indices(obj.pos)
        self.cells[ix][iy][iz] = obj

    def get_indices(self, pos):
        """
        Returns the cell index for the position `pos`.
        NOTE: periodic boundry conditions.
        """
        pos_wrapped = [wrap(x, L) for x, L in zip(pos, self.L)]
        indices = [np.digitize(x, b) for x, b in zip(pos_wrapped, self.bins)]
        # subtract 1 from each index because np starts counting from 1
        indices = [n-1 for n in indices]
        return indices

    def get_coordinates(self, index):
        i, j = index
        c1 = (self.bins[0,i], self.bins[1,j])
        c2 = (self.bins[0,i+1], self.bins[1,j+1])
        return np.array([c1, c2])


# -------- Main  ------- #

if __name__ == '__main__':
    grid = Grid(n=10, L=[800,600,0])
