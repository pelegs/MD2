#!/usr/bin/env python3
# -*- coding: iso-8859-15 -*-

import numpy as np
from numpy.linalg import norm

# ---------------------------------------- #

def norm_sqr(v):
    return np.sum([x**2 for x in v])

def dist_sqrt(v1, v2):
    return norm_sqr(v1-v2)

def dist(v1, v2):
    return norm(v1-v2)

# ---------------------------------------- #

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

# ---------------------------------------- #

if __name__ == '__main__':
    print('hi')
