#!/usr/bin/env python3
# -*- coding: iso-8859-15 -*-

import numpy as np
from sys import argv


infile = argv[1]
outfile = argv[2]
dt = float(argv[3])

# Create main data file
data = np.load(infile)

# Create .xyz file
xyz_filename = '{}.xyz'.format(outfile)
with open(xyz_filename, 'w') as f:
    for s, step in enumerate(data):
        f.write('{}\n'.format(len(step)))
        f.write('Step = {:05d}, t={:0.5f}\n'.format(s, s*dt))
        for atom in step:
            f.write('He {}\n'.format(' '.join(map(str, atom[:3]))))

# Create kinetic energy and momentum file
# NOTE: Assuming mass=1 for all atoms
measure_filename = '{}.data'.format(outfile)
with open(measure_filename, 'w') as f:
    f.write('Ek P\n')
    for step in data:
        total_ek = np.sum([0.5*np.linalg.norm(v)**2 for v in step[:,3:]])
        cm_momentum = np.sum([v for v in step[:,3:]], axis=0)
        f.write('{} {}\n'.format(
            total_ek,
            ' '.join(map(str, cm_momentum))
        ))
