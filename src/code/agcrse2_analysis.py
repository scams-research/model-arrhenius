# Copyright (c) Andrew R. McCluskey
# Distributed under the terms of the MIT License
# author: Andrew R. McCluskey (arm61)

import numpy as np
import scipp as sc
import MDAnalysis as mda
from kinisi.analyze import ConductivityAnalyzer
from mda_helper import find_step_skip, flatten_universe_list, universe_slice

temperature = snakemake.params['temp']
length = snakemake.params['length']

ngp = {300: 70, 350: 50, 400: 30, 500: 25, 600: 22.5, 700: 20}

dt_skip = ngp[temperature]
files = [f'src/data/raw/agcrse2/agcrse2_T{temperature}_{i}.lammpstrj' for i in range(1, 22)]

step_skip = find_step_skip(files[0])
u_list = [mda.Universe(f, topology_format='LAMMPSDUMP', atom_style='id type x y z') for f in files]
u = flatten_universe_list(u_list, 'src/data/raw/agcrse2/structure.data')

max_slice = int((dt_skip + length) / (step_skip * 0.001))
if max_slice > len(u.trajectory):
    raise ValueError('Not enough data.')

short_u = universe_slice(u, 'src/data/raw/agcrse2/structure.data', max_slice)

rng = np.random.RandomState(42)
np.random.seed(42)

d = ConductivityAnalyzer.from_universe(short_u,
                                       specie = None,
                                       species_indices = sc.arange(dim='indices', start=384, stop=576, step=1),
                                       time_step = 0.001 * sc.Unit('ps'),
                                       step_skip = step_skip,
                                       system_particles = len(short_u.select_atoms('type 2')),
                                       ionic_charge = 1 * sc.Unit('e'),
                                       dt=sc.arange(dim='time interval',
                                                          start=0.001 * step_skip * sc.Unit('ps'),
                                                          stop=0.001 * step_skip * sc.Unit('ps') * len(short_u.trajectory),
                                                          step=0.001 * step_skip * sc.Unit('ps') * 5)
                                       )
d.conductivity(dt_skip * sc.Unit('ps'), temperature=temperature * sc.Unit('K'))

d.to_hdf5(f'src/data/reduced/agcrse2/Dmsd{temperature}_{length}.h5')