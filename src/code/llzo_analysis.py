# Copyright (c) Andrew R. McCluskey
# Distributed under the terms of the MIT License

# author: Andrew R. McCluskey (arm61)

import MDAnalysis as mda
import numpy as np
import scipp as sc
from kinisi.analyze import ConductivityAnalyzer
from tqdm import tqdm


def universe_slice(u: mda.Universe) -> mda.Universe:
    """
    Renames atom types and sets dimensions for LLZO data.

    :param u: Input Universe
    :returns: Processed Universe
    """
    for t in u.trajectory:
        for i in u.atoms:
            if i.type == 'LI':
                i.type = 'Li'
            if i.type == 'Li1':
                i.type = 'Li'
            if i.type == 'Li2':
                i.type = 'Li'
            if i.type == 'L':
                i.type = 'La'
            if i.type == 'Z':
                i.type = 'Zr'
        t.dimensions = np.array([26.175106662244392197, 26.175106662244392197, 26.175106662244392197, 90, 90, 90])
    return u


step_skip = 100  # sampling rate
timestep = 5.079648e-4 * sc.Unit('ps')

temp = snakemake.params['temp']

ul = mda.Universe(f'src/data/raw/llzo/t{temp}/traj_{temp}.gro', f'src/data/raw/llzo/t{temp}/traj_{temp}.xtc')
u = universe_slice(ul)
da_params = {'specie': 'Li', 'time_step': timestep, 'step_skip': step_skip}

new_dt = sc.arange(dim='time interval', start=timestep * step_skip,
                   stop=timestep * step_skip * len(u.trajectory),
                   step=timestep * step_skip * 5)
da_params['dt'] = new_dt
da_params['system_particles'] = len(u.select_atoms('type Li'))
da_params['ionic_charge'] = 1 * sc.Unit('e')

rng = np.random.RandomState(0)
np.random.seed(42)

d = ConductivityAnalyzer.from_universe(u, **da_params)
d.conductivity(10 * sc.Unit('ps'), temperature=temp * sc.Unit('K'))
d.to_hdf5(f'src/data/reduced/llzo/t{temp}/diffusion.h5')