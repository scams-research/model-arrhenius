# Copyright (c) Andrew R. McCluskey
# Distributed under the terms of the MIT License
# author: Andrew R. McCluskey (arm61)

import numpy as np
import scipp as sc
from kinisi.arrhenius import Arrhenius, VogelFulcherTammann
from kinisi.analyze import ConductivityAnalyzer

length = snakemake.params['length']
seed = snakemake.params['seed']

np.random.seed(seed)

temp = np.array([300, 350, 400, 500, 600, 700])
D = {'mean': [], 'var': []}
for i in temp:
    d = ConductivityAnalyzer.from_hdf5(f'src/data/reduced/agcrse2/Dmsd{i}_{length}.h5')
    sigmaT = sc.to_unit(d.sigma, 'mS/cm') * (i * sc.Unit('K'))
    D['mean'].append(sc.mean(sigmaT).value)
    D['var'].append(sc.var(sigmaT, ddof=1).value)
    unit = sc.mean(sigmaT).unit

td = sc.DataArray(
    data=sc.array(dims=['temperature'],
                  values=D['mean'],
                  variances=D['var'],
                  unit=unit),
    coords={
        'temperature': sc.Variable(dims=['temperature'],
                                   values=temp,
                                   unit='K')
    })

s = Arrhenius(td, bounds=((0 * sc.Unit('eV'), 1 * sc.Unit('eV')), (1e4 * td.data.unit, 1e8 * td.data.unit)))
s.nested_sampling()
s_samples = np.array([i.values for i in s.flatchain.values()])

t = VogelFulcherTammann(td, 
                        bounds=((0 * sc.Unit('eV'), 1 * sc.Unit('eV')), 
                                (1e4 * td.data.unit, 1e8 * td.data.unit), 
                                (0 * sc.Unit('K'), 300 * sc.Unit('K'))))
t.nested_sampling()
t_samples = np.array([i.values for i in t.flatchain.values()])

bf = 2 * (t.data_group['logz'] - s.data_group['logz'])
np.savez(f'src/data/reduced/agcrse2/bayes_{length}_{seed}.npz',
         s_samples=s_samples,
         ev=np.array([s.data_group['logz'].value, s.data_group['logz'].variance]),
         t_samples=t_samples,
         sev=np.array([t.data_group['logz'].value, t.data_group['logz'].variance]),
         bf=np.array([bf.value, bf.variance]))