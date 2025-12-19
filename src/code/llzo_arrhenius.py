# Copyright (c) Andrew R. McCluskey
# Distributed under the terms of the MIT License
# author: Andrew R. McCluskey (arm61)

import numpy as np
import scipp as sc
from scipp.constants import k
from kinisi.analyze import ConductivityAnalyzer
from kinisi.arrhenius import Arrhenius
from kinisi.samples import Samples

temp = np.arange(500, 900, 100)
D = {'mean': [], 'var': []}
for i in temp:
    d = ConductivityAnalyzer.from_hdf5(f"src/data/reduced/llzo/t{i}/diffusion.h5")
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

s = Arrhenius(td, bounds=((0 * sc.Unit('eV'), 0.5 * sc.Unit('eV')), (1e5 * td.data.unit, 1e7 * td.data.unit)))
s.mcmc()
samples = np.array([i.values for i in s.flatchain.values()])

extra_samples = s.extrapolate(extrapolated_temperature=300 * sc.Unit('K'))
np.savez('src/data/reduced/llzo/arrhenius_samples.npz', samples=samples, extrapolated_temperature=extra_samples.values)
