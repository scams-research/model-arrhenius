import sys
import itertools
import numpy as np
import scipp as sc
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from kinisi.analyze import ConductivityAnalyzer

from _fig_params import ALPHABET, CREDIBLE_INTERVALS, mid_points, PAGE_WIDTH, GOLDEN_RATIO
import _fig_params as fp
import paths

D = []
temp = np.array([500, 600, 700, 800])
for i in temp:
    d = ConductivityAnalyzer.from_hdf5(paths.data / f"reduced/llzo/t{i}/diffusion.h5")
    D.append(d)

size_adapt = 1
figsize = (8.06, 2)
label_offset = [0, 1 + 0.1 * size_adapt]
fig = plt.figure(figsize=figsize)
gs = gridspec.GridSpec(1, 4, figure=fig, wspace=0.6)

axes = []
axes.append(fig.add_subplot(gs[0, 0]))
axes.append(fig.add_subplot(gs[0, 1]))
axes.append(fig.add_subplot(gs[0, 2]))
axes.append(fig.add_subplot(gs[0, 3]))

for i, ax in enumerate(axes):
    line1 = ax.errorbar(D[i].dt.values, D[i].mscd.values, D[i].mscd.variances ** 0.5, color=fp.colors[0])
    line2 = ax.fill_between(D[i].dt.values,
                            *np.percentile(D[i].distributions, CREDIBLE_INTERVALS[0][:-1], axis=1),
                            alpha=CREDIBLE_INTERVALS[0][-1],
                            color=fp.colors[2],
                            lw=0,
                            zorder=10)
    for j in CREDIBLE_INTERVALS[1:]:
        ax.fill_between(D[i].dt.values,
                        *np.percentile(D[i].distributions, j[:-1], axis=1),
                        alpha=j[-1],
                        color=fp.colors[2],
                        lw=0,
                        zorder=10)
    ax.set_xlabel(r"$\Delta t$ / ps")
    ax.set_ylabel(r"$\langle\Delta\mathbf{r}_{\mathrm{c}}(t)^2_{\mathrm{Ag}^+}\rangle$ / Å$^2$")
    ax.set_xlim(0, None)
    ax.set_ylim(0, 35)
    ax.set_yticks([0, 10, 20, 30])
    ax.text(1, 0.05, f'$T={temp[i]}$ K', horizontalalignment='right', verticalalignment='bottom', transform=ax.transAxes)

plt.savefig(paths.figures / f"llzo_plots.pdf", transparent=True)
plt.close()