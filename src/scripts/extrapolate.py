import numpy as np
import scipp as sc
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from kinisi.analyze import ConductivityAnalyzer
from kinisi.arrhenius import arrhenius
from matplotlib import rcParams

from _fig_params import CREDIBLE_INTERVALS, ALPHABET, PAGE_WIDTH, COLUMN_WIDTH, GOLDEN_RATIO
import _fig_params as fp
import paths

rcParams.update({"axes.spines.top": True, "axes.spines.right": True})

long = {'mean': [], 'var': []}
temp = np.arange(500, 900, 100)
for i in temp:
    d = ConductivityAnalyzer.from_hdf5(paths.data / f"reduced/llzo/t{i}/diffusion.h5")
    sigmaT = sc.to_unit(d.sigma, 'mS/cm') * (i * sc.Unit('K'))
    long['mean'].append(sc.mean(sigmaT).value)
    long['var'].append(sc.var(sigmaT, ddof=1).value)
    unit = sc.mean(sigmaT).unit

samples = np.load(paths.data / "reduced/llzo/arrhenius_samples.npz")['samples']
extrapolated_temperature = np.load(paths.data / "reduced/llzo/arrhenius_samples.npz")['extrapolated_temperature']

print(unit)

figsize = (COLUMN_WIDTH, COLUMN_WIDTH/GOLDEN_RATIO)
fig = plt.figure(figsize=figsize)
gs = gridspec.GridSpec(1, 1, figure=fig, wspace=0.2, hspace=0.5)

axes = []
labels = []

axes.append(fig.add_subplot(gs[0:2, :]))
axes[-1].errorbar(1 / temp * 1000,
                  long['mean'],
                  np.sqrt(long['var']),
                  marker=".",
                  ls="",
                  color='k',
                  label=r'$\sigma_{\mathrm{Li}^+}(T)T$')
temp_smooth = np.arange(1000 / 3.8, 1000 / 1.15, 1)
k = arrhenius(temp_smooth[:, np.newaxis], *samples)
axes[-1].fill_between(1 / temp_smooth * 1000,
                      *np.percentile(k, CREDIBLE_INTERVALS[0][:-1], axis=1),
                      color=fp.colors[0],
                      lw=0,
                      label=r'$p[\mathbf{\theta}_{\alpha}\,|\,\sigma_{\mathrm{Li^+}}(T)]$')
for j in CREDIBLE_INTERVALS[1:]:
    axes[-1].fill_between(1 / temp_smooth * 1000,
                          *np.percentile(k, j[:-1], axis=1),
                          alpha=j[-1],
                          color=fp.colors[0],
                          lw=0)

axes[-1].axvline(1000 / 300, color='k', ls='--', lw=1, ymin=0.15, ymax=0.5)
axes[-1].set_yscale("log", base=np.e)
axes[-1].set_ylim(np.e ** np.log(4000), np.e ** np.log(700000))
axes[-1].set_yticks([5000, 50000, 500000])
axes[-1].set_yticklabels([f"{np.e ** np.log(val / 1e4):.1f}" for val in [5000, 50000, 500000]])
axes[-1].set_xlabel("$T^{-1}$ / K$^{-1}$")
axes[-1].set_ylabel(r"$\sigma_{\mathrm{Li}^+}(T)T$ / $10^4$ mScm$^{-1}$K")
axes[-1].set_xlim(1000 / 800 - 0.075, 1000 / 275 + 0.075)
# axes[-1].set_xticks([1.2, 1.4, 1.6, 1.8, 2.0])
ax2 = axes[-1].secondary_xaxis('top', functions=(lambda x: 1000 / x, lambda x: 1000 / x))
ax2.set_xlabel("$T$ / 10$^{-3}$ K")
ax2.set_xticks([700, 600, 500, 400, 300])
labels.append(ALPHABET[0])

axes[-1].legend(ncol=1, loc='upper right')

rcParams.update({"axes.spines.top": False, "axes.spines.right": False, "axes.spines.left": False})

y, x = np.histogram(extrapolated_temperature / 300, density=True, bins=fp.NBINS)
inset_ax = fig.add_axes([0.17, 0.3, 0.3, 0.3])
inset_ax.stairs(y, x, fill=True, color=fp.colors[0])
# axes[-1].set_ylabel(r'$p[E_a\,|\,\sigma_{\mathrm{Li}^+}(T)]$ / eV$^{-1}$')
inset_ax.set_yticks([])
inset_ax.set_xlabel(r"$\sigma_{\mathrm{Li}^+}(300$ K$)$ / mScm$^{-1}$")
inset_ax.patch.set_alpha(0.0) 
labels.append(ALPHABET[1])

# axes.append(fig.add_subplot(gs[2, 1]))

# y, x = np.histogram(samples[1] / 1e5, density=True, bins=fp.NBINS)
# axes[-1].stairs(y, x, fill=True, color=fp.colors[0])
# # axes[-1].set_ylabel(r"$p[A\,|\,\sigma_{\mathrm{Li}^+}(T)]$ / $10^7$ mS$^{-1}$cmK$^{-1}$")
# axes[-1].set_xlabel(r"$A\,|\,\sigma_{\mathrm{Li}^+}(T)$ / $10^5$ mScm$^{-1}$K")
# axes[-1].set_yticks([])
# axes[-1].set_xticks([20, 30, 40])
# labels.append(ALPHABET[2])

# fig.align_ylabels(axes)

for i, ax in enumerate(axes):
    print(ax.get_window_extent().x1 - ax.get_window_extent().x0)
    print(ax.get_window_extent().y1 - ax.get_window_extent().y0)

fig.set_size_inches(*figsize)
plt.savefig(paths.figures / "extrapolate.pdf", bbox_inches="tight", transparent=True)
plt.close()

