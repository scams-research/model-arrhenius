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

print(unit)

figsize = (COLUMN_WIDTH, COLUMN_WIDTH/GOLDEN_RATIO * 3)
fig = plt.figure(figsize=figsize)
gs = gridspec.GridSpec(3, 2, figure=fig, wspace=0.2, hspace=0.4, height_ratios=[1, 0.8, 0.6])

axes = []
labels = []

axes.append(fig.add_subplot(gs[0, :]))
axes[-1].errorbar(1 / temp * 1000,
                  long['mean'],
                  np.sqrt(long['var']),
                  marker=".",
                  ls="",
                  color='k',
                  label=r'$\sigma_{\mathrm{Li}^+}T$')
temp_smooth = np.arange(1000 / 2.1, 1000 / 1.15, 1)
k = arrhenius(temp_smooth[:, np.newaxis], *samples)
axes[-1].fill_between(1 / temp_smooth * 1000,
                      *np.percentile(k, CREDIBLE_INTERVALS[0][:-1], axis=1),
                      color=fp.colors[0],
                      lw=0,
                      label=r'$p[\mathbf{\theta}\,|\,\sigma_{\mathrm{Li^+}}]$')
for j in CREDIBLE_INTERVALS[1:]:
    axes[-1].fill_between(1 / temp_smooth * 1000,
                          *np.percentile(k, j[:-1], axis=1),
                          alpha=j[-1],
                          color=fp.colors[0],
                          lw=0)

axes[-1].set_yscale("log", base=np.e)
axes[-1].set_ylim(np.e ** np.log(130000), np.e ** np.log(520000))
axes[-1].set_yticks([200000, 300000, 400000, 500000])
axes[-1].set_yticklabels([f"{np.e ** np.log(val / 1e5):.0f}" for val in [200000, 300000, 400000, 500000]])
axes[-1].set_xlabel("$T^{-1}$ / 10$^{-3}$ K$^{-1}$")
axes[-1].set_ylabel(r"$\sigma_{\mathrm{Li}^+}T$ / $10^5$ mScm$^{-1}$K")
axes[-1].set_xlim(1000 / 800 - 0.075, 1000 / 500 + 0.075)
axes[-1].set_xticks([1.2, 1.4, 1.6, 1.8, 2.0])
ax2 = axes[-1].secondary_xaxis('top', functions=(lambda x: 1000 / x, lambda x: 1000 / x))
ax2.set_xlabel("$T$ / K")
ax2.set_xticks(temp)
labels.append(ALPHABET[0])

axes[-1].legend(ncol=1, loc='upper right', frameon=False)

rcParams.update({"axes.spines.top": False, "axes.spines.right": False, "axes.spines.left": False})

axes.append(fig.add_subplot(gs[2, 0]))

y, x = np.histogram(samples[0], density=True, bins=fp.NBINS)
axes[-1].stairs(y, x, fill=True, color=fp.colors[0])
# axes[-1].set_ylabel(r'$p[E_a\,|\,\sigma_{\mathrm{Li}^+}]$ / eV$^{-1}$')
axes[-1].set_xticks([0.11, 0.14])
axes[-1].set_xlabel(r'$E_a\,|\,\sigma_{\mathrm{Li}^+}T$ / eV')
axes[-1].set_yticks([])
# axes[-1].set_aspect('equal')
labels.append(ALPHABET[2])

axes.append(fig.add_subplot(gs[2, 1]))

y, x = np.histogram(samples[1] / 1e6, density=True, bins=fp.NBINS)
axes[-1].stairs(y, x, fill=True, color=fp.colors[0])
# axes[-1].set_ylabel(r"$p[A\,|\,\sigma_{\mathrm{Li}^+}]$ / $10^7$ mS$^{-1}$cmK$^{-1}$")
axes[-1].set_xlabel(r"$A\,|\,\sigma_{\mathrm{Li}^+}T$ / $10^6$ mScm$^{-1}$K")
axes[-1].set_xticks([2, 3, 4])
# axes[-1].set_aspect('equal')
axes[-1].set_yticks([])
labels.append(ALPHABET[3])

rcParams.update({"axes.spines.top": False, "axes.spines.right": False, "axes.spines.bottom": True, "axes.spines.left": True})

axes.append(fig.add_subplot(gs[1, :]))

H, xedges, yedges = np.histogram2d(samples[0], samples[1] / 1e6, density=True, 
                                   bins=(np.linspace(0.105, 0.145, fp.NBINS), np.linspace(1.9, 3.9, fp.NBINS)))
axes[-1].contourf(xedges[:-1] + np.diff(xedges), yedges[:-1] + np.diff(yedges), H.T, levels=20, cmap='Blues')
axes[-1].set_box_aspect(1)
axes[-1].text(-0.23, 0.5, r"$A\,|\,\sigma_{\mathrm{Li}^+}T$ / $10^6$ mScm$^{-1}$K", va='center', rotation=90, transform=axes[-1].transAxes)
axes[-1].set_xticks([0.12, 0.14])
axes[-1].set_yticks([2, 3])
axes[-1].set_xlabel(r'$E_a\,|\,\sigma_{\mathrm{Li}^+}T$ / eV')
labels.append(ALPHABET[1])

fig.align_ylabels(axes)

for i, ax in enumerate(axes[-2:]):
    print(ax.get_window_extent().x1 - ax.get_window_extent().x0)
    print(ax.get_window_extent().y1 - ax.get_window_extent().y0)

x_correction = [5, 5, 5, 70]
y_correction = [-155, -15, -15, -15]
for i, ax in enumerate(axes):
    x = ax.get_window_extent().x0 + x_correction[i]
    y = ax.get_window_extent().y1 + y_correction[i]
    x, y = fig.transFigure.inverted().transform([x, y])
    f = fig.text(x, y, labels[i], ha='left', va='bottom')

fig.set_size_inches(*figsize)
plt.savefig(paths.figures / "llzo.pdf", bbox_inches="tight", transparent=True)
plt.close()

