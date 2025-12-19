import numpy as np
import scipp as sc
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from uncertainties import ufloat
from kinisi.analyze import ConductivityAnalyzer
from kinisi.arrhenius import Arrhenius, VogelFulcherTammann

from _fig_params import CREDIBLE_INTERVALS, mid_points, COLUMN_WIDTH, GOLDEN_RATIO
import _fig_params as fp
import paths

lengths = [40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140]
temp = np.array([300, 350, 400, 500, 600, 700])
D0 = {'mean': [], 'var': []}
D1 = {'mean': [], 'var': []}
for i in temp:
    d0 = ConductivityAnalyzer.from_hdf5(paths.data / f"reduced/agcrse2/Dmsd{i}_{lengths[0]}.h5")
    d1 = ConductivityAnalyzer.from_hdf5(paths.data / f"reduced/agcrse2/Dmsd{i}_{lengths[-1]}.h5")
    sigmaT0 = sc.to_unit(d0.sigma, 'mS/cm') * (i * sc.Unit('K'))
    sigmaT1 = sc.to_unit(d1.sigma, 'mS/cm') * (i * sc.Unit('K'))
    D0['mean'].append(sc.mean(sigmaT0).value)
    D0['var'].append(sc.var(sigmaT0, ddof=1).value)
    D1['mean'].append(sc.mean(sigmaT1).value)
    D1['var'].append(sc.var(sigmaT1, ddof=1).value)
    unit = sigmaT0.unit


n = np.zeros((10, len(lengths)))
for i, l in enumerate(lengths):
    for j in range(n.shape[0]):
        n[j, i] = np.load(paths.data / f"reduced/agcrse2/bayes_{l}_{j}.npz")['bf'][0] / 2

td0 = sc.DataArray(
    data=sc.array(dims=['temperature'],
                  values=D0['mean'],
                  variances=D0['var'],
                  unit=unit),
    coords={
        'temperature': sc.Variable(dims=['temperature'],
                                   values=temp,
                                   unit='K')
    })

td1 = sc.DataArray(
    data=sc.array(dims=['temperature'],
                  values=D1['mean'],
                  variances=D1['var'],
                  unit=unit),
    coords={
        'temperature': sc.Variable(dims=['temperature'],
                                   values=temp,
                                   unit='K')
    })

figsize = (COLUMN_WIDTH, COLUMN_WIDTH / GOLDEN_RATIO * 3) 
fig = plt.figure(figsize=figsize)
gs = gridspec.GridSpec(3, 2, figure=fig, wspace=0.1, hspace=0.6)
axes = []
titles = []

s = Arrhenius(td0, bounds=((0 * sc.Unit('eV'), 1 * sc.Unit('eV')), (1e4 * td0.data.unit, 1e8 * td0.data.unit)))
t = VogelFulcherTammann(td0, 
                        bounds=((0 * sc.Unit('eV'), 1 * sc.Unit('eV')), 
                                (1e4 * td0.data.unit, 1e8 * td0.data.unit), 
                                (0 * sc.Unit('K'), 300 * sc.Unit('K'))))
samples1 = np.zeros((n.shape[0], 2))
samples2 = np.zeros((n.shape[0], 3))
for i in range(n.shape[0]):
    samples1[i] = np.load(paths.data / f"reduced/agcrse2/bayes_{lengths[0]}_{i}.npz")['s_samples'].mean(1)
    samples2[i] = np.load(paths.data / f"reduced/agcrse2/bayes_{lengths[0]}_{i}.npz")['t_samples'].mean(1)

smooth = np.linspace(295, 729, 1000)
axes.append(fig.add_subplot(gs[0, 0]))
axes[-1].errorbar(1 / td0.coords['temperature'].values * 1000,
                  td0.data.values,
                  td0.data.variances ** 0.5,
                  marker=".",
                  ls="",
                  c='k',
                  zorder=10)
ss = np.load(paths.data / f"reduced/agcrse2/bayes_{lengths[0]}_0.npz")['s_samples']
k = s.function(smooth[:, np.newaxis], *ss)
axes[-1].fill_between(1 / smooth * 1000,
                      *np.percentile(k, CREDIBLE_INTERVALS[0][:-1], axis=1),
                      alpha=CREDIBLE_INTERVALS[0][-1],
                      color=fp.colors[0],
                      lw=0,
                      label=r'$p[\theta_{\alpha} \, | \, (\sigma^{*}_{\mathrm{Ag^+}}\,|\,\mathbf{x})]$')
for j in CREDIBLE_INTERVALS[1:]:
    axes[-1].fill_between(1 / smooth * 1000, *np.percentile(k, j[:-1], axis=1), alpha=j[-1], color=fp.colors[0], lw=0)
axes[-1].set_yscale("log", base=np.e)
axes[-1].set_xlabel("$T^{-1}$ / 10$^{-3}$ K$^{-1}$")
axes[-1].set_ylabel(r"$\hat{\sigma}_{\mathrm{Ag}^+}$ / mS cm$^{-1}$ K")
axes[-1].set_ylim(5e3, 5e5)
axes[-1].set_yticks([1e4, 1e5])
axes[-1].set_yticklabels([f"$10^{np.log10(val):.0f}$" for val in [1e4, 1e5]])
# axes[-1].xaxis.set_major_locator(plt.MaxNLocator(3)) 
titles.append(r"a - $t_{\mathrm{diff}} = $40 ps")

axes.append(fig.add_subplot(gs[0, 1]))
ss = np.load(paths.data / f"reduced/agcrse2/bayes_{lengths[0]}_0.npz")['t_samples']
axes[-1].errorbar(1 / td0.coords['temperature'].values * 1000,
                  td0.data.values,
                  td0.data.variances ** 0.5, marker=".", ls="", c='k', zorder=10)
k = t.function(smooth[:, np.newaxis], *ss)
axes[-1].fill_between(1 / smooth * 1000, *np.percentile(k, CREDIBLE_INTERVALS[0][:-1], axis=1), alpha=CREDIBLE_INTERVALS[0][-1], color=fp.colors[1], lw=0,
                      label=r'$p[\theta_{\beta} \, | \, (\sigma_{\mathrm{Ag^+}}\,|\,\mathbf{x})]$')
for j in CREDIBLE_INTERVALS[1:]:
    axes[-1].fill_between(1 / smooth * 1000, *np.percentile(k, j[:-1], axis=1), alpha=j[-1], color=fp.colors[1], lw=0)
axes[-1].set_xlabel("$T^{-1}$ / 10$^{-3}$ K$^{-1}$")
axes[-1].set_yscale("log", base=np.e)
axes[-1].set_ylim(5e3, 5e5)
axes[-1].set_yticks([])
# axes[-1].xaxis.set_major_locator(plt.MaxNLocator(3))
titles.append(r"b - $t_{\mathrm{diff}} = $40 ps")

s = Arrhenius(td1, bounds=((0 * sc.Unit('eV'), 1 * sc.Unit('eV')), (1e4 * td1.data.unit, 1e8 * td1.data.unit)))
t = VogelFulcherTammann(td1, 
                        bounds=((0 * sc.Unit('eV'), 1 * sc.Unit('eV')), 
                                (1e4 * td0.data.unit, 1e8 * td1.data.unit), 
                                (0 * sc.Unit('K'), 300 * sc.Unit('K'))))

samples1 = np.zeros((n.shape[0], 2))
samples2 = np.zeros((n.shape[0], 3))
for i in range(n.shape[0]):
    samples1[i] = np.load(paths.data / f"reduced/agcrse2/bayes_{lengths[-1]}_{i}.npz")['s_samples'].mean(1)
    samples2[i] = np.load(paths.data / f"reduced/agcrse2/bayes_{lengths[-1]}_{i}.npz")['t_samples'].mean(1)

axes.append(fig.add_subplot(gs[1, 0]))
axes[-1].errorbar(1 / td1.coords['temperature'].values * 1000, td1.data.values, td1.data.variances ** 0.5, marker=".", ls="", c='k', zorder=10,
                  label=r'Estimated $\sigma_{\mathrm{Ag}^+}$')
ss = np.load(paths.data / f"reduced/agcrse2/bayes_{lengths[-1]}_0.npz")['s_samples']
k = s.function(smooth[:, np.newaxis], *ss)
axes[-1].fill_between(1 / smooth * 1000,
                      *np.percentile(k, CREDIBLE_INTERVALS[0][:-1], axis=1),
                      alpha=CREDIBLE_INTERVALS[0][-1],
                      color=fp.colors[0],
                      lw=0)
for j in CREDIBLE_INTERVALS[1:]:
    axes[-1].fill_between(1 / smooth * 1000, *np.percentile(k, j[:-1], axis=1), alpha=j[-1], color=fp.colors[0], lw=0)
axes[-1].set_xlabel("$T^{-1}$ / 10$^{-3}$ K$^{-1}$")
axes[-1].set_ylabel(r"$\hat{\sigma}_{\mathrm{Ag}^+}$ / mS cm$^{-1}$ K")
axes[-1].set_yscale("log", base=np.e)
axes[-1].set_ylim(5e3, 5e5)
axes[-1].set_yticks([1e4, 1e5])
axes[-1].set_yticklabels([f"$10^{np.log10(val):.0f}$" for val in [1e4, 1e5]])
# axes[-1].xaxis.set_major_locator(plt.MaxNLocator(3))
titles.append(r"c - $t_{\mathrm{diff}} = $140 ps")

axes.append(fig.add_subplot(gs[1, 1]))
axes[-1].errorbar(1 / td1.coords['temperature'].values * 1000,
                  td1.data.values,
                  td1.data.variances ** 0.5, marker=".", ls="", c='k', zorder=10)
ss = np.load(paths.data / f"reduced/agcrse2/bayes_{lengths[-1]}_0.npz")['t_samples']
k = t.function(smooth[:, np.newaxis], *ss)
for j in CREDIBLE_INTERVALS:
    axes[-1].fill_between(1 / smooth * 1000, *np.percentile(k, j[:-1], axis=1), alpha=j[-1], color=fp.colors[1], lw=0)
axes[-1].set_xlabel("$T^{-1}$ / 10$^{-3}$ K$^{-1}$")
axes[-1].set_yscale("log", base=np.e)
axes[-1].set_ylim(5e3, 5e5)
axes[-1].set_yticks([])
# axes[-1].set_ylabel("$\hat{D}^*_{\mathrm{Ag}^+}$/cm$^2$s$^{-1}$")
# axes[-1].xaxis.set_major_locator(plt.MaxNLocator(3))
titles.append(r"d - $t_{\mathrm{diff}} = $140 ps")

axes.append(fig.add_subplot(gs[2, :]))
axes[-1].axhline(5, c=fp.colors[2], ls='--', label=r'$\ln(B_{{\beta\alpha}})=5$')
axes[-1].errorbar(lengths,
                  n.mean(0),
                  n.std(0, ddof=1) * 1.96,
                  marker=".",
                  ls="",
                  zorder=10,
                  c=fp.colors[3])
axes[-1].set_xlabel("Diffusive simulation length / ps")
axes[-1].set_ylabel(r"$\ln(B_{\beta\alpha})$")
# axes[-1].set_ylim(-10, None)
axes[-1].set_xticks([40, 90, 140])
axes[-1].yaxis.set_major_locator(plt.MaxNLocator(3))
titles.append(r"e - Bayes' factor vs. diffusive simulation time")

x, _ = mid_points(axes[0])
y = fig.axes[0].get_window_extent().y1 + 60
x, y = fig.transFigure.inverted().transform([x, y])
fig.text(x, y, 'Arrhenius', ha='center', fontweight='bold')
x, _ = mid_points(axes[1])
y = fig.axes[1].get_window_extent().y1 + 60
x, y = fig.transFigure.inverted().transform([x, y])
fig.text(x, y, 'VTF', ha='center', fontweight='bold')

fig.align_ylabels(axes)

x_correction = [39, 10, 39, 10, 39]
for i, ax in enumerate(axes):
    x = ax.get_window_extent().x0 - x_correction[i]
    y = ax.get_window_extent().y1 + 20
    x, y = fig.transFigure.inverted().transform([x, y])
    fig.text(x, y, titles[i], ha='left')

# plt.setp(axes[1].get_yticklabels(), visible=False)
# plt.setp(axes[3].get_yticklabels(), visible=False)
# axes[1].minorticks_off()
# axes[3].minorticks_off()

handles = []
labels = []
for ax in axes[:4]:
    h, l = ax.get_legend_handles_labels()
    handles += h
    labels += l

plt.figlegend(handles[:2], labels[:2],
              loc='upper center',
              bbox_to_anchor=(0.5, 0.975),
              ncol=2)
handles = []
labels = []
for ax in axes[4:]:
    h, l = ax.get_legend_handles_labels()
    handles += h
    labels += l


axes[-1].legend(handles, labels,
              loc='upper left',
              ncol=1)

fig.set_size_inches(*figsize)
plt.savefig(paths.figures / "agcrse2.pdf", bbox_inches="tight", transparent=True)
plt.close()

# l1 = np.load(paths.data / f"reduced/bayes_{lengths[0]}_0.npz")['samples'][0]
# v = ufloat(l1.mean(), l1.std() * 1.96)
# write_out = open(paths.output / "agcrse2_standard_ea.txt", "w")
# write_out.write(r"$\SI{" + f"{v:.3fL}" + r"}{\electronvolt}$")
# write_out.close()
# l1 = np.load(paths.data / f"reduced/bayes_{lengths[-1]}_0.npz")['ssamples'][0]
# v = ufloat(l1.mean(), l1.std() * 1.96)
# write_out = open(paths.output / "agcrse2_super_ea.txt", "w")
# write_out.write(r"$\SI{" + f"{v:.3fL}" + r"}{\electronvolt}$")
# write_out.close()