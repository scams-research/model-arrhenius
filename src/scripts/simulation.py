import sys
import itertools
import numpy as np
import scipp as sc
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from kinisi.analyze import ConductivityAnalyzer
from kinisi.arrhenius import Arrhenius, VogelFulcherTammann

from _fig_params import ALPHABET, CREDIBLE_INTERVALS, mid_points, PAGE_WIDTH, GOLDEN_RATIO
import _fig_params as fp
import paths

for length in [40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140]:
    D = []
    temp = np.array([300, 350, 400, 500, 600, 700])
    for i in temp:
        d = ConductivityAnalyzer.from_hdf5(paths.data / f"reduced/agcrse2/Dmsd{i}_{length}.h5")
        D.append(d)

    samples1 = np.load(paths.data / f"reduced/agcrse2/bayes_{length}_0.npz")['s_samples']
    samples2 = np.load(paths.data / f"reduced/agcrse2/bayes_{length}_0.npz")['t_samples']

    td0 = sc.DataArray(
    data=sc.array(dims=['temperature'],
                  values=[sc.mean(sc.to_unit(D[i].sigma, 'mS/cm') * (t * sc.Unit('K'))).value for i, t in enumerate(temp)],
                  variances=[sc.var(sc.to_unit(D[i].sigma, 'mS/cm') * (t * sc.Unit('K')), ddof=1).value for i, t in enumerate(temp)],
                  unit=D[0].sigma.unit * sc.Unit('K')),
    coords={
        'temperature': sc.Variable(dims=['temperature'],
                                   values=temp,
                                   unit='K')
    })

    s = Arrhenius(td0, bounds=((0 * sc.Unit('eV'), 1 * sc.Unit('eV')), (1e4 * td0.data.unit, 1e8 * td0.data.unit)))
    t = VogelFulcherTammann(td0, 
                        bounds=((0 * sc.Unit('eV'), 1 * sc.Unit('eV')), 
                                (1e4 * td0.data.unit, 1e8 * td0.data.unit), 
                                (0 * sc.Unit('K'), 300 * sc.Unit('K'))))

    size_adapt = 1
    figsize = (8.06, 8)
    label_offset = [0, 1 + 0.1 * size_adapt]
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(4, 8, figure=fig, wspace=1.2, hspace=1)

    axes = []
    lines = []
    titles = []

    for i in range(3):
        jj = 0
        for j in range(4):
            axes.append(fig.add_subplot(gs[i, jj:jj + 2]))
            jj += 2
    i = 0
    for ax in axes[::2]:
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
        ax.set_ylim(0, None)
        i += 1
    lines.append(line1)
    lines.append(line2)
    i = 0
    for ax in axes[1::2]:
        y, x = np.histogram(sc.to_unit(D[i].sigma, 'mS/cm').values, density=True, bins=fp.NBINS)
        hist = ax.stairs(y * 1e2, x, fill=True, color=fp.colors[1])
        if y.max() < 1.5e-2:
            ax.set_yticks([0, 1])
        ax.set_xlabel(r"$\hat{\sigma}_{\mathrm{Ag}^+}$ / mS cm$^{-1}$")
        ax.set_ylabel(r"$p(D^*)$ / $10^{-2}$ mS cm$^{-1}$")
        i += 1
    lines.append(hist)
    for i, ax in enumerate(axes):
        if i % 2 == 0:
            titles.append(r"$\mathbf{" + f"{ALPHABET[i]}" +
                            r"}$ - $T = {" + f"{temp[int(i/2)]}" +
                            r"}$ K")
        else:
            titles.append(r"$\mathbf{" + f"{ALPHABET[i]}" + r"}$ - $T = {" +
                            f"{temp[int(i/2)]}" + r"}$ K")
    axes.append(fig.add_subplot(gs[3, :3]))
    axes[-1].errorbar(1 / td0.coords['temperature'].values * 1000, td0.data.values, td0.data.variances ** 0.5, marker='.', ls='', c=fp.colors[1])
    x = np.linspace(295, 729, 1000)
    k = s.function(x[:, np.newaxis], *samples1)
    line1 = axes[-1].fill_between(1 / x * 1000,
                                *np.percentile(k, CREDIBLE_INTERVALS[0][:-1], axis=1),
                                alpha=CREDIBLE_INTERVALS[0][-1],
                                color=fp.colors[3],
                                lw=0)
    for j in CREDIBLE_INTERVALS[1:]:
        axes[-1].fill_between(1 / x * 1000, *np.percentile(k, j[:-1], axis=1), alpha=j[-1], color=fp.colors[3], lw=0)
    lines.append(line1)
    axes[-1].set_yscale('log', base=np.e)
    axes[-1].set_ylim(5e3, 5e5)
    axes[-1].set_yticks([1e4, 1e5])
    axes[-1].set_yticklabels([f"$10^{np.log10(val):.0f}$" for val in [1e4, 1e5]]) 
    axes[-1].set_xlabel(r"$1000T^{-1}$/K$^{-1}$")
    axes[-1].set_ylabel(r"$\hat{\sigma}_{\mathrm{Ag}^+}(T)T$ / mS cm$^{-1}$ K")
    titles.append(r"$\mathbf{" + f"{ALPHABET[12]}" + r"}$ - Arrhenius model")

    axes.append(fig.add_subplot(gs[3, 3:-2]))
    axes[-1].errorbar(1 / td0.coords['temperature'].values * 1000, td0.data.values, td0.data.variances ** 0.5, marker='.', ls='', c=fp.colors[1])    
    x = np.linspace(temp.min(), temp.max(), 1000)
    k = t.function(x[:, np.newaxis], *samples2)
    line2 = axes[-1].fill_between(1 / x * 1000,
                                *np.percentile(k, CREDIBLE_INTERVALS[0][:-1], axis=1),
                                alpha=CREDIBLE_INTERVALS[0][-1],
                                color=fp.colors[4],
                                lw=0)
    for j in CREDIBLE_INTERVALS[1:]:
        axes[-1].fill_between(1 / x * 1000, *np.percentile(k, j[:-1], axis=1), alpha=j[-1], color=fp.colors[4], lw=0)
    lines.append(line2)
    axes[-1].set_xlabel(r"$1000T^{-1}$/K$^{-1}$")
    axes[-1].set_yscale('log', base=np.e)
    axes[-1].set_ylim(5e3, 5e5)
    axes[-1].set_yticks([1e4, 1e5])
    axes[-1].set_yticks([])
    titles.append(r"$\mathbf{" + f"{ALPHABET[13]}" + r"}$ - VTF model")

    axes.append(fig.add_subplot(gs[3, -2:]))
    hist1 = axes[-1].hist(samples1[0], density=True, color=fp.colors[3], bins=fp.NBINS)
    hist2 = axes[-1].hist(samples2[0], density=True, color=fp.colors[4], bins=fp.NBINS)
    axes[-1].set_xlabel('$E_a$/eV')
    axes[-1].set_ylabel('$p(E_a)$/eV')
    titles.append(r"$\mathbf{" + f"{ALPHABET[14]}" + r"}$ - $E_a$ distribution")

    labels = [
        r'Observed $\langle\Delta\mathbf{r}_c(t)^2_{\mathrm{Ag}^+}\rangle$', 'Einstein Relation', r'$\hat{D}^*$', 'Arrhenius Relation',
        'VTF Equation'
    ]

    fig.align_ylabels(axes)

    def flip(items, ncol):
        return itertools.chain(*[items[i::ncol] for i in range(ncol)])


    legend = plt.figlegend(lines, labels, loc='upper center', bbox_to_anchor=(0.5, 0.05), ncol=5)

    x, _ = mid_points(axes[0])
    y = fig.axes[0].get_window_extent().y1 + 40
    x, y = fig.transFigure.inverted().transform([x, y])
    fig.text(x, y, r'$\langle\Delta\mathbf{r}_c(t)^2_{\mathrm{Ag}^+}\rangle$', ha='center', fontweight='bold')
    x, _ = mid_points(axes[1])
    y = fig.axes[1].get_window_extent().y1 + 40
    x, y = fig.transFigure.inverted().transform([x, y])
    fig.text(x, y, r'$\hat{D}^*_{\mathrm{Ag}^+}$', ha='center', fontweight='bold')
    x, _ = mid_points(axes[2])
    y = fig.axes[2].get_window_extent().y1 + 40
    x, y = fig.transFigure.inverted().transform([x, y])
    fig.text(x, y, r'$\langle\Delta\mathbf{r}_c(t)^2_{\mathrm{Ag}^+}\rangle$', ha='center', fontweight='bold')
    x, _ = mid_points(axes[3])
    y = fig.axes[3].get_window_extent().y1 + 40
    x, y = fig.transFigure.inverted().transform([x, y])
    fig.text(x, y, r'$\hat{D}^*_{\mathrm{Ag}^+}$', ha='center', fontweight='bold')

    x_correction = [28, 28, 28, 28] * 3 + [30, 15, 30]
    for i, ax in enumerate(axes):
        x = ax.get_window_extent().x0 - x_correction[i]
        y = ax.get_window_extent().y1 + 10
        x, y = fig.transFigure.inverted().transform([x, y])
        f = fig.text(x, y, titles[i], ha='left')

    plt.savefig(paths.figures / f"agcrse2_{length}.pdf", transparent=True)
    plt.close()