import numpy as np
from scipy.stats import normaltest
from kinisi.analyze import ConductivityAnalyzer
import paths

SIGMA = 0.05

samples = np.load(paths.data / "reduced/llzo/arrhenius_samples.npz")['samples']

write_out = open(paths.output / "activation_energy_llzo.txt", "w")
ea = np.array(samples[0])
write_out.write(r"$\SI{" + f"{ea.mean():.3f}" + r"\pm" + f"{ea.std() * 1.96:.3f}" + r"}{\electronvolt}$")
write_out.close()

write_out = open(paths.output / "preexp_llzo.txt", "w")
a = np.array(samples[1] * 1e-6)
ci = np.percentile(a, [2.5, 50, 97.5])
write_out.write(r"$\left(" + f"{ci[1]:.3f}" + r"^{+" + f"{ci[2] - ci[1]:.3f}" + r"}_{-" + f"{ci[1] - ci[0]:.3f}" + r"}\right)\times10^{6}\,\si{\centi\meter^2\second^{-1}}$")
write_out.close()

samples = np.load(paths.data / "reduced/llzo/arrhenius_samples.npz")['extrapolated_temperature']

write_out = open(paths.output / "d_300.txt", "w")
a = np.array(samples / 300)
ci = np.percentile(a, [2.5, 50, 97.5])
write_out.write(r"$\big(" + f"{ci[1]:.1f}" + r"^{+" + f"{ci[2] - ci[1]:.1f}" + r"}_{-" + f"{ci[1] - ci[0]:.1f}" + r"}\big)\,\si{\milli\siemens\centi\meter^{-1}}$")
write_out.close()

lengths = [40, 140]
n = np.zeros((10, len(lengths)))
for i, l in enumerate(lengths):
    for j in range(n.shape[0]):
        n[j, i] = np.load(paths.data / f"reduced/agcrse2/bayes_{l}_{j}.npz")['bf'][0] / 2
        

write_out = open(paths.output / "bayes_40.txt", "w")
write_out.write(r"$\num{" + f"{np.mean(n[:, 0]):.1f}" + r"\pm" + f"{np.std(n[:, 0]) * 1.96:.1f}" + r"}$")
write_out.close()

write_out = open(paths.output / "bayes_140.txt", "w")
write_out.write(r"$\num{" + f"{np.mean(n[:, -1]):.1f}" + r"\pm" + f"{np.std(n[:, -1]) * 1.96:.1f}" + r"}$")
write_out.close()

temps = [300, 350, 400, 500, 600, 700]
for i in temps:
	d = ConductivityAnalyzer.from_hdf5(paths.data / f"reduced/agcrse2/Dmsd{i}_{lengths[0]}.h5")
	write_out = open(paths.output / f"dt{i}.txt", "w")
	write_out.write(r"\num{" + f"{d.dt.min().value * 5:.2f}" + r"}")
	write_out.close()
