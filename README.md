# Bayesian Methods for the Investigation of Temperature-Dependence in Conductivity

<p align="justify">
Temperature-dependent transport data, including diffusion coefficients and ionic conductivities, are routinely analysed by fitting empirical models such as the Arrhenius equation. 
These fitted models yield parameters such as the activation energy, and can be used to extrapolate to temperatures outside the measured range.
Researchers frequently face challenges in this analysis: quantifying the uncertainty of fitted parameters, assessing whether the data quality is sufficient to support a particular empirical model, and using these models to predict behaviour at extrapolated temperatures. 
Bayesian methods offer a coherent framework that addresses all of these challenges. 
This tutorial introduces the use of Bayesian methods for analysing temperature-dependent transport data, covering parameter estimation, model selection, and extrapolation with uncertainty propagation, with illustrative examples from molecular dynamics simulations of superionic materials.
</p>

---

<p align="center">
<a href="https://github.com/arm61/linearization-issues/actions/workflows/build.yml">
<img src="https://github.com/scams-research/model-arrhenius/actions/workflows/build.yml/badge.svg" alt="Article status"/>
</a>
<a href="https://github.com/scams-research/model-arrhenius/raw/main-pdf/arxiv.tar.gz">
<img src="https://img.shields.io/badge/article-tarball-blue.svg?style=flat" alt="Article tarball"/>
</a>
<a href="https://github.com/scams-research/model-arrhenius/raw/main-pdf/ms.pdf">
<img src="https://img.shields.io/badge/article-pdf-blue.svg?style=flat" alt="Read the article"/>
<!-- </a>
<a href="https://doi.org/10.5281/zenodo.7949905">
<img src="https://zenodo.org/badge/DOI/10.5281/zenodo.7949905.svg"/>
</a> -->
<!-- <a href="https://doi.org/10.26434/chemrxiv-2023-44b29">
<img src="https://img.shields.io/badge/ChemRxiv-10.26434%2Fchemrxiv--2023--44b29-orange.svg"/>
</a> -->
<br><br>
<a href="https://orcid.org/0000-0003-3381-5911">Andrew R. McCluskey</a>&ast;
<a href="https://orcid.org/0000-0001-9722-5676">Samuel W. Coles</a>
<a href="https://orcid.org/0000-0002-3056-8233">Benjamin J. Morgan</a>&dagger;<br>
&ast;<a href="mailto:andrew.mccluskey@bristol.ac.uk">andrew.mccluskey@bristol.ac.uk</a>
&dagger;<a href="mailto:b.j.morgan@bath.ac.uk">b.j.morgan@bath.ac.uk</a>
</p>

---

This is the electronic supplementary information (ESI) associated with the publication "Bayesian Methods for the Investigation of Temperature-Dependence in Conductivity". 
This ESI uses [`showyourwork`](https://show-your.work) to provide a completely reproducible and automated analysis, plotting, and paper generation workflow. 
To run the workflow and generate the paper locally using the cached data run the following: 
```
git clone git@github.com:scams-research/model-arrhenius.git
cd model-arrhenius 
pip install showyourwork
showyourwork build 
```
Full details of the workflow can be determined from the [`Snakefile`](https://github.com/scams-research/model-arrhenius/blob/main/Snakefile) and the [`showyourwork.yml`](https://github.com/scams-research/model-arrhenius/blob/main/showyourwork.yml).