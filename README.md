# PTAsimulations_from_SMBHBpop

In this repository we stored all the main scripts used for the project described in the paper:

["Testing strengths, limitations and biases of current Pulsar Timing Arrays detection analyses on realistic data"](https://arxiv.org/abs/2309.13117).

The correspondent datasets can be found at  [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10276364.svg)](https://doi.org/10.5281/zenodo.10276364)

If you make use of any of these scripts/data, please cite:
```
@article{Valtolina:2023axp,
    author = "Valtolina, Serena and Shaifullah, Golam and Samajdar, Anuradha and Sesana, Alberto",
    title = "{Testing strengths, limitations and biases of current Pulsar Timing Arrays detection analyses on realistic data}",
    eprint = "2309.13117",
    archivePrefix = "arXiv",
    primaryClass = "astro-ph.HE",
    month = "9",
    year = "2023"
}
```

## Repository content details

Two tutorials are included in this repo (jupyter notebooks):

- **mcmc_and_reweighting.ipynb** : how to run an inference run with mcmc looking for a common uncorrelated red noise; then we derive the posteriors for GWB amplitude and slope through the re-weighting method described in [Hourihane et al. 2023](https://arxiv.org/abs/2212.06276)
- **HD_frequentist_recovery.ipynb** : how to recover the HD correlation between timing residuals using the Optimal Statistic (see [Anholm et al. 2008](https://arxiv.org/abs/0809.0701)) tools (frequentist approach); the average over angular separation bins takes into account for the covariance between pulsar pairs ([Allen & Romano 2023](https://arxiv.org/abs/2208.07230))

Here's the list of scripts contained in this repository (and the reference to the result produced for our analysis):

- **combined_posterior.py** : script to compute the combined posterior, Bayesian combination weighted over prior distributions *(Fig. 6)*
- **Lrecovery.py** : likelihood evaluations over prior distribution for one free parameter (the results are then used to build the *fast* PPplot)
- **myHD_universe_avg.py** : HD correlation recovery averaged over many GWB realisations *(Fig. 4 and 8)*
- **optimal_statistic_covariances.py** : OS functions that include the weighting over covariance matrices ([source](https://gitlab.com/IPTA/ng_15yr_gwb_analysis_code/-/blob/main/prelim_code/optimal_statistic_covariances.py))
- **PPplot_fast_script.py** : PP plot from likelihood evaluations (inference run over one parameter only) *(dashed lines in Fig. 2 and 5)*
- **PPplot_script.py** : PP plot from mcmc posteriors *(solid lines in Fig. 2 and 5)*

*NOTE:* all the paths defined in the scripts of this repository refer to the datasets published in Zenodo. By saving *SMBHB_datasets.zip* in the same directory of the scripts, it is possible (after the 'unzip') to run all of them without changes.

Here's the list of dictionaries:

- **ml_params.json** : maximum likelihood noise parameters dictionary; the values comes from EPTA DR2new from inference runs where a common red noise was also included
- **white_noise_dict.json** : white noise dictionary (EFAC and EQUAD for each pulsar)
- **complete_noise_dict.json** : combination of *ml_params.json* and *white_noise_dict.json*

## Before running

Those scripts use the EPTA branch of ENTERPRISE (Enhanced Numerical Toolbox Enabling a Robust PulsaR Inference SuitE), a pulsar timing analysis code aimed at noise analysis, gravitational-wave searches, and timing model analysis.
```
git clone https://gitlab.in2p3.fr/epta/enterprise_extensions/
git clone https://gitlab.in2p3.fr/epta/enterprise/
```


