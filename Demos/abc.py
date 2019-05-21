# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.2'
#       jupytext_version: 1.1.2
#   kernelspec:
#     display_name: transmission
#     language: python
#     name: transmission3
# ---

# %% markdown [markdown]
# # Demonstration of parameter estimation using Approximate Bayesian Computation
#
#
#
#
#


# %%
from __future__ import print_function, division

import itertools

import msprime as ms
import transmission as txmn
import numpy as np
import pandas as pd

# %matplotlib inline


# %% markdown [markdown]
#  First, simulate a data set. The `sim()` function in Transmission is the
#  workhorse function that simulates a geneaology given a set of parameters.
#  Generally, it is called by `ms_simulate()` and there is no reason to call
#  it directly, but here we can use it for a proof of concept.

# %%
eta = 0
tau = 0.75
rho = 0.55


prior_seed = 3
random_seed = 3
host_theta = 1.5
npop = 10
nchrom = 10
host_Nm = 0.56  # host F_(ST-mt) = 0.64
num_replicates = 24

population_config = [ms.PopulationConfiguration(nchrom)
                     for _ in range(npop)]
populations = np.repeat(range(npop), nchrom)

simulated_target = txmn.sim(
    (eta, tau, rho),
    host_theta=host_theta,
    host_Nm=host_Nm,
    population_config=population_config,
    populations=populations,
    stats=("fst_mean", "fst_sd", "pi_h"),
    num_replicates=num_replicates,
    random_seed=random_seed
)
                          

# %%
simulated_target

# %%
