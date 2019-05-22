---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.1'
      jupytext_version: 1.1.2
  kernelspec:
    display_name: transmission
    language: python
    name: transmission3
---

<!-- #region markdown {} -->
# Demonstration of parameter estimation using Approximate Bayesian Computation





<!-- #endregion -->


```python
from __future__ import print_function, division

import itertools

import msprime as ms
import transmission as txmn
import numpy as np
import pandas as pd

%matplotlib inline
```


<!-- #region markdown {} -->
 First, simulate a data set. The `sim()` function in Transmission is the
 workhorse function that simulates a geneaology given a set of parameters.
 Generally, it is called by `ms_simulate()` and there is no reason to call
 it directly, but here we can use it for a proof of concept.
<!-- #endregion -->

```python
eta = 0
tau = 0.75
rho = 0.55


prior_seed = 3
random_seed = 3
host_theta = 1.5
npop = 10
nchrom = 24
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
    stats=("fst_mean", "fst_sd", "pi_h", "num_sites"),
    num_replicates=num_replicates,
    random_seed=random_seed
)
                          
simulated_target.dtype = np.dtype({"names": ("fst_mean", "fst_sd", "pi_h",
                                             "num_sites", "eta", "tau", "rho"),
                                   "formats": ['f8' for _ in range(7)]})
target_df = pd.DataFrame.from_records(simulated_target)
target_df
```

```python

```
