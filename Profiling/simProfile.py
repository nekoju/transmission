from __future__ import print_function, division

from joblib import Memory
import msprime as ms
import transmission as txmn
import numpy as np
import pandas as pd

if __name__ == "__main__":
    memory = Memory("./Cache", verbose=0)
    eta = 0.15  # Exponent of 10 representing the multiplicative difference
    # between the host's mutation rate and the symbiont's.
    tau = 0.75  # The fraction of new infections that result from vertical
    # transmission.
    rho = 0.55  # The fraction of the population that is female.

    prior_seed = 3
    random_seed = 3
    host_theta = 1  # Estimated from the host mitochondria.
    npop = 10  # Number of populations
    nchrom = 10  # Number of chromosomes sampled from each population.
    host_Nm = 2
    num_replicates = 25

    # Create populations using msprime API
    population_config = [
        ms.PopulationConfiguration(nchrom) for _ in range(npop)
    ]
    # Gives population identity (0 -- npop - 1) to each sampled
    # chromosome, 0, 1, 2, ...
    populations = np.repeat(range(npop), nchrom)

    # The following takes a minute or so. We are generating a target using the
    # above parameters.
    simulated_target = txmn.sim(
        (eta, tau, rho),
        host_theta=host_theta,
        host_Nm=host_Nm,
        population_config=population_config,
        populations=populations,
        stats=("fst_mean", "fst_sd", "pi_h"),
        num_replicates=num_replicates,
        random_seed=random_seed,
    )
