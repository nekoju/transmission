#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Usage ./multi.py num_cores num_repeats num_simulations

from textwrap import dedent
import timeit
import sys


if __name__ == "__main__":

    out_time = timeit.timeit(
        dedent(
            """
            txmn.generate_priors(
                nchrom=nchrom,
                num_populations=npop,
                host_theta=host_theta,
                host_Nm=host_Nm,
                num_simulations=num_simulations,
                num_replicates=num_replicates,
                num_cores=num_cores,
                prior_seed=3,
                progress_bar=True,
                random_seed=random_seed,
            )
            """
        ),
        setup=dedent(
            """
            import transmission as txmn
            import msprime as ms
            

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
            num_simulations = {num_simulations}
            num_cores = {num_cores}

            # Create populations using msprime API
            population_config = [
                ms.PopulationConfiguration(nchrom) for _ in range(npop)
            ]
            # Gives population identity (0 -- npop - 1) to each sampled
            # chromosome, 0, 1, 2, ...
            """.format(
                **{
                    "num_cores": int(sys.argv[1]),
                    "num_simulations": int(sys.argv[3]),
                }
            )
        ),
        number=int(sys.argv[2]),
    )
    print(out_time)
