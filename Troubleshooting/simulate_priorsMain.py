import functools
from math import ceil
from multiprocessing import get_context, cpu_count

import msprime as ms
import numpy as np
import pandas as pd
from tqdm.autonotebook import tqdm
import txmn as txmn


if __name__ == "__main__":

    eta = 0.15  # Exponent of 10 representing the multiplicative difference
                # between the host's mutation rate and the symbiont's.
    tau = 0.75  # The fraction of new infections that result from vertical
                # txmn.
    rho = 0.55  # The fraction of the population that is female.

    prior_seed = 3
    random_seed = 3
    host_theta = 1        # Estimated from the host mitochondria.
    num_populations = 10  # Number of populations
    nchrom = 10           # Number of chromosomes sampled from each population.
    host_Nm = 2
    num_replicates = 25
    npop = 10
    nsamp_populations = None
    num_simulations = 100
    prior_params = {"eta": (0, 0.1), "tau": (1, 1), "rho": (10, 10)}
    stats = ("fst_mean", "fst_sd", "pi_h")
    h_opts = {}
    kwargs = {}
    num_cores = "auto"
    progress_bar = True

    populations = np.repeat(np.arange(num_populations), nchrom)
    population_config = tuple(
        ms.PopulationConfiguration(nchrom) for _ in np.arange(num_populations)
    )
    nsamp_populations = (
        num_populations if not nsamp_populations else nsamp_populations
    )
    if prior_seed:
        np.random.seed(prior_seed)
    if isinstance(prior_params["eta"], float):
        eta = np.full((num_simulations,), prior_params["eta"])
    elif isinstance(prior_params["eta"], tuple):
        eta = np.random.normal(
            prior_params["eta"][0], prior_params["eta"][1], num_simulations
        )
    else:
        raise Exception("eta must be tuple or float")
    tau = np.random.beta(
        prior_params["tau"][0], prior_params["tau"][1], size=num_simulations
    )
    rho = np.random.beta(
        prior_params["rho"][0], prior_params["rho"][1], size=num_simulations
    )
    params = np.array([eta, tau, rho]).T
    simpartial = functools.partial(
        txmn.sim,
        host_theta=host_theta,
        host_Nm=host_Nm,
        population_config=population_config,
        num_replicates=num_replicates,
        populations=populations,
        average_reps=True,
        stats=stats,
        h_opts=h_opts,
        **kwargs
    )
    structure = {
        "names": stats + tuple(prior_params.keys()),
        "formats": tuple(np.repeat("f8", len(stats) + len(prior_params))),
    }
    if num_cores:
        # for multiprocessing, None automatically detects number of cores.
        # Switching here allows autodetection.
        if num_cores == "auto":
            num_cores = cpu_count()
        chunksize = ceil(num_simulations / num_cores)
        with get_context("spawn").Pool(processes=num_cores) as pool:
            if progress_bar:
                out = np.array(
                    list(
                        tqdm(
                            pool.imap_unordered(
                                simpartial, params, chunksize=chunksize
                            ),
                            total=len(params),
                        )
                    )
                )
            else:
                out = np.array(
                    list(
                        pool.imap_unordered(
                            simpartial, params, chunksize=chunksize
                        )
                    )
                )
    else:
        out = np.apply_along_axis(simpartial, 1, params)
    out = np.core.records.fromarrays(out.T, dtype=structure)

    priors = pd.DataFrame.from_records(
        txmn.generate_priors(
            nchrom=nchrom,
            num_populations=npop,
            host_theta=host_theta,
            host_Nm=host_Nm,
            num_simulations=40,
            num_replicates=num_replicates,
            prior_seed=3,
            progress_bar=True,
            random_seed=random_seed
        )
    )
    print(priors)
