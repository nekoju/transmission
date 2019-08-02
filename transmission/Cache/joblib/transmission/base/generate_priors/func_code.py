# first line: 49
def generate_priors(
    nchrom,
    num_populations,
    host_theta,
    host_Nm,
    num_simulations,
    stats=("fst_mean", "fst_sd", "pi_h"),
    prior_params={"eta": (0, 0.1), "tau": (1, 1), "rho": (1, 1)},
    nsamp_populations=None,
    num_replicates=1,
    num_cores="auto",
    prior_seed=None,
    average_reps=True,
    progress_bar=False,
    h_opts={},
    **kwargs
):
    """
    Generate random sample summary using msprime for the specified prior
    distributions of tau (vertical transmission rate) and rho (sex ratio).

    From a supplied dict of prior (hyper)parameters, generate_priors generates a
    num_simulations-tall array of parameters with which to estimate the
    posterior distributions of tau and rho using abc methods. It then draws
    num_replicates coalescent samples for each parameter combination and
    outputs a summary table consisting of mean fst, fst standard deviation, pi,
    and the original simulated parameter values for each metasimulation.

    Args:
        nchrom (int): The number of chromosomes to sample from each population
        num_populations (int): The number of populations to include in the
            migration matrix.
        host_theta (float): The host's haploid theta (2 * Ne * mu * L)
            estimated from mitochondrial or nuclear data.
        host_Nm (float): The symmetric migration parameter (2 * Ne * m),
            estimated from host mitochondrial or nuclear fst.
        num_simulations (int): The number of tau-rho pairs of parameters
            to draw from their priors, and hence the number of metasimulations
            to include in the output.
        stats (tuple): The statistics to be returned by the simulation.
            May take values "fst_mean" for mean locus Fst, "fst_sd" for locus
            standard deviation, or "pi_h", "pi_nei", or "pi_tajima". See
            Sample.pi() documentation for details.
        prior_params (dict): A dict containing tuples specifiying the prior
            distribution parameters for eta, tau, and rho. That is, the
            mutation rate multiplier, vertical transmission frequency, and
            sex ratio. Optionally one may provide a scalar value for eta to
            fix the mutation rate multiplier.
        num_replicates (int): Number of msprime replicates to run for each
            simulated parameter pair and the number of simulations in each
            metasimulation.
        average_reps (bool): Whether to average replicates. False will return
            raw simulations.
        progress_bar (bool): Display a progress bar for prior simulations.
            Redundant if num_cores is set to None.
        num_cores (int or str): The number of cores to use for computation.
            "auto" will automatically detect using multiprocessing.cpu_count().
            None will use a single thread not routed through
            multiprocessing.Pool, primarily for debugging and profiling.
            An int value will specify the number of cores to use.
        prior_seed (int): The seed used to draw samples for tau and rho.
            Setting a seed will allow repeatability of results.
        h_opts (dict): Extra options for Sample.h() in the form of a kwargs
            dictionary.
        **kwargs (): Extra arguments for ms.simulate(), Sample.pi(), and
            Sample.fst().
    """

    # Number of output statistics plus parameters eta, tau, and rho.

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
        _sim,
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
        chunksize = 10 
        pool = get_context("spawn").Pool(processes=num_cores)
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
        pool.close()
        pool.join()
    else:
        out = np.apply_along_axis(simpartial, 1, params)
    out = np.core.records.fromarrays(out.T, dtype=structure)
    return out
