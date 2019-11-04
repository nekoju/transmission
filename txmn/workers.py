import msprime as ms
import numpy as np

from txmn.classes import MetaSample


def _sim(
    params,
    host_theta,
    host_Nm,
    population_config,
    populations,
    stats,
    num_replicates,
    keep_populations=None,
    migration=None,
    average_reps=True,
    theta_source="mt",
    h_opts={},
    **kwargs
):
    """
    Runs actual simulation with ms. Intended as helper for generate_priors().

    At top level for picklability (for multiprocessing).

    Args:
        params (tuple): eta, tau, rho for simulation.
        host_theta (float): Estimate of host theta (2*Ne*mu*L) for host.
        host_Nm (float): Estimate of host migration parameter (Ne*m).
        population_config (list): List of population configurations for
            msprime.
        populations (np.ndarray): A nchrom np.ndarray indicating to which
            population each chromosome belongs.
        stats (tuple): The requested statistics to be calculated. May contain
            any of "fst_mean", "fst_sd", "pi_h", "pi_nei", "pi_tajima",
            "theta_w", or "num_sites".
        migration (np.ndarray): A np.ndarray with proportional migration
            rate among demes. These values will be multiplied by Ne*m for the
            actual migration rates.
        average_reps (bool): Whether to return averaged replicates. False will
            return raw summaries.
        num_replicates (int): Number of msprime replicates to run for each
            simulated parameter pair and the number of simulations in each
            metasimulation.
        h_opts (dict): Extra arguments for Sample.h(), formatted as kwargs
            dictionary.
        **kwargs (): Extra arguments for msprime.simulate().
    """

    eta, tau, rho = params
    a = rho if theta_source == "mt" else 1
    A = tau ** 2 * (3 - 2 * tau) * (1 - rho)
    B = 2 * rho * (1 - rho) * (A + rho)
    symbiont_Nm = np.true_divide(host_Nm, 2 * B)
    symbiont_theta = np.true_divide(10 ** eta * host_theta, 2 * a * B)
    num_populations = len(population_config)
    migration_null = np.full(
        (num_populations, num_populations),
        np.true_divide(
            # num_populations - 1 is multiplied by 2 to preserve the value of
            # 2*Nm during the simulation as ms does.
            symbiont_Nm,
            ((num_populations - 1) * 2),
        ),
    )
    np.fill_diagonal(migration_null, 0)
    if migration is None:
        migration = migration_null
    else:
        migration = migration_null * migration
    tree = ms.simulate(
        Ne=0.5,  # again, factor of 1/2 to preserve ms behavior
        num_replicates=num_replicates,
        migration_matrix=migration,
        population_configurations=population_config,
        mutation_rate=symbiont_theta / 2,  # factor of 1/2 to preserve ms beh.
        **kwargs
    )
    treesample = MetaSample(tree, populations, keep_populations)
    out = (
        np.zeros((num_replicates, len(stats) + len(params)))
        if not average_reps
        else np.zeros((1, len(stats) + len(params)))
    )
    if len(set(("fst_mean", "fst_sd")).intersection(set(stats))) > 0:
        fst_summ = treesample.fst(
            average_sites=True,
            average_reps=average_reps,
            summary=("mean", "sd"),
            h_opts=h_opts,
        )
    if np.isnan(np.array(list(fst_summ.values()))).any():
        breakpoint()
    for statidx, stat in enumerate(stats):
        if stat == "pi_h":
            out[:, statidx] = (
                treesample.pi(pi_method="h", h_opts=h_opts, **kwargs)
                if not average_reps
                else np.mean(
                    treesample.pi(pi_method="h", h_opts=h_opts, **kwargs)
                )
            )
        elif stat == "pi_nei":
            out[:, statidx] = (
                treesample.pi(pi_method="nei", **kwargs)
                if not average_reps
                else np.mean(treesample.pi(pi_method="nei", **kwargs))
            )
        elif stat == "pi_tajima":
            out[:, statidx] = (
                treesample.pi(pi_method="tajima", **kwargs)
                if not average_reps
                else np.mean(treesample.pi(pi_method="tajima", **kwargs))
            )
        elif stat == "num_sites":
            out[:, statidx] = (
                treesample.segsites()
                if not average_reps
                else np.mean(treesample.segsites())
            )
        elif stat == "theta_w":
            out[:, statidx] = (
                treesample.theta_w()
                if not average_reps
                else np.mean(treesample.theta_w())
            )
        elif stat == "fst_mean":
            out[:, statidx] = fst_summ["mean"]
        elif stat == "fst_sd":
            out[:, statidx] = fst_summ["sd"]
    out[:, -3] = eta
    out[:, -2] = tau
    out[:, -1] = rho
    return out if not average_reps else out[0]
