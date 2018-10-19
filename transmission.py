from __future__ import print_function, division
import collections
import functools
from multiprocessing import Pool
import pdb

import numpy as np
import msprime as ms


class Abc(object):
    """
    Interface to R package abc for estimating tau and rho posteriors.
    """

    def __init__(self, stats_observed, stats_simulated, params=None):
        """
        stats_observed(np.ndarray): 1 x nstat array of calculated summary
            statistics from observed sample.
        stats_simulated(np.ndarray): num_iterations x nstat array of summary
            statistics from simulated data.
        params (np.ndarray): num_iterations x 2 array of tau, rho values
            simulated from their prior distributions. If None, will be read
            from stats_simulated["rho"] and stats_simulated["rho"]
        """


class Sample(object):
    """
    Coalescent simulation output representing a population class and methods.
    Input argument should be a np.ndarray with 2 dimensions detailed
    in __init__ or an msprime.TreeSequence object.

    Attributes:
        segsites: the number of segregating sites in the sample
        nchrom: the number of chromosomes sampled
        gtmatrix: the original matrix supplied, consisting of [0, 1]
        type: the input type, either "TreeSequence" or "ndarray"
        popdata: the original object passed to instantiate the class
    """

    def __init__(self, popdata):
        """
        Args:
            popdata (np.ndarray OR msprime.TreeSequence): A 2 x 2 np.ndarray
                with individual chromosomes
                represented by rows and snp positions represented by columns.
                Alternatively, an object of class "TreeSequence" can be
                provided and the relevant attributes will be inferred.
                Note that this it the transpose of
                msprime.TreeSequence.genotype_matrix(). This accomodates
                generators output by msprime.simulate() num_replicates
                argument and accepts such an iterable  while returning
                a tuple of "TreeSequence" or "np.ndarray"
                objects.
        """

        # Make self.popdata iterable for consistency with multiple replicates.
        if isinstance(popdata, collections.Iterator):
            self.popdata = tuple(popdata)
        elif (isinstance(popdata, list) or isinstance(popdata, tuple)):
            self.popdata = popdata
        else:
            self.popdata = tuple(popdata for _ in (0, ))
        if type(self.popdata[0]).__name__ == "TreeSequence":
            self.nchrom = self.popdata[0].num_samples
            self.type = "TreeSequence"
        else:
            self.nchrom = self.popdata[0].shape[0]
            self.type = "ndarray"
        self.populations = np.zeros(self.nchrom)
        self.pop_sample_sizes = np.array(self.nchrom)

    def __str__(self):
        print(
            "a {reps} replicate {nchrom} chromosome sample with "
            "{segsites} segregating sites".format(
                nchrom=self.nchrom,
                segsites=self.segsites(),
                reps=len(self.popdata)
                )
            )

    def gtmatrix(self):
        """
        Sample.gtmatrix() returns a
        self.nchrom X self.sites matrix of sample genotypes.
        """
        out = []
        for rep in self.popdata:
            if self.type == "TreeSequence":
                out.append(rep.genotype_matrix().T)
            else:
                out.append(rep)
        return out if len(out) > 1 else out[0]

    def h(self, replace=False, average=False, bias=True,
          by_population=False, **kwargs):
        """
        Calculate heterozygosity over sites in sample.
        Returns a list of np.ndarrays of dimensions npop X segsites,
        or, if average==True, a npop X num_replicates np.ndarray.
        Args:
            replace (bool): calculate heterozygosity with replacement or not
            average (bool): whether to average H values for all sites
            bias (bool): whether to apply bias-correction to calculations
            by_population (bool): whether to calculate heterozygosity by
                population, or alternatively, as a metapopulation (H_T).
        """

        # Check if multiple populations are provided or desired
        populations = (self.populations
                       if by_population
                       else np.zeros(self.nchrom, dtype=int))
        # populations for portability to MetaSample
        out = []
        # Each row in harray is a population; each column a snp.
        for repidx, replicate in enumerate(self.popdata):
            if replicate.num_sites != 0:
                num_mutants = self.num_mutants(populations, replicate)
                sample_sizes = (self.pop_sample_sizes
                                if by_population
                                else np.array([self.nchrom]))
                if replace:
                    parray = np.true_divide(num_mutants, sample_sizes.T)
                    harray = 2 * parray * (1 - parray)
                # heterozygosity should probably be calculated without
                # replacement especially with small samples
                else:
                    harray = (
                        2 * np.true_divide(
                            num_mutants, sample_sizes.T
                            ) *
                        np.true_divide(
                            (sample_sizes.T - num_mutants),
                            (sample_sizes - 1).T
                            )
                        )
                if bias:
                    harray = (np.true_divide(
                        sample_sizes, (sample_sizes - 1)
                        ).T * harray)
                if average:
                    out.append(np.mean(harray, axis=1))
                else:
                    out.append(harray)
            elif by_population:
                out.append(np.full(self.npop, float("nan")))
            else:
                out.append(np.full(1, float('nan')))
        if not average:
            return tuple(out)
        elif not by_population:
            return np.array(out).T[0]
        else:
            return np.array(out).T

    def num_mutants(self, populations, popdata=None):
        """
        Returns the number of mutant alleles observed at each site.
        Used by other methods.

        Args:
            populations (list, tuple, or np.ndarray): Ints corresponding to
                which population each chromosome belongs.
            popdata (TreeSequence object or np.ndarray): Single replicate for
                which to compute num_mutants.
        """
        popdata = (self.popdata
                   if not popdata
                   else tuple(popdata for _ in (0, )))
        popset = set(populations)
        out = []
        for repidx, replicate in enumerate(popdata):
            num_mutants = np.zeros(
                (len(popset), replicate.num_sites),
                dtype=int
                )
            if self.type == "TreeSequence":
                for siteidx, site in enumerate(replicate.variants()):
                    for pop in popset:
                        num_mutants[pop, siteidx] = np.count_nonzero(
                            site.genotypes[populations == pop]
                            )
            else:
                for pop in popset:
                    num_mutants[pop] = np.count_nonzero(
                        self.gtmatrix[np.nonzero(populations == pop)], axis=0
                        )
            out.append(num_mutants)
        # Return tuple if multiple replicates, otherwise single matrix.
        if len(popdata) > 1:
            return tuple(out)
        else:
            return out[0]

    def pi(self, pi_method="nei", *args, **kwargs):
        """
        Calculates different metrics of nucleotide diversity.

        Args:
            pi_method (str): what method to use to calculate pi.
            available methods are 'nei', 'tajima', and 'h'.
                nei: \sum_{ij} x_i x_j pi_{ij} where x_y is the
                    frequency of the yth sequence and pi_{ij} is the
                    number of pairwise differences
                tajima: \sum_{ij} k_{ij} / n choose 2
                    k_{ij} is the number of pairwise differences between the
                    ith and jth samples. This is the implementation used
                    originally in Tajima's D.
                h: Simply the sum of heterozygosities across all sites.
                    Takes arguments to self.h() as optional arguments,
                    namely 'bias' and 'replace'.
        """

        if pi_method is "nei":
            # Hash rows of genotype matrix to reduce time comparing.
            hashes = np.apply_along_axis(
                    lambda row: hash(tuple(row)), 1, self.gtmatrix()
                    )
            seqs = dict.fromkeys(set(hashes))
            # Loop over seqs keys to calculate sequence frequncies
            # as well as create dict items for sequences themselves.
            for seqid in seqs:
                seqs[seqid] = {}
                seqid_array = np.full(self.nchrom, seqid)
                # Get frequencies.
                if self.nchrom == 0:
                    seqs[seqid]["p"] = 0.0
                else:
                    seqs[seqid]["p"] = (np.count_nonzero(
                                        seqid_array == hashes) /
                                        self.nchrom)
                # Associate sequences with hashes.
                for i in np.arange(self.nchrom):
                    if seqid == hash(tuple(self.gtmatrix()[i, ])):
                        seqs[seqid]["seq"] = np.array(self.gtmatrix()[i, ])
                        break
            # Calculate nucleotide diversity.
            nucdiv = 0
            for i in seqs:
                for j in seqs:
                    if i != j:
                        nucdiv += (seqs[i]['p'] * seqs[j]['p'] *
                                   np.count_nonzero(
                                       seqs[i]["seq"] != seqs[j]["seq"]
                                       )
                                   )
            return nucdiv
        elif pi_method is "tajima":
            k = 0
            # count number of pairwise differences for all unique comparisons.
            gtmatrix = self.gtmatrix()
            for i in np.arange(self.nchrom - 1):
                for j in np.arange(i, self.nchrom):
                    k += np.count_nonzero(
                        gtmatrix[i, ] != gtmatrix[j, ]
                        )
            # Formula: \sum{k_ij} / (nchrom choose 2)
            return k / ((self.nchrom - 1) * self.nchrom / 2)
        elif pi_method is "h":
            return np.sum(self.h(**kwargs))
        else:
            print("Unsupported method, {}".format(pi_method))

    def polymorphic(self, threshold=0, output=("num", "which")):
        """
        Returns polymorphism attributes as dict or value.
        Attributes are number of polymorphic sites and index positions of
        polymorphic sites.

        Args:
            threshold (int): Number of derived alleles
                above which polymorphism is counted, e.g.,
                threshold=1 excludes singletons.
            output (tuple or str): Whether to output number of
                polymorphic sites or index positions of polymorphic sites,
                or both. Possible values are "num" returning number of sites,
                "which" returning indices,
                or a tuple of both (default) returning a dict of both.
        """
        num_mutants = self.num_mutants()
        polymorphic = np.array(
            [np.full(self.segsites, threshold) < num_mutants,
             num_mutants < np.full(self.segsites, self.nchrom - threshold)]
            ).all(axis=0)
        if "which" in output:
            which = np.nonzero(polymorphic)
        if "num" in output:
            num = np.count_nonzero(polymorphic)
        if type(output) == tuple:
            return {"which": which, "num": num}
        elif output == "num":
            return num
        elif output == "which":
            return which
        else:
            print(
                "Invalid output, {output} specified."
                "Specify 'num' or 'which'.".format(output=output)
                )

    def segsites(self):
        out = np.zeros(len(self.popdata), dtype=int)
        for repidx, replicate in enumerate(self.popdata):
            if self.type == "TreeSequence":
                out[repidx] = replicate.num_sites
            else:
                out[repidx] = replicate.shape[1]
        return out


class MetaSample(Sample):
    """
    Class representing a metapopulation sample and associated methods.
    Input argument should be a nchrom X segsites np.ndarray with an
    array listing to which population each chromosome belongs,
    or an msprime.TreeSequence object.
    """

    def __init__(self, popdata, populations, force_meta=False):
        super(MetaSample, self).__init__(popdata)
        if (len(set(populations)) == 1
            or (self.type == "TreeSequence"
                and self.popdata[0].num_populations == 1
                and not force_meta)):
            raise Exception(
                "Only 1 population provided. "
                "Use force_meta=True for MetaSample or use Sample."
                )
        else:
            self.npop = (self.popdata[0].num_populations
                         if self.type == "TreeSequence"
                         else len(set(populations)))
        self.pop_sample_sizes = np.array(
            [[np.count_nonzero(np.full(self.nchrom, x) == populations)
              for x in set(populations)]]
            )
        self.populations = populations

    def fst(self, fst_method="gst",
            summary="mean", average_final=False, average_sites=True,
            **kwargs):
        """
        Returns specified fst statistic.
        Args:
            fst_method (str): Only "gst" is supported right now. This returns
            the classic (ht - hs) / ht statistic.
            average_sites (bool): Whether to average the heterozygosities
                across sites when calculating. If false, will return a
                num_replicates list of self.segsites()-long arrays.
            summary (str or tuple): The summary statistics returned for
                replicates. Can take arguments "mean" or "sd" for standard
                deviation. Returns a dictionary of the supplied statistics.
                Redundant if not average_sites.
            average_final (bool): Whether to average final fst values. If true,
                returns 1 fst value, otherwise 1 value for each replicate in an
                array. Redundant if not average_sites.
            **kwargs (): Optional arguments for self.h().
        """
        if isinstance(summary, str):
            summary = (summary, )
        if fst_method == "gst":
            ind = np.where(self.segsites() != 0)[0]
            h_by_site = tuple([self.h(by_population=True, **kwargs)[i]
                               for i in ind])
            hs = tuple(
                    [np.average(
                     x, axis=0, weights=self.pop_sample_sizes[0])
                     for x in h_by_site]
                )
            ht = tuple(
                    [x[0] for x in tuple([self.h(**kwargs)[i]
                     for i in ind])]
                )
            fst = tuple([(1 - np.true_divide(x, y))
                         for x, y
                         in zip(hs, ht)])
            if average_sites:
                stats = []
                if type(summary) == str:
                    summary = tuple(summary)
                if "mean" in summary:
                    if self.segsites().any():
                        stats.append([np.average(x) for x in fst])
                    else:
                        stats.append(0.0)
                if "sd" in summary:
                    if self.segsites().any():
                        stats.append([np.std(x) for x in fst])
                    else:
                        stats.append(0.0)
                if average_final:
                    try:
                        stats = [np.average(
                                 x, weights=self.segsites()[ind])
                                 for x in stats]
                    except Exception:
                        stats = np.array([0.0 for _ in stats])
                return dict(zip(summary, stats))
            else:
                return fst
        else:
            raise Exception("invalid method {}".format(fst_method))


def beta_nonst(alpha, beta, a=0, b=1, n=1):
    """
    Draws instances of a random variable from a nonstandard beta distribution
    over the interval [a, b]. By default returns standard beta variables.
    Args:
        alpha (float): The first beta shape parameter. A nonnegative float.
        beta (float): The second beta shape parameter.
        a (float): The beginning of the interval of the distribution.
        b (float): The end of the interval of the distribution.
        n (int): The number of observations to be drawn.
    """

    return np.random.beta(alpha, beta, n) * (b - a) + a


def ms_simulate(nchrom, num_populations, host_theta, M, num_simulations,
                stats=("fst_mean", "fst_sd", "pi_h"),
                # 3.766e-3
                prior_params={"sigma": 1., "tau": (1, 1), "rho": (1, 1)},
                nsamp_populations=None, nrep=1, num_cores="auto",
                prior_seed=None, **kwargs):
    """
    Generate random sample summary using msprime for the specified prior
    distributions of tau (vertical transmission rate) and rho (sex ratio).

    From a supplied array of prior (hyper)parameters, ms_simulate generates a
    num_simulations-tall array of parameters with which to estimate the
    posterior distributions of tau and rho using abc methods. It then draws
    nrep coalescent samples for each parameter combination and outputs a
    summary table consisting of mean fst, fst standard deviation, pi, and the
    original simulated parameter values for each metasimulation.

    Args:
        nchrom (int): The number of chromosomes to sample from each population
        num_populations (int): The number of populations to include in the
            migration matrix.
        host_theta (float): The host's haploid theta (2 * Ne * mu * L)
            estimated from mitochondrial or nuclear data.
        M (float): The symmetric migration parameter (2 * Ne * m), estimated
            from host mitochondrial or nuclear fst.
        num_simulations (int): The number of tau-rho pairs of parameters
            to draw from their priors, and hence the number of metasimulations
            to include in the output.
        prior_params (dict): A dict containing tuples specifiying the prior
            distribution parameters for sigma, tau, and rho. That is, the
            mutation rate multiplier, vertical transmission frequency, and
            sex ratio. Optionally one may provide a scalar value for sigma to
            fix the mutation rate multiplier.
            Default mutation rate multiplier estimate from values found in:
            Haag-Liautard et al. 2008. Direct Estimation of the mitochondrial
            DNA mutation rate in Drosophila melanogaster. PLOS Biology.
            Sproufske et al. 2018. High mutation rates limit evolutionary
            adaptation in Escherichia coli. PLOS Genetics.
        nrep (int): Number of msprime replicates to run for each simulated
            parameter pair and the number of simulations in each
            metasimulation.
        num_cores (int or str): The number of cores to use for computation.
            "auto" will automatically detect using multiprocessing.cpu_count().
            None will use a single thread not routed through
            multiprocessing.Pool, primarily for debugging and profiling.
            An int value will specify the number of threads to use.
        prior_seed (int): The seed used to draw samples for tau and rho.
            Setting a seed will allow repeatability of results.
        **kwargs (): Extra arguments for ms.simulate(), Sample.pi(), and
            Sample.fst().
    """

    # Number of output statistics plus parameters tau and rho.

    populations = np.repeat(np.arange(num_populations), nchrom)
    population_config = tuple([ms.PopulationConfiguration(nchrom)
                               for _ in np.arange(nrep)])
    nsamp_populations = (num_populations
                         if not nsamp_populations
                         else nsamp_populations)
    migration = np.full((num_populations, num_populations),
                        np.true_divide(M, ((num_populations - 1) * 4)))
    for i in np.arange(num_populations):
        migration[i, i] = 0
    if prior_seed:
        np.random.seed(prior_seed)
    if isinstance(prior_params["sigma"], float):
        sigma = np.full((num_simulations, ), prior_params["sigma"])
    else:
        sigma = np.exp(
            np.random.normal(
                prior_params["sigma"][0],
                prior_params["sigma"][1],
                num_simulations
                )
            )
    tau = beta_nonst(prior_params["tau"][0], prior_params["tau"][1],
                     n=num_simulations)
    rho = beta_nonst(prior_params["rho"][0], prior_params["rho"][1], a=0, b=2,
                     n=num_simulations)
    theta = host_theta * sigma * np.true_divide(
        rho,
        tau ** 2 * (3 - 2 * tau) * (2 - rho) + rho
        )
    params = zip(theta, sigma, tau, rho)
    simpartial = functools.partial(
        sim, migration=migration,
        population_config=population_config,
        populations=populations, stats=stats, **kwargs
        )
    structure = {"names": tuple(list(stats) + prior_params.keys()),
                 "formats": tuple(
                     np.repeat("f8", len(stats) + len(prior_params))
                     )}
    if num_cores:
        if num_cores == "auto":
            num_cores = None
        pool = Pool(processes=num_cores)
        out = np.array(pool.map(simpartial, params))
        out = np.core.records.fromarrays(out.T, dtype=structure)
    else:
        out = np.apply_along_axis(simpartial, 1, params)
        out = np.core.records.fromarrays(out.T, dtype=structure)
    return out


def sim(params, migration, population_config, populations, stats, **kwargs):
    """
    Runs actual simulation with ms. Intended as helper for ms_simulate().

    At top level for picklability (for multiprocessing).

    Args:
        params (tuple): theta, tau, rho for simulation.
        migration (np.ndarray): The migration matrix. Off-diagonals are
            calculated as M / ((d - 1) * 4).
        popoulation_config (list): List of population configurations for
            msprime.
        populations (np.ndarray): A nchrom np.ndarray indicating to which
            population each chromosome belongs.
        **kwargs (): Extra arguments for msprime.simulate(), Sample.pi(),
            and Sample.h().
    """
    theta, sigma, tau, rho = params
    tree = ms.simulate(
                migration_matrix=migration,
                population_configurations=population_config,
                mutation_rate=theta / 4,
                **kwargs
            )
    treesample = MetaSample(tree, populations)
    # - 1 for calculated theta
    out = np.zeros((len(stats) + len(params) - 1, ))
    for i, stat in enumerate(stats):
        if len(set(("fst_mean", "fst_sd")).intersection(set(stats))) > 0:
            fst_summ = treesample.fst(average_sites=True, average_final=True,
                                      summary=("mean", "sd"), **kwargs)
        if stat == "pi_h":
            out[i] = treesample.pi(method="h", **kwargs)
        elif stat == "pi_nei":
            out[i] = treesample.pi(method="nei", **kwargs)
        elif stat == "pi_tajima":
            out[i] = treesample.pi(method="tajima", **kwargs)
        elif stat == "fst_mean":
            out[i] = fst_summ["mean"]
        elif stat == "fst_sd":
            out[i] = fst_summ["sd"]
    out[-3] = sigma
    out[-2] = tau
    out[-1] = rho
    return out


def main():
    # d = 2
    # n = 4
    # population_configurations = [ms.PopulationConfiguration(n)
    #                              for _ in range(d)]
    # migration = np.full((d, d), 2.5)
    # for i in range(d):
    #     for j in range(d):
    #         if i == j:
    #             migration[i, j] = 0
    # test = tuple(ms.simulate(
    #     population_configurations=population_configurations,
    #     migration_matrix=migration,
    #     # reset to 1 / 4
    #     # no_snps, set to 0.25 / 4
    #     mutation_rate=1 / 4,
    #     num_replicates=2,
    #     # reset to 3
    #     # 6 gives one popln with no snps and one with
    #     random_seed=3
    #     ))
    # # a = (testsample.h(average=False, by_population=True))
    # b = testsample.h()
    print(ms_simulate(4, 2, 1, 1, 2, nsamp_populations=None, nrep=2,
                      random_seed=3, num_cores=None))
    print(ms_simulate(4, 2, 1, 1, 2, nsamp_populations=None, nrep=2,
                      random_seed=3, num_cores=4)["fst_mean"])


if __name__ == "__main__":
    main()
