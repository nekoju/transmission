import collections.abc

import numpy as np
from rpy2 import robjects as robjects
from rpy2.rinterface import NULL
from rpy2.robjects import numpy2ri as numpy2ri, vectors as vectors
from rpy2.robjects.packages import importr


class Sample(object):

    """
    Coalescent simulation output representing a population class and methods.
    Input argument should be a np.ndarray with 2 dimensions detailed
    in __init__ or an msprime.TreeSequence object.

    Attributes:
        nchrom: the number of chromosomes sampled
        gtmatrix: the original matrix supplied, consisting of [0, 1]
        type: the input type, either "TreeSequence" or "ndarray"
        popdata: the original object passed to instantiate the class

    """

    def __init__(self, popdata):
        """
        Args:
            popdata (np.ndarray OR msprime.TreeSequence): A 2D np.ndarray
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
        if isinstance(popdata, collections.abc.Iterator):
            self.popdata = tuple(popdata)
        elif isinstance(popdata, list) or isinstance(popdata, tuple):
            self.popdata = popdata
        else:
            self.popdata = tuple(popdata for _ in (0,))
        if type(self.popdata[0]).__name__ == "TreeSequence":
            self.nchrom = self.popdata[0].num_samples
            self.type = "TreeSequence"
        else:
            self.nchrom = self.popdata[0].shape[0]
            self.type = "ndarray"
        self.populations = np.zeros(self.nchrom, dtype=int)
        self.npop = 1
        self.pop_sample_sizes = np.array(self.nchrom)
        self.num_replicates = len(self.popdata)
        self.keep_populations = set(self.populations)

    def __str__(self):
        print(
            "a {reps} replicate {nchrom} chromosome sample with "
            "{segsites} segregating sites".format(
                nchrom=self.nchrom,
                segsites=self.segsites(),
                reps=len(self.popdata),
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
                out.append(
                    rep.genotype_matrix().T[
                        np.isin(self.populations, tuple(self.keep_populations))
                    ]
                )
            else:
                out.append(rep)
        return out if len(out) > 1 else out[0]

    def h(
        self, average=False, bias=True, by_population=False, **polymorphic_opts
    ):

        """
        Calculate heterozygosity over sites in sample.
        Returns a list of np.ndarrays of dimensions npop X segsites,
        or, if average==True, a npop X num_replicates np.ndarray.
        Args:
            average (bool): whether to average H values for all sites
            bias (bool): whether to apply bias-correction to calculations
            by_population (bool): whether to calculate heterozygosity by
                population, or alternatively, as a metapopulation (H_T).
            polymorphic_opts (dict): Extra arguments for Sample.polymorphic().
        """

        # Check if multiple populations are provided or desired
        populations = (
            self.populations
            if by_population
            else np.zeros(self.nchrom, dtype=int)
        )
        # populations for portability to MetaSample
        # Each row in harray is a population; each column a snp.
        sample_sizes = (
            self.pop_sample_sizes[
                0,
                np.isin(
                    list(set(self.populations)), list(self.keep_populations)
                ),
            ]
            if by_population
            else np.sum(
                self.pop_sample_sizes.reshape(1, -1)[
                    0,
                    np.isin(
                        list(set(self.populations)),
                        list(self.keep_populations),
                    ),
                ]
            )
        ).reshape(1, -1)
        out = []
        for repidx, replicate in enumerate(self.popdata):
            if replicate.num_sites != 0:
                num_mutants = self.num_mutants(
                    by_population=by_population, popdata=replicate
                )[0]
                parray = np.true_divide(num_mutants, sample_sizes.T)
                harray = 2 * parray * (1 - parray)
                if bias:
                    harray = (
                        np.true_divide(sample_sizes, (sample_sizes - 1)).T
                        * harray
                    )
                if average:
                    out.append(np.mean(harray, axis=1))
                else:
                    out.append(harray)
            elif by_population:
                out.append(np.full(self.npop, 0.0))
            else:
                out.append(np.full(1, 0.0))
        if not average:
            return tuple(out)
        else:
            return np.array(out).T

    def num_mutants(self, by_population=False, popdata=None):
        """
        Returns the number of mutant alleles observed at each site.
        Used by other methods.

        Args:
            populations (numpy.ndarray): Ints corresponding to
                which population each chromosome belongs.
            popdata (TreeSequence object or np.ndarray): Single replicate for
                which to compute num_mutants.
        """
        # Allows computation of single replicate snps.
        popdata = self.popdata if not popdata else (popdata,)
        out = []
        for repidx, replicate in enumerate(popdata):
            num_mutants = (
                np.zeros(
                    (len(self.keep_populations), replicate.num_sites),
                    dtype=int,
                )
                if by_population
                else np.zeros((1, replicate.num_sites), dtype=int)
            )
            if self.type == "TreeSequence":
                for siteidx, site in enumerate(replicate.variants()):
                    if by_population:
                        for popidx, pop in enumerate(self.keep_populations):
                            num_mutants[popidx, siteidx] = np.count_nonzero(
                                site.genotypes[self.populations == pop]
                            )
                    else:
                        num_mutants[0, siteidx] = np.count_nonzero(
                            site.genotypes[
                                np.isin(
                                    self.populations,
                                    np.array(list(self.keep_populations)),
                                )
                            ]
                        )
            else:
                for pop in (
                    self.keep_populations
                    if by_population
                    else set(self.populations)
                ):
                    num_mutants[pop] = np.count_nonzero(
                        self.gtmatrix[np.nonzero(self.populations == pop)],
                        axis=0,
                    )
            out.append(num_mutants)
        # Return tuple if multiple replicates, otherwise single matrix.
        return tuple(out)

    def pi(self, pi_method="h", h_opts={}, polymorphic_opts={}, **kwargs):
        """
        Calculates different metrics of nucleotide diversity.

        Args:
            pi_method (str): what method to use to calculate pi.
            available methods are 'nei', 'tajima', and 'h'.
                nei: sum_{ij} x_i x_j pi_{ij} where x_y is the
                    frequency of the yth sequence and pi_{ij} is the
                    number of pairwise differences
                tajima: sum_{ij} k_{ij} / n choose 2
                    k_{ij} is the number of pairwise differences between the
                    ith and jth samples. This is the implementation used
                    originally in Tajima's D.
                h: Simply the sum of heterozygosities across all sites.
                    Takes arguments to self.h() as optional arguments,
                    namely 'bias'.
            h_opts (dict): Extra arguments for self.h() in the form of a
                kwargs dictionary.
            polymorphic_opts (dict): Extra arguments for Sample.polymorphic().
        """

        out = np.zeros((len(self.popdata),))
        if pi_method == "nei":
            for repidx, rep in enumerate(self.gtmatrix()):
                # Hash rows of genotype matrix to reduce time comparing.
                hashes = np.apply_along_axis(
                    lambda row: hash(tuple(row)), 1, rep
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
                        seqs[seqid]["p"] = (
                            np.count_nonzero(seqid_array == hashes)
                            / self.nchrom
                        )
                    # Associate sequences with hashes.
                    for i in np.arange(self.nchrom):
                        if seqid == hash(tuple(rep[i,])):
                            seqs[seqid]["seq"] = np.array(rep[i,])
                            break
                # Calculate nucleotide diversity.
                nucdiv = 0
                for i in seqs:
                    for j in seqs:
                        if i != j:
                            nucdiv += (
                                seqs[i]["p"]
                                * seqs[j]["p"]
                                * np.count_nonzero(
                                    seqs[i]["seq"] != seqs[j]["seq"]
                                )
                            )
                out[repidx] = nucdiv
        elif pi_method == "tajima":
            for repidx, rep in enumerate(self.gtmatrix()):
                k = 0
                # count number of pairwise differences for unique comparisons.
                for i in np.arange(self.nchrom - 1):
                    for j in np.arange(i, self.nchrom):
                        k += np.count_nonzero(rep[i,] != rep[j,])
                # Formula: \sum{k_ij} / (nchrom choose 2)
                out[repidx] = k / ((self.nchrom - 1) * self.nchrom / 2)
        elif pi_method == "h":
            out = np.array(
                [np.sum(x) for x in self.h(**h_opts, **polymorphic_opts)]
            )
        else:
            raise NameError("Unsupported method, {}".format(pi_method))
        return out if out.size > 1 else out[0]

    def polymorphic(
        self, by_population=False, threshold=0, **polymorphic_opts
    ):
        """
        Returns polymorphism attributes as dict or value.
        Attributes are number of polymorphic sites and index positions of
        polymorphic sites.

        Args:
            by_population (bool): Whether to calculate by population.
            threshold (int): Number of derived alleles
                above which polymorphism is counted, e.g.,
                threshold=1 excludes singletons.
        Returns:
            A tuple of dicts equal in length to the number of replicates.
            Each dict contains "which", the index positions of polymorphic
            loci in the form np.ndarray([pop]), np.ndarray([site]),
            and "num", which contains the number of polymorphic sites in
            each population.
        """
        for repidx, rep in enumerate(
            self.num_mutants(by_population=by_population)
        ):
            num = np.count_nonzero(rep > threshold, axis=1)
            # np.nonzero returns an array of *coordinates*.
            # Transpose for ordered pairs on rows.
            which = np.nonzero(rep > threshold)
            out = {"which": which, "num": num}
            yield out

    def segsites(self):
        out = np.zeros(len(self.popdata), dtype=int)
        for repidx, replicate in enumerate(self.popdata):
            if self.type == "TreeSequence":
                out[repidx] = replicate.num_sites
            else:
                out[repidx] = replicate.shape[1]
        return out[0] if out.size == 1 else out

    def theta_w(self, by_population=False, polymorphic_opts={}):
        """
        Calculate the Watterson estimator.

        Args:
            by_population (bool): Whether to compute theta_w per-population.
            polymorphic_opts (dict): Extra arguments for Sample.polymorphic().
        Returns:
            If by_population: npop X num_replicates np.ndarray
            If not by_population: 1 X num_replicates np.ndarray
        """
        # Check if multiple populations are provided or desired
        populations = (
            self.populations
            if by_population
            else np.zeros((self.nchrom), dtype=int)
        )
        # populations for portability to MetaSample
        nchrom = np.bincount(
            populations[np.isin(populations, tuple(self.keep_populations))]
        )
        nchrom = nchrom[np.nonzero(nchrom)]
        harmonic_num = np.vectorize(lambda x: np.sum(1 / np.arange(1, x)))(
            nchrom
        ).reshape(len(self.keep_populations) if by_population else 1, 1)
        num_polymorphic = np.zeros(
            (
                len(self.keep_populations) if by_population else 1,
                self.num_replicates,
            )
        )
        for repidx, rep in enumerate(
            self.polymorphic(by_population, **polymorphic_opts)
        ):
            num_polymorphic[:, repidx] = rep["num"]

        return num_polymorphic / harmonic_num

    def filter_variants(self):
        """
        Returns a generator yielding only sites that are polymorphic across
        all populations.
        """

        def polymorphic_sites(variants, polymorphic):
            """
            For a replicate, filters out sites that are not polymorphic.
            
            Args:
                variants (msprime TreeSequence)
                polymorphic (np.ndarray): The indices of polymorphic sites.
            """
            for siteidx, site in enumerate(variants):
                if siteidx in polymorphic:
                    site.genotypes = site.genotypes[
                        np.isin(self.populations, self.keep_populations)
                    ]
                    yield site

        polymorphic = (x["which"][1] for x in self.polymorphic())
        for rep in self.popdata:
            yield polymorphic_sites(rep.variants(), next(polymorphic))


class MetaSample(Sample):

    """

    Class representing a metapopulation sample and associated methods.

    Input argument should be a nchrom X segsites np.ndarray with an
    array listing to which population each chromosome belongs,
    or an msprime.TreeSequence object.

    Attributes:
        npop (int): The number of populations in the sample.
        pop_sample_sizes (np.ndarray): The number of individuals sampled from
            each population.
        populations (np.ndarray): A 1-dimensional np.ndarray of ints
            specifying to which population each sample belongs.
        keep_populations (np.ndarray): Specifys which populations to keep
            in analysis. Allows simulation of sampling only a subset of
            actual populations. Defaults to None, which retains all
            populations.

    """

    def __init__(
        self, popdata, populations, keep_populations=None, force_meta=False
    ):
        """
        Args:
            popdata (np.ndarray OR msprime.TreeSequence): A 2D np.ndarray
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
            populations (np.ndarray): A 1D np.ndarray of ints specifying to
                which population each sample belongs.
            keep_populations (np.ndarray): Specifys which populations to keep
                in analysis. Allows simulation of sampling only a subset of
                actual populations. Defaults to None, which retains all
                populations.
            force_meta (bool): Whether to force a single population to become
                a metapopulation rather than coercion to Population class.
        """
        super().__init__(popdata)
        if len(set(populations)) == 1 or (
            self.type == "TreeSequence"
            and self.popdata[0].num_populations == 1
            and not force_meta
        ):
            raise Exception(
                "Only 1 population provided. "
                "Use force_meta=True for MetaSample or use Sample."
            )
        else:
            self.npop = (
                self.popdata[0].num_populations
                if self.type == "TreeSequence"
                else len(set(populations))
            )
        self.pop_sample_sizes = np.array(
            [
                [
                    np.count_nonzero(np.full(self.nchrom, x) == populations)
                    for x in set(populations)
                ]
            ]
        )
        self.populations = populations
        self.keep_populations = (
            set(populations) if keep_populations is None else keep_populations
        )
        self.nchrom = np.count_nonzero(
            np.isin(self.populations, tuple(self.keep_populations))
        )

    def fst(
        self,
        fst_method="gst",
        summary="mean",
        average_reps=False,
        average_sites=True,
        h_opts={},
        polymorphic_opts={},
    ):
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
            average_reps (bool): Whether to average final fst values. If true,
                returns 1 fst value, otherwise 1 value for each replicate in an
                array. Redundant if not average_sites.
            h_opts (dict): Extra arguments for Sample.h(), in the form of a
                kwargs dictionary.
            polymorphic_opts (dict): Extra arguments for Sample.polymorphic().
        """
        if isinstance(summary, str):
            summary = (summary,)
        if fst_method == "gst":
            ind = np.where(self.segsites() != 0)[0]
            h = self.h(by_population=True, **h_opts, **polymorphic_opts)
            h_by_site = (h[i] for i in ind)
            hs = (
                np.average(
                    x,
                    axis=0,
                    weights=self.pop_sample_sizes[
                        0, tuple(self.keep_populations)
                    ],
                )
                for x in h_by_site
            )
            ht = self.h(**h_opts, **polymorphic_opts)
            ht = (x.reshape(-1) for x in ht)
            fst = tuple((1 - x / y) for x, y in zip(hs, ht))
            if average_sites:
                stats = []
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
                if average_reps:
                    try:
                        stats = [
                            np.average(x, weights=self.segsites()[ind])
                            for x in stats
                        ]
                    except Exception:
                        stats = np.array([0.0 for _ in stats])
                return dict(zip(summary, stats))
            else:
                # Clean up NaNs
                return tuple(x[~np.isnan(x)] for x in fst)
        else:
            raise Exception("invalid method {}".format(fst_method))


class Abc(object):

    """
    Interface to R package abc for estimating tau and rho posteriors.
    """

    def __init__(
        self,
        target,
        param,
        sumstat,
        tol=0.1,
        method="loclinear",
        transf={"tau": "logit", "rho": "logit", "eta": None},
        logit_bounds={"tau": (0, 1), "rho": (0, 1), "eta": (1, 1)},
        **kwargs
    ):
        """
        target(np.ndarray): 0 x nstat array of calculated summary
            statistics from observed sample. If a structured array is provided,
            column names will be returned.
        param (np.ndarray): num_iterations x 3 array of eta, tau, rho values
            simulated from their prior distributions. If a structured array is
            provided, column names will be returned.
        sumstat(np.ndarray): num_iterations x num_simulations array of summary
            statistics from simulated data. If a structured array is provided,
            column names will be returned.
        tol (float): Accepted proportion of sumstats accepted.
            tol * sumstat.shape[0] must be >= 2
        method (str): The ABC algorithm to use. May take "loclinear",
            "neuralnet", "rejection", or "ridge" as values. "loclinear" is
            chosen as the default as rejection does not allow transformation
            of parameters.
        transf (str): A dictionary of parameter: transformation pairs to use
            for parameter values. Each may take "log", "logit", or None.
        logit_bounds (np.ndarray): A dictionary with elements "tau", "rho",
            and "eta" each representing a tuple of lower and upper bounds
            for logit transformation to use if transf="logit" is specified.
        **kwargs: Additional arguments for r function 'abc'. Must be formatted
            for passing to R as described in the documentation for rpy2.
        """

        if tol * sumstat.shape[0] <= 1:
            raise ValueError("tol * sumstat.shape[0] must be >= 2")
        logit_bounds_matrix = (
            np.array([logit_bounds[x] for x in param.dtype.names])
            .view(type=np.ndarray)
            .reshape((len(param.dtype.names), -1))
        )
        # Change transformation from None to "none" for R.
        for key, value in transf.items():
            if not value:
                transf[key] = "none"
        if param.dtype.names:
            transf_list = [transf[x] for x in param.dtype.names]
        else:
            raise AttributeError("No names for parameter array")

        numpy2ri.activate()
        importr("abc")
        importr("base")
        self.abc = robjects.r.abc(
            target=target,
            param=Abc.rmatrix(param),
            sumstat=Abc.rmatrix(sumstat),
            tol=vectors.FloatVector([tol]),
            method=vectors.StrVector([method]),
            transf=(
                vectors.StrVector(["none"])
                if not transf
                else vectors.StrVector(transf_list)
            ),
            logit_bounds=Abc.rmatrix(logit_bounds_matrix),
            **kwargs
        )
        if method != "rejection":
            self.adj_values = np.core.records.fromarrays(
                np.array(self.abc.rx2("adj.values")).T, names=param.dtype.names
            )
            self.weights = np.array(self.abc.rx2("weights"))
            self.residuals = np.array(self.abc.rx2("residuals"))
        if method == "neuralnet":
            self.lambda_0 = np.array(self.abc.rx2("lambda"))
        if self.abc.rx2("na.action") != NULL:
            self.na_action = robjects.vectors.BoolVector(
                self.abc.rx2("na.action")
            )
        if not self.abc.rx2("region") != NULL:
            self.region = robjects.vectors.BoolVector(self.abc.rx2("region"))
        self.unadj_values = np.core.records.fromarrays(
            np.array(self.abc.rx2("unadj.values")).T, names=param.dtype.names
        )
        self.ss = np.core.records.fromarrays(
            np.array(self.abc.rx2("ss")).T, names=sumstat.dtype.names
        )
        self.dist = np.array(self.abc.rx2("dist"))
        self.call = self.abc.rx2("call")
        self.transf = np.array(self.abc.rx2("transf"))
        self.logit_bounds = np.array(self.abc.rx2("logit.bounds"))
        self.method = self.abc.rx2("method")
        self.kernel = self.abc.rx2("kernel")
        self.numparam = np.array(self.abc.rx2("numparam"))
        self.numstat = np.array(self.abc.rx2("numstat"))
        self.aic = np.array(self.abc.rx2("aic"))
        self.bic = np.array(self.abc.rx2("bic"))
        self.names = {
            "paramater_names": self.abc.rx("names").rx2("parameter.names"),
            "statistics_names": self.abc.rx("names").rx2("statistics.names"),
        }
        numpy2ri.deactivate()

    def summary(self):
        return robjects.r.summary(self.abc, print=False)

    @staticmethod
    def rmatrix(rec_array):
        """
        Return an r matrix, preserving column and row names, from a
            numpy record array.

        Args:
            rec_array (np.ndarray): A n X m numpy structured array.
        """

        if isinstance(rec_array, np.recarray):
            # Extract data from record array, i.e. make into standard
            # ndarray.
            rec_array_in = rec_array.view(np.float64).reshape(
                (rec_array.shape + (-1,) if rec_array.shape else (1, -1))
            )
        else:
            rec_array_in = rec_array
        if rec_array.shape:
            out = robjects.r.matrix(
                rec_array_in,
                nrow=rec_array.shape[0],
                ncol=(
                    len(rec_array.dtype.names)
                    if rec_array.dtype.names
                    else rec_array.shape[1]
                ),
                dimnames=robjects.r.list(
                    NULL,
                    (
                        vectors.StrVector(rec_array.dtype.names)
                        if rec_array.dtype.names
                        else NULL
                    ),
                ),
            )
        else:
            out = robjects.vectors.FloatVector(rec_array_in)
            out.names = robjects.vectors.StrVector(rec_array.dtype.names)
        return out
