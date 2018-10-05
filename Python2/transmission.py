from __future__ import print_function, division
import collections
import numpy as np
import msprime as ms
import pdb


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
        return out

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
        # popset = set(populations)
        out = []
        # Each row in harray is a population; each column a snp.
        for repidx, replicate in enumerate(self.popdata):
            # harray = (np.zeros(len(popset), dtype=int)
            #           if average
            #           else np.zeros(
            #             (len(popset), self.segsites()[repidx]),
            #             dtype=int)
            #           )
            # k=sum of number of mutants per site, klist=list of those sums
            num_mutants = self.num_mutants(populations, replicate)
            sample_sizes = (self.pop_sample_sizes
                            if by_population
                            else np.array([self.nchrom]))
            if replace:
                parray = np.true_divide(num_mutants, sample_sizes.T)
                harray = 2 * parray * (1 - parray)
            # heterozygosity should probably be calculated without replacement,
            # especially with small samples
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
                harray = (np.true_divide(sample_sizes, (sample_sizes - 1)).T
                          * harray)
            if average:
                out.append(np.mean(harray, axis=1))
            else:
                out.append(harray)
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

    def pi(self, method="nei", *args, **kwargs):
        """
        Calculates different metrics of nucleotide diversity.

        Args:
            method (str): what method to use to calculate pi.
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

        if method is "nei":
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
        elif method is "tajima":
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
        elif method is "h":
            return np.sum(self.h(**kwargs))
        else:
            print("Unsupported method, {}".format(method))

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

    def fst(self, method="gst",
            average_h=True, average_final=False, average_sites=True,
            **kwargs):
        """
        Returns specified fst statistic.
        Args:
            method (str): Only "gst" is supported right now. This returns the
            classic (ht - hs) / ht statistic.
            average_sites (bool): Whether to average the heterozygosities
                across sites when calculating. If false, will return a
                num_replicates list of self.segsites()-long arrays.
            average_h (bool): Whether to average heterozygisities before taking
                their ration to return fst. If true, will return 1 fst value.
                Redundant if average_sites.
            average_final (bool): Whether to average final fst values. If true,
                returns 1 fst value, otherwise 1 value for each replicate in an
                array. Redundant if average_h.
            **kwargs (): Optional arguments for self.h().
        """
        if method == "gst":
            if average_sites:
                h_by_rep = self.h(average=True, by_population=True,
                                  **kwargs)
                hs_by_rep = np.average(h_by_rep, axis=0,
                                       weights=self.pop_sample_sizes[0])
                ht_by_rep = self.h(average=True, **kwargs)
                if average_h:
                    hs = np.average(hs_by_rep, weights=self.segsites())
                    ht = np.average(ht_by_rep, weights=self.segsites())
                    return 1 - hs / ht
                else:
                    fst = 1 - np.true_divide(hs_by_rep, ht_by_rep)
                    if average_final:
                        return np.average(fst, weights=self.segsites())
                    else:
                        return fst
            else:
                h_by_site = self.h(by_population=True, **kwargs)
                hs = tuple([np.average(
                                x, axis=0, weights=self.pop_sample_sizes[0]
                                )
                            for x in h_by_site])
                ht = tuple(
                        [x[0] for x in self.h(**kwargs)]
                    )
                # returns num_replicates tuple of self.segsites() np.arrays
                return tuple([(1 - x / y) for x, y in zip(hs, ht)])
        else:
            raise Exception("invalid method {}".format(method))


def main():
    d = 2
    n = 4
    population_configurations = [ms.PopulationConfiguration(n)
                                 for _ in range(d)]
    migration = np.full((d, d), 2.5)
    for i in range(d):
        for j in range(d):
            if i == j:
                migration[i, j] = 0
    test = tuple(ms.simulate(
        population_configurations=population_configurations,
        migration_matrix=migration,
        mutation_rate=1 / 4,
        num_replicates=2,
        random_seed=3
        ))
    testsample = MetaSample(test, populations=np.repeat(np.arange(d), 4))
    # a = (testsample.h(average=False, by_population=True))
    # b = testsample.h()

    print(testsample.fst(
        average_sites=True, average_h=False, average_final=True,
        bias=False, replace=True)
        )
    # print(a, b, end='\n')
    # for array in a:
    #     print(array.shape)
    # print(testsample.fst())
    # print(testsample.fst(average=False, return_scalar=True))
    # print(testsample.fst(average=False))


if __name__ == "__main__":
    main()
