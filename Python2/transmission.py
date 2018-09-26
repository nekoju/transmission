from __future__ import print_function, division
import numpy as np
import msprime as ms
import pdb


class Sample():
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
                msprime.TreeSequence.genotype_matrix()

        """
        if type(popdata).__name__ == "TreeSequence":
            self.segsites = popdata.num_sites
            self.nchrom = popdata.num_samples
            self.type = "TreeSequence"
        else:
            self.segsites = popdata.shape[1]
            self.nchrom = popdata.shape[0]
            self.type = "ndarray"
        self.popdata = popdata

    def __str__(self):
        print(
            "a {nchrom} chromosome sample with "
            "{segsites} segregating sites".format(
                nchrom=self.nchrom, segsites=self.segsites
                )
            )

    def gtmatrix(self):
        """
        Sample.gtmatrix() returns a
        self.nchrom X self.sites matrix of sample genotypes.
        """
        if self.type == "TreeSequence":
            return self.popdata.genotype_matrix().T
        else:
            return self.popdata

    def h(self, replace=False, average=False, bias=True, **kwargs):
        """
        Calculate heterozygosity over sites in sample.
        Args:
            replace (bool): calculate heterozygosity with replacement or not
            average (bool): whether to average H values for all sites
            bias (bool): whether to apply bias-correction to calculations
        """
        # k=sum of number of mutants per site, klist=list of those sums
        klist = self.num_mutants()
        if replace:
            plist = klist / self.nchrom
            hlist = 2 * plist * (1 - plist)
        # heterozygosity should probably be calculated without replacement,
        # especially with small samples
        else:
            hlist = (2 * klist / self.nchrom *
                     (self.nchrom - klist) / (self.nchrom - 1))
        if bias:
            hlist = self.nchrom / (self.nchrom - 1) * hlist
        if average:
            return np.mean(hlist)
        else:
            return hlist

    def num_mutants(self):
        """
        Returns the number of mutant alleles observed at each site.
        Used by other methods.
        """

        if self.type == "TreeSequence":
            num_mutants = np.array(
                [np.count_nonzero(x.genotypes)
                 for x in self.popdata.variants()]
                )
        else:
            num_mutants = np.count_nonzero(self.gtmatrix, axis=0)
        return num_mutants

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


class MetaSample(Sample):
    """
    Class representing a metapopulation sample and associated methods.
    Input argument should be a npop X nchrom X segsites np.ndarray detailed
    in the docstring for __init__, a nchrom X segsites np.ndarray with an
    array listing to which population each chromosome belongs,
    or an msprime.TreeSequence object.
    """

    def __init__(self, popdata, ):
        super(MetaSample, self).__init__(self, popdata)


def main():
                # test = ms.simulate(sample_size=10, mutation_rate=1)
    # matrix = Sample(test.genotype_matrix())
    # print(matrix.gtmatrix)
    # print("h =", matrix.h(average=True))
    # print("pi =", matrix.pi())
    test2 = np.array([[0, 0, 0], [1, 1, 1], [1, 1, 0], [1, 1, 0]])
    testsample = Sample(test2)
    # pi should be 1.125
    print(testsample.pi())
    print(testsample.pi("tajima"))
    print(testsample.pi("h"))
    test3 = ms.simulate(sample_size=10, mutation_rate=1)
    testsample3 = Sample(test3)
    print(testsample3.gtmatrix())
    print(testsample.pi())

    print(testsample3.polymorphic(output="num"))
    print(testsample3.polymorphic(output="which"))


if __name__ == "__main__":
    main()
