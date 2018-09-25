from __future__ import print_function, division
import numpy as np
import msprime as ms


class Sample():
    """
    Coalescent simulation output representing a population class and methods.
    Input argument should be a np.ndarray with 2 dimensions detailed
    in __init__.

    Attributes:
        segsites: the number of segregating sites in the sample
        self.nchrom: the number of chromosomes sampled
        gtmatrix: the original matrix supplied, consisting of [0, 1]
    """

    def __init__(self, input):
        """
        Args:
            input (np.ndarray OR msprime.TreeSequence): A 2 x 2 np.ndarray
                with individual chromosomes
                represented by rows and snp positions represented by columns.
                Alternatively, an object of class "TreeSequence" can be
                provided and the relevant attributes will be inferred.
                Note that this it the transpose of
                msprime.TreeSequence.genotype_matrix()

        """
        if type(input).__name__ == "TreeSequence":
            self.segsites = input.num_sites
            self.nchrom = input.num_samples
            self.type = "TreeSequence"
        else:
            self.segsites = input.shape[1]
            self.nchrom = input.shape[0]
            self.type = "np.ndarray"
        self.input = input

    def gtmatrix(self):
        """
        Sample.gtmatrix() returns a
        self.nchrom X self.sites matrix of sample genotypes.
        """
        if self.type == "TreeSequence":
            return self.input.genotype_matrix().T
        else:
            return self.input

    def h(self, replace=False, average=False, bias=True, **kwargs):
        """
        Calculate heterozygosity over sites in sample.
        Args:
            replace (bool): calculate heterozygosity with replacement or not
            average (bool): whether to average H values for all sites
            bias (bool): whether to apply bias-correction to calculations
        """
        # k=sum of number of mutants per site, klist=list of those sums
        if type == "TreeSequence":
            klist = np.array(
                    [np.count_nonzero(x.genotypes) for x in input.variants()]
                    )
        else:
            klist = np.count_nonzero(self.gtmatrix, axis=0)
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
                    lambda row: hash(tuple(row)), 1, self.gtmatrix
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
                for i in np.arange(self.gtmatrix.shape[0]):
                    if seqid == hash(tuple(self.gtmatrix[i, ])):
                        seqs[seqid]["seq"] = np.array(self.gtmatrix[i, ])
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
            for i in np.arange(self.gtmatrix.shape[0] - 1):
                for j in np.arange(i, self.gtmatrix.shape[0]):
                    k += np.count_nonzero(
                        self.gtmatrix[i, ] != self.gtmatrix[j, ]
                        )
            # Formula: \sum{k_ij} / (nchrom choose 2)
            return np.sum(k) / ((self.nchrom - 1) * self.nchrom / 2)
        elif method is "h":
            return np.sum(self.h(**kwargs))
        else:
            print("Unsupported method")


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
    test3 = ms.simulate(sample_size=10, mutation_rate=1 / 4)
    testsample3 = Sample(test3)
    print(testsample3.gtmatrix)
    print(testsample.pi())


if __name__ == "__main__":
    main()
