from __future__ import print_function, division
import numpy as np


class Pop():
    """
    Coalescent simulation output representing a population class and methods.
    Input argument should be a np.ndarray with 2 dimensions detailed
    in __init__.

    Attributes:
        segsites: the number of segregating sites in the sample
        self.nchrom: the number of chromosomes sampled
        gtmatrix: the original matrix supplied, consisting of [0, 1]
    """

    def __init__(self, matrix):
        """
        Args:
            matrix (np.ndarray): a 2 x 2 np.ndarray with individual chromosomes
                represented by rows and snp positions represented by columns.
                Note that this it the transpose of
                msprime.TreeSequence.genotype_matrix()

        """
        self.gtmatrix = matrix
        self.segsites = self.gtmatrix.shape[1]
        self.nchrom = self.gtmatrix.shape[0]

    def h(self, replace=False, average=False, bias=True):
        """
        Calculate heterozygosity over sites in sample.
        Args:
            replace (bool): calculate heterozygosity with replacement or not
            average (bool): whether to average H values for all sites
            bias (bool): whether to apply bias-correction to calculations
        """
        # k=sum of number of mutants per site, klist=list of those sums
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

    def pi(self, method="nei"):
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
                h: simply the sum of heterozygosities across all sites
        """

        if method is "nei":
            hashes = [hash(self.gtmatrix[i, ]) for i in range(self.nchrom)]
            seqs = dict(set(hashes))
            for seq in seqs.keys:
                seqs[seq]["p"] = np.count_nonzero(seq == hashes) / self.nchrom
                for i in range(self.gtmatrix.shape[0]):
                    if seq == hash(self.gtmatrix[i, ]):
                        seqs[seq]["seq"] = self.gtmatrix[i, ]
                        break
            nucdiv = 0
            for i in seqs:
                for j in seqs:
                    if i != j:
                        nucdiv += (seqs[i]['p'] * seqs[j]['p'] *
                                   np.count_nonzero(
                                       seqs[i][seq] == seqs[j][seq])
                                   )
            return nucdiv
        elif method is "tajima":

            def pairwise(seq, array):
                return np.sum([np.count_nonzero(seq == array[i, ])
                              for i
                              in range(array.shape[0])]
                              )

            k = np.apply_along_axis(lambda row: pairwise(row, self.gtmatrix),
                                    1,
                                    self.gtmatrix)
            return k / ((self.nchrom - 1) * self.nchrom / 2)
        elif method is "h":
            np.sum(self.h())
        else:
            print("Unsupported method")
