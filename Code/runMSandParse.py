# runs MS and parses output a specified number of times using a parameter file
# writes fst mean and variance to outfile csv with headers
## nchrom is PER POPULATION
# usage retrieved from help flag

from __future__ import print_function
import sys
import re
import shlex
import argparse
from subprocess import Popen, PIPE
import pdb
import itertools as it
from multiprocessing import Pool


class Copier(object):
    def __init__(self, arguments):
        self.args = arguments
    def __call__(self, pars):
        getSummary(pars, self.args)


def getH(population):
    p = [weightedMean(x) for x in zip(*population)]
    H = [2 * freq * (1 - freq) for freq in p]
    return [p, H]


def weightedMean(values, weights = None):
    if weights:
        return sum([values[i] * weights[i] 
            for i in xrange(len(values))]) / sum(weights)
    else:
        return sum(values) / len(values)


def weightedVar(values, weights):
    V1 = float(sum(weights))
    xBar = weightedMean(values, weights)
    try:
        return sum([weights[i] * (values[i] - xBar) ** 2 
            for i in xrange(len(values))]) / (V1 - 1)
    except:
        return 0


def map_level(f, item, level):
    if level == 0:
        return f(item)
    else:
        return [map_level(f, x, level - 1) for x in item]


def getFst(H, segsites, nchrom, npop):
    HBySite = [zip(*x) for x in zip(*H)]
    HBarPerSite = map_level(weightedMean, HBySite, 2)
    pBar = HBarPerSite[0]
    Ht = [(nchrom * npop) / (nchrom * npop - 1)* 2 * p * (1 - p) for p in pBar]
    Hs = [nchrom / (nchrom - 1) * x for x in HBarPerSite[1]]
    HtBar = weightedMean(Ht)
    HsBar = weightedMean(Hs)
    return 1 - HsBar / HtBar


def getLines(infile):
    total = 0
    for line in infile:
        if line.startswith("ms"):
            total += 1
    return total


def getOutput(cmd):
    args = shlex.split(cmd)
    proc = Popen(args, stdout = PIPE, stderr = PIPE)
    out, err = proc.communicate()
    exitcode = proc.returncode
    return exitcode, out, err


def getPi(populations):
    A = tuple([tuple(x) for x in it.chain(*populations)])
    aHash = tuple([hash(x) for x in A])
    segsites = range(len(A[0]))
    terms = []
    for i in range(len(A)):
        p_i = float(len([1.0 for x in aHash if x == aHash[i]])) / float(len(A))
        for j in range(i + 1, len(populations)):
            p_j = float(len([1.0 for x in aHash if x == aHash[j]])) \
                    / float(len(A))
            ids = tuple(it.chain(*[zip(segsites, A[x]) for x in (i, j)]))
            pi_ij = float(len(set(ids)) - len(segsites)) / float(len(segsites))
            terms.append(p_i * p_j * pi_ij)
    return 2 * sum(terms)


def getSummary(row, args):
    H = []
    pi = []
    fst = []
    segsites = []
    population = []
    npop = 0
    tau, rho = row
    const = rho / (tau ** 2 * (3 - 2 * tau) * (2 - rho) + rho)
    reps = {"total": args.nchrom * args.npop, "nsamp": args.nsamp,
         "theta": args.theta * const, "npop": args.npop, 
         "mig": args.mig * const, "nchrom": args.nchrom}
    cmd = " ".join(("ms %(total)d %(nsamp)d -t %(theta)f -I %(npop)d" % reps, 
     " ".join(str(args.nchrom) for x in range(args.npop)),
     "%(mig)f" % reps))
    msOut = getOutput(cmd)[1].split("\n")

    for line in msOut:
        if ind.match(line):
            population.append([float(x) for x in line.strip()])
            i += 1
            if i == args.nchrom:
                i = 0
                npop +=1
                H.append(getH(population))
                populations.append(population)
                population = []
        elif line.strip() == "//":
            if H:
                fst.append(getFst(H, segsites, args.nchrom, npop))
                pi.append(getPi(populations))
            i = 0 
            npop = 0
            H = []
            populations = []
        elif line.startswith("segsites"):
            segsites.append(float(line.lstrip("segsites: ")))
    fst.append(getFst(H, segsites, args.nchrom, npop))
    pi.append(getPi(populations))
    meanVar = "%f , %f, %f\n" % (
            weightedMean(fst, segsites),
            weightedVar(fst, segsites),
            weightedMean(pi),
            tau,
            rho)
    return meanVar
           
            
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("params", 
            help = "The path of a .csv file containing float values of tau and rho,\
                    one line per MS iteration. Headers will be stripped")
    parser.add_argument("-n",
            "--nchrom",
            type = int,
            help = "Number of chromosomes sampled from each population")
    parser.add_argument("-p",
            "--npop",
            type = int,
            help = "Number of populations per sample")
    parser.add_argument("-s",
            "--nsamp",
            type = int,
            help = "Number of samples per iteration")
    parser.add_argument("-t",
            "--theta",
            type = float,
            help = "Value of theta (2*Ne*u) for the host population")
    parser.add_argument("-M",
            "--mig",
            type = float,
            help = "Population migration parameter (2*Ne*m) for each simulation"
            )
    parser.add_argument("-c",
            "--ncore", 
            nargs = "?",
            default = "1",
            type = int,
            help = "Number of processing cores")
    parser.add_argument("outfile",
            help = "Output file")
    args = parser.parse_args()

    params = []
    with open(args.params) as pars:
        for line in pars:
            try:
                params.append([float(x) for x in 
                    filter(None, line.strip().split(","))])
            except:
                pass

    popSizes = ["%d" % args.nchrom for x in range(args.npop)]

    pool = Pool(args.ncore)
    pdb.set_trace()
    summaries = pool.map(Copier(args), params)

    with open(args.outfile, "w") as out:
        out.write("mean.fst, var.fst, pi, tau, rho\n")
        out.writelines(", ".join(str(i) for i in j) + "\n" for j in summaries)


if __name__ == "__main__":
    ind = re.compile('[01]+$')
    main()




