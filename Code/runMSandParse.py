# runs MS and parses output a specified number of times using a parameter file
# writes fst mean and variance to outfile csv with headers
## nchrom is PER POPULATION
# usage retrieved from help flag

from __future__ import print_function
import sys
import statistics as stat
import re
import shlex
import argparse
from subprocess import Popen, PIPE
import pdb


def getH(population):
    p = [stat.mean(x) for x in zip(*population)]
    H = [2 * freq * (1 - freq) for freq in p]
    return [p, H]


def weightedMean(values, weights):
    return sum([values[i] * weights[i] for i in xrange(len(values))]) / sum(weights)


def weightedVar(values, weights):
    V1 = sum(weights)
    xBar = weightedMean(values, weights)
    return sum([weights[i] * (values[i] - xBar) ** 2 for i in xrange(len(values))]) / (V1 - 1)


def map_level(f, item, level):
    if level == 0:
        return f(item)
    else:
        return [map_level(f, x, level - 1) for x in item]


def getFst(H, segsites, nchrom, npop):
    HBySite = [zip(*x) for x in zip(*H)]
    HBarPerSite = map_level(stat.mean, HBySite, 2)
    pBar = HBarPerSite[0]
    Ht = [(nchrom * npop) / (nchrom * npop - 1)* 2 * p * (1 - p) for p in pBar]
    Hs = [nchrom / (nchrom - 1) * x for x in HBarPerSite[1]]
    HtBar = stat.mean(Ht)
    HsBar = stat.mean(Hs)
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

    ind = re.compile('[01]+$')
    popSizes = ["%d" % args.nchrom for x in range(args.npop)]
    total = len(params)
    count = 0
    fst = []
    segsites = []

    with open(args.outfile, "w") as out:
        out.write("mean.fst , var.fst\n")
        for row in params:
            tau, rho = row
            const = rho / (tau ** 2 * (3 - 2 * tau) * (2 - rho) + rho)
            reps = {"total": args.nchrom * args.npop, "nsamp": args.nsamp,
                 "theta": args.theta * const, "npop": args.npop, 
                 "mig": args.mig * const, "nchrom": args.nchrom}
            cmd = " ".join(("ms %(total)d %(nsamp)d -t %(theta)f -I %(npop)d" % reps, 
             " ".join(str(args.nchrom) for x in range(args.npop)),
             "%(mig)f" % reps))
            msOut = getOutput(cmd)[1].split("\n")

            population = []
            H = []
            npop = 0
            count += 1

            print("\r%f percent complete" % ((float(count) - 1) /\
                float(total) * 100) , end = "")
            if fst:
                meanVar = "%f , %f\n" % (
                        weightedMean(fst, segsites), weightedVar(fst, segsites))
                out.write(meanVar)
                segsites = []
                fst = []

            for line in msOut:
                if ind.match(line):
                    population.append([int(x) for x in line.strip()])
                    i += 1
                    if i == args.nchrom:
                        i = 0
                        npop +=1
                        H.append(getH(population))
                        population = []
                elif line.strip() == "//":
                    if H:
                        fst.append(getFst(H, segsites, args.nchrom, npop))
                    i = 0 
                    npop = 0
                    H = []
                elif line.startswith("segsites"):
                    segsites.append(int(line.lstrip("segsites: ")))

        fst.append(getFst(H, segsites, args.nchrom, npop))
        meanVar = "%f , %f\n" % (
                weightedMean(fst, segsites), weightedVar(fst, segsites))
        out.write(meanVar)
        print("\r%f percent complete" % (100))
        print("\n", end = "")


if __name__ == "__main__":
    main()




