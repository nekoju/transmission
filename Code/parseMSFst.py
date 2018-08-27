# parses a multiply parameterized .ms file with equally sized populations and returns
# fst mean and variance to outfile csv with headers
## nchrom is PER POPULATION
# usage: python ./parseMS.py infile.ms nchrom outfile.csv

from __future__ import print_function
import sys
import statistics as stat
import re
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


def main():
    ind = re.compile('[01]+$')
    nchrom = int(sys.argv[2])

    with open(sys.argv[1], "r") as infile, open(sys.argv[3], "w") as out:

        out.write("mean.fst , var.fst\n")
        population = []
        fst = []
        segsites = []
        H = []
        npop = 0
        count = 0
        total = getLines(infile)

        infile.seek(0)
        for line in infile:
            if ind.match(line):
                population.append([int(x) for x in line.strip()])
                i += 1
                if i == nchrom:
                    i = 0
                    npop +=1
                    H.append(getH(population))
                    population = []
            elif line.strip() == "//":
                if H:
                    fst.append(getFst(H, segsites, nchrom, npop))
                i = 0 
                npop = 0
                H = []
            elif line.startswith("ms"):
                count += 1
                print("\r",
                        (float(count) - 1) / float(total) * 100,
                        "% complete",
                        sep = "",
                        end = "")
                if fst:
                    meanVar = ''.join(
                            (str(x) for x in (
                            weightedMean(fst, segsites),
                            " , ",
                            weightedVar(fst, segsites),
                            "\n")))
                    out.write(meanVar)
                    segsites = []
                    H = []
                    fst = []
            elif line.startswith("segsites"):
                segsites.append(int(line.lstrip("segsites: ")))

        fst.append(getFst(H, segsites, nchrom, npop))
        meanVar = "%f , %f\n" % (
                weightedMean(fst, segsites), weightedVar(fst, segsites))
        out.write(meanVar)
        print("\n", end = "")


if __name__ == "__main__":
    main()




