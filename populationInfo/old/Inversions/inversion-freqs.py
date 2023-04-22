import sys
from collections import defaultdict as d


def sync2freqh(x):
    ''' convert string in SYNC format to dictionary of freqencies where x is a string in sync format'''
    from collections import defaultdict as d
    if x == ".:.:.:.:.:." or x == "0:0:0:0:0:0":
        return ({"A": "na", "T": "na", "C": "na", "G": "na"}, 0)
        return "na", "na"
    nuc = ["A", "T", "C", "G"]
    counts = [int(X) for X in x.split(":")[:4]]
    if sum(counts) == 0:
        return ({"A": 0.0, "T": 0.0, "C": 0.0, "G": 0.0}, 0)
    CO = {X: Y for X, Y in zip(*[nuc, counts])}
    #print(CO, list(zip(*[nuc,counts])))
    h = d(float)
    for k, v in CO.items():
        h[k] = v / float(sum(CO.values()))
    return h, sum(CO.values())


def meanstdv(x):
    ''' calculate mean, stdev and standard error : x must be a list of numbers'''
    from math import sqrt
    if "na" in x:
        x = [X for X in x if X != "na"]
    n, mean, std, se = len(x), 0, 0, 0
    if len(x) == 0:
        return "na", "na", "na"
    for a in x:
        mean = mean + a
    mean = mean / float(n)
    if len(x) > 1:
        for a in x:
            std = std + (a - mean)**2
        std = sqrt(std / float(n - 1))
        se = std / sqrt(n)
    else:
        std = 0
        se = 0
    return mean, std, se


invmarker = open(sys.argv[1], "r")
data = open(sys.argv[2], "r")
names = sys.argv[3].split(",")
invh = d(list)

for l in invmarker:
    if l.startswith("inversion") or l.startswith("#"):
        continue
    a = l.split()
    invh[a[1] + ":" + a[2]] = [a[0], a[3]]

invdata = d(lambda: d(list))
Inv = []
for l in data:
    a = l.rstrip().split()
    pops = a[3:]
    # print(len(pops),len(names))
    I, A = invh[a[0] + ":" + a[1]]
    Inv.append(I)
    for i in range(len(names)):
        # print(pops[i],A,sync2freqh(pops[i]))
        invdata[names[i]][I].append(sync2freqh(pops[i])[0][A])

Inversions = sorted(list(set(Inv)))

print("Sample\t" + "\t".join(Inversions))
for I, v in sorted(invdata.items()):
    AF = []
    for P, V in sorted(v.items()):
        AF.append(str(round(meanstdv(V)[0], 2)))
    print(I + "\t" + "\t".join(AF))
