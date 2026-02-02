import sys

if len(sys.argv) != 4:
    print("Usage: python script_name.py fileWL fileGenes inputFile")
    sys.exit(1)

fileWL = sys.argv[1]
fileGenes = sys.argv[2]
inputFile = sys.argv[3]

WL = {}
geneID = {}
U = {}

def getTag(tag, line):
    tagOut = line
    if tag not in line:
        return None
    tagOut = line.split(tag, 1)[1].split("\t")[0]
    return tagOut

with open(fileWL) as wl_file:
    for ii, line in enumerate(wl_file, 1):
        WL[line.split()[0]] = ii

with open(fileGenes) as genes_file:
    for ii, line in enumerate(genes_file, 1):
        geneID[line.split()[0]] = ii

nTot = 0

with open(inputFile) as input_file:
    for line in input_file:
        GX = getTag("GX:Z:", line)
        if GX == "0" or GX == "-":
            continue

        CB = getTag("CB:Z:", line)
        if CB == "0" or CB == "-":
            continue

        UB = getTag("UB:Z:", line)
        if UB == "0" or UB == "-":
            continue

        cb = WL.get(CB)
        ge = geneID.get(GX)

        if cb is not None and ge is not None:
            U.setdefault(cb, {}).setdefault(ge, {}).setdefault(UB, 0)
            U[cb][ge][UB] += 1
            nTot += 1

for cb in U:
    for ge in U[cb]:
        print(ge, cb, len(U[cb][ge]))

