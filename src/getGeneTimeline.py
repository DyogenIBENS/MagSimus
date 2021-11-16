#!/usr/bin/python

__doc__ = """
        Returs for each ancestral gene the vectors of evolution (duplication and losses) on each branch
"""

import sys
import collections

import utils.myTools as myTools
import utils.myPhylTree as myPhylTree
import utils.myFile as myFile

arguments = myTools.checkArgs(
    [("speciesTree.phylTree", file),
     ("geneEvents", file)],
    [],
    __doc__
)

speciesTree = myPhylTree.PhylogeneticTree(arguments["speciesTree.phylTree"])

f = myFile.openFile(arguments["geneEvents"], "r")

ininames = eval(f.readline())
print >> sys.stderr, "ini %s contains %d genes" % (speciesTree.root, len(ininames))
alldata = {}


# Recursive function to load informations from geneEvents
def load(father):
    global alldata
    alldata[father] = []
    for (child, bLength) in speciesTree.items.get(father, []):
        t = f.readline().split("\t")
        # trans: {..., oldGeneName1: [newGeneName1], ..., oldGeneName2:
        # [newGeneName2: [newGeneName2a, newGeneName2b]]}
        # oldGeneName1 has been conserved in one copy from 'father' to 'child'
        # oldGeneName2 has been duplicated and has two gene duplicates in 'child'
        trans = eval(t[0])
        # deleted: set([..., geneName, ...])
        deleted = eval(t[1])
        # gained: set([..., geneName, ...])
        gained = eval(t[2])
        alldata[father].append((child, trans, deleted, gained))
        speciesTree.printBranchName(child, sys.stderr)
        # foo = "# %s -> %s branch #" % (father, child)
        # bar = "#" * len(foo)
        # print >> sys.stderr, bar
        # print >> sys.stderr, foo
        # print >> sys.stderr, bar
        nbGenesFather = len(trans.keys()) + len(deleted)
        nbGenesSon = sum([len(hGs) for hGs in trans.values()]) + len(gained)
        nbDuplicates = sum([len(x)-1 for x in trans.values()])
        print >> sys.stderr, "%d genes -> %d genes" % (nbGenesFather, nbGenesSon)
        # FIXME, check why it seems not the good value...
        # print >> sys.stderr, "%d kept genes (no dup included)" % sum([1 for hGs in trans.values() if len(hGs) == 1])
        print >> sys.stderr, "%d gene duplicates (initial gene excluded), [rate=%.2f dup/My]" % (nbDuplicates, float(nbDuplicates)/ float(bLength))
        print >> sys.stderr, "%d new genes (gene births + duplicates of newly born genes) [rate=%.2f birth/My]" % (len(gained), float(len(gained)) / bLength)
        print >> sys.stderr, "%d gene loss (gene in ancestor that are not in the child and have no descendant in the child) [rate=%.2f loss/My]" % (len(deleted), float(len(deleted)) / bLength)
        # print >> sys.stderr, "%d herited (dup included), %d deleted (rate=%.2f), %d appeared (rate=%.2f)" %\
        #                     (sum([len(hGs) for hGs in trans.values()]), len(deleted), float(len(deleted)) / bLength, len(gained), float(len(gained)) / bLength)
        # print >> sys.stderr, "%d Duplications (initial gene excluded), rate = %.2f" % (nbDuplicates, float(nbDuplicates)/ float(bLength))
        load(child)
load(speciesTree.root)
f.close()
#alldata is in the form of {Amniota:[(Theria,trans,deleted,gained),(Sauria,trans,deleted,gained)], Theria:[(Euarchontoglires,trans,deleted,gained),(Monodelphis.domestica,trans,deleted,gained)], ..., Monodelphis.domestica:[]}

# indices of branches
ib = speciesTree.indBranches
# at initiatio, every values of vectors are set to None (meanning that the gene
# has not yet been seen in the branch)
Dvecs = collections.defaultdict(lambda: [None]*len(ib))
Lvecs = collections.defaultdict(lambda: [None]*len(ib))


def newValue(formerValue, nbEvents):
    assert nbEvents is not None
    # nbEvt = 0, it means that the gene is present in the branch
    # nbEvt > 0, means that the gene undergo at least one event
    if formerValue is None or formerValue == 0:
        if nbEvents == 0:
            return 0
        elif nbEvents > 0:
            # at least one gene in the family has an event on this branch
            return 1
        else:
            raise ValueError('nbEvents cannot be negative')
    else:
        assert formerValue == 1
        return 1


# Recursive function to count the nb of events on all branchs
def countEvents(node, g, Lvecg, Dvecg):
    for (child, trans, deleted, gained) in alldata[node]:
        ibc = ib[child]
        if (g in deleted):  #FIXME and ((child not in speciesTree.lstEsp2X) or arguments["count2XSpeciesLoss"]):
            Lvecg[ibc] = newValue(Lvecg[ibc], 1)
            Dvecg[ibc] = newValue(Dvecg[ibc], 0)
        elif g in trans:
            Lvecg[ibc] = newValue(Lvecg[ibc], 0)
            Dvecg[ibc] = newValue(Dvecg[ibc], int(len(trans[g]) > 1))
            for g2 in trans[g]:
                countEvents(child, g2, Lvecg, Dvecg)
                if (g2 != g) and (len(alldata[child]) > 0):  # (len(alldata[child]) > 0) <=> if child is not a leaf
                    countEvents(child, g2, Lvecs[g2], Dvecs[g2])
        else:
            # FIXME one possible reason is an incoherency between gene tree and
            # genes in modern genes.
            raise ValueError("Gene %s is neither in the 'deleted' set nor in the 'trans' set")

# start to count for genes in the initial genome
for g in ininames:
    countEvents(speciesTree.root, g, Lvecs[g], Dvecs[g])
assert set(Dvecs.keys()) == set(Lvecs.keys()),\
    (len(set(Dvecs.keys()).difference(set(Lvecs.keys()))),
     len(set(Lvecs.keys()).difference(set(Dvecs.key()))))
# then for genes that originates in children genomes
for node in alldata:
    for (child, trans, deleted, gained) in alldata[node]:
        if len(alldata[child]) > 0:  # (len(alldata[child]) > 0) <=> if child is not a leaf
            for g in gained:
                assert g not in Lvecs
                assert g not in Dvecs
                countEvents(child, g, Lvecs[g], Dvecs[g])

assert set(Dvecs.keys()) == set(Lvecs.keys()),\
    (len(set(Dvecs.keys()).difference(set(Lvecs.keys()))),
     len(set(Lvecs.keys()).difference(set(Dvecs.key()))))
for g in Dvecs:
    print myFile.myTSV.printLine([g, Dvecs[g], Lvecs[g]])
