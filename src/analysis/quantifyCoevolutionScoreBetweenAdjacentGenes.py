#! /usr/bin/python
# -*- coding: utf-8 -*-
# requires python >2.6

# boilerplate to allow relative imports even if this file is run as main
# Here we want to import ../libs/myClustering
# http://stackoverflow.com/questions/2943847/nightmare-with-relative-imports-how-does-pep-366-work
if __name__ == "__main__" and __package__ is None:
    import sys, os
    # The following assumes the script is in the top level of the package
    # directory.  We use dirname() to help get the parent directory
    # to add to sys.path, so that we can import the current package.  This
    # is necessary since when invoked directly, the 'current' package is
          # not automatically imported.
    parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    print >> sys.stderr, parent_dir # prints /<path>/MagSimus/src/
    sys.path.insert(1, parent_dir)
    from libs import myClustering
    del sys, os

import utils.myTools as myTools
import utils.myLightGenomes as myLightGenomes
import utils.myPhylTree as myPhylTree
import utils.myMapping as myMapping

import sys
from collections import Counter
import random

__doc__ = """ Analyse coevolution scores of adjacent genes in extant genomes """


def printCoevolutionStatsAboutAdjacentGenesIn(genome, clusterator, coevolScoreThreshold=0):
    assert isinstance(genome, myLightGenomes.LightGenome)
    # Use all the branch indices for coevolutionScore calculations.
    bIndices = speciesTree.indBranches.values()
    indicesDup = bIndices
    indicesDel = [val + len(bIndices) for val in bIndices]
    # The events vector is equal to [DupG + DelG], and has a total of 8+8=16
    # indices.
    allindices = frozenset(indicesDup + indicesDel)
    #compute the coevolution scores of adjacent genes
    listOfCoevolutionScores = []
    for chrom in genome:
        for (g1, g2) in myTools.myIterator.slidingTuple(genome[chrom]):
            v1 = clusterator.byname[g1.n]
            v2 = clusterator.byname[g2.n]
            (coevolScore, n, s) = clusterator.coevoluTionScore(v1, v2, allindices)
            listOfCoevolutionScores.append((coevolScore, n, s))
    if coevolScoreThreshold > 0:
        listOfCoevolutionScores = [(coevolScore_, n_, s_) for (coevolScore_, n_, s_) in listOfCoevolutionScores if coevolScore_ > coevolScoreThreshold]
    statisticCoevolScore = dict(Counter([coevolScore_ for (coevolScore_, _, _) in listOfCoevolutionScores]))
    statisticN = dict(Counter([n_ for (_, n_, _) in listOfCoevolutionScores]))
    statisticS = dict(Counter([s_ for (_, _, s_) in listOfCoevolutionScores]))
    statisticCoevolScore = sorted(statisticCoevolScore.iteritems(), key=lambda x: x[0])
    statisticN = sorted(statisticN.iteritems(), key=lambda x: x[0])
    statisticS = sorted(statisticS.iteritems(), key=lambda x: x[0])
    for (name, stats) in [('CoevolScore', statisticCoevolScore), ('n', statisticN), ('s', statisticS)]:
        print >> sys.stderr, name + ' value\t' + '\t'.join(["%4.2f" % value for (value, _) in stats])
        print >> sys.stderr, name + ' count\t' + '\t'.join(["%2.2f" % count for (_, count) in stats])
    return listOfCoevolutionScores

arguments = myTools.checkArgs([("realModernGenome", file),
                               ("ancGenes", file),
                               ("genetimeline", file),
                               ("speciesTree", file)],
                              [],
                              __doc__)

# Load data from files
modernGenome = myLightGenomes.LightGenome(arguments["realModernGenome"])
ancGenes = myLightGenomes.Families(arguments["ancGenes"])
speciesTree = myPhylTree.PhylogeneticTree(arguments["speciesTree"])
# In geneTimeline the gene names are strings
clusterator = myClustering.Clusterator(arguments["genetimeline"])

nbGenes = sum([len(modernGenome[chrom]) for chrom in modernGenome.keys()])
nbChroms = len(modernGenome.keys())
print >> sys.stderr, "The extant genome contains %s genes on %s chromosomes" % (nbGenes, nbChroms)

# Change data format
# Change modern names by ancGenes whenever it is possible, other gene names are set to None

# Remove genes that do not come from genes of the ancGenes
newModernGenome = myMapping.labelWithFamNames(modernGenome, ancGenes)
(newModernGenome, _, _) = myMapping.remapFilterGeneContent(newModernGenome, set([None]))
nbChroms = len(newModernGenome.keys())
nbGenes = sum([len(newModernGenome[chrom]) for chrom in newModernGenome.keys()])
print >> sys.stderr, "The extant genome after filtering 'InAncGenes' contains %s genes on %s chromosomes" % (nbGenes, nbChroms)

listOfCoevolutionScores = printCoevolutionStatsAboutAdjacentGenesIn(newModernGenome, clusterator, coevolScoreThreshold=0)
#printCoevolutionStatsAboutAdjacentGenesIn(newModernGenome, coevolScoreThreshold=0.0001)
print >> sys.stderr, "###################################### Null hypothesis #########################################"

#Now the Null hypothesis:
#Shuffle newModernGenome
# do not do a set since there are duplicates (ancGenes names!)
geneList = [g.n for chrom in newModernGenome for g in newModernGenome[chrom]]
random.shuffle(geneList)
chromLengths = [len(newModernGenome[chrom]) for chrom in newModernGenome]
shuffledGenome = myLightGenomes.LightGenome()
for (i, chromLength) in enumerate(chromLengths):
    for _ in range(chromLength):
        shuffledGenome[i].append(myLightGenomes.OGene(geneList.pop(), _))
assert len(geneList) == 0
listOfCoevolutionScoresNullHypothesis = printCoevolutionStatsAboutAdjacentGenesIn(shuffledGenome, clusterator, coevolScoreThreshold=0)

#plot the histogram
import matplotlib as mpl # this line alow not to use the X server
mpl.use('Agg') # same
from matplotlib import pyplot as plt
plt.figure()
X = [coevolScore for (coevolScore,n,s) in listOfCoevolutionScores if coevolScore > 0.001]
XNullHypothesis = [coevolScore for (coevolScore,n,s) in listOfCoevolutionScoresNullHypothesis if coevolScore > 0.001]
n,bins,patches = plt.hist([X,XNullHypothesis], 10, histtype='bar', label=['coevolScoresReal','coevolScoresNullHypothesis'])
plt.legend(loc='upper center')
plt.xlabel('Coevolution Score Value')
plt.ylabel('Nb of Coevolution Scores')
plt.title('Distribution of Coevolution Scores')
#n,bins,patches = plt.hist(XNullHypothesis,bins, histtype='bar')
plt.savefig(sys.stdout, format='svg')
