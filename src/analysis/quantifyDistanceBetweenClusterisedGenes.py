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
    print >> sys.stderr, parent_dir  # prints /<path>/MagSimus/src/
    sys.path.insert(1, parent_dir)
    from libs import myClustering
    del sys, os

import utils.myTools as myTools
import utils.myLightGenomes as myLightGenomes
import utils.myPhylTree as myPhylTree
import utils.myMapping as myMapping

import sys
import pickle
import collections

__doc__ = """ Analyse distances between coevolving genes in modern genomes """

#TODO, allow to start from another ancestor than Amniota for the clustering

arguments = myTools.checkArgs( [("realModernGenome", file),
                                ("ancGenes",file),
                                ("geneevents", file),
                                ("genetimeline", file),
                                ("speciesTree",file)],
                              [("clusteredGenesRatio",float,0.99),
                               ("forceClustering",bool,True),
                               ("in:AllBlocks", str, None),
                               ("out:AllBlocks", str, None),
                               ("minCoevolScore", float, 1)],
                              __doc__)

#Check the coherency of args
assert (arguments["in:AllBlocks"] != None and  arguments["out:AllBlocks"] == None) or arguments["in:AllBlocks"] == None

#Load data from files
modernGenome = myLightGenomes.LightGenome(arguments["realModernGenome"])
ancGenes = myLightGenomes.Families(arguments["ancGenes"])
speciesTree = myPhylTree.PhylogeneticTree(arguments["speciesTree"])
clusteredGenesRatio = arguments['clusteredGenesRatio']
forceClustering = arguments['forceClustering']
clusterator = myClustering.Clusterator(arguments["genetimeline"]) # dans genetimeline les noms des genes sont des strings
file = open(arguments['geneevents'],'r')
geneNamesAmniota = eval(file.readline()) # Lit la premiere ligne de geneEvents : un set(['gn1','gn2',...]), avec les genes de Amniota
file.close()
# Verifie que geneNames est strictement inclus dans byname
assert set(str(x) for x in geneNamesAmniota).issubset(clusterator.byname), set(str(x) for x in geneNamesAmniota).difference(clusterator.byname)

#Change data format
#Change modern names by ancGenes whenever it is possible, other gene names are set to None
#Understand which genes are subject to be clusterized by coevolution in the human genome
newModernGenome = myMapping.labelWithFamNames(modernGenome, ancGenes)
# newModernGenome is the genome rewriten with ancGenes whenever it is possible
# otherwise gene names are set to None.
newModernGenome.computeDictG2Ps()
newModernGenomeByGeneNames = newModernGenome.g2ps
# FIXME: for each family name 'fn', newModernGenomeByGeneNames[fn] is a list of
# tuples: [..., (chrom, index), ...] of the modern genes of this family in newModernGenome
# FIXME: do we need this step ?
# newModernGenome = myMapping.remapFilterGeneContent(newModernGenome, set([None]))

# Use all the branch indices for coevolutionScore calculations.
bIndices = speciesTree.indBranches.values()
indicesDup = bIndices
indicesDel = [val + len(bIndices) for val in bIndices]
# The events vector is equal to [DupG + DelG], and has a total of 8+8=16
# indices.
allindices = frozenset(indicesDup + indicesDel)

#Clusterise genes of ancGenesInModernGenome, or load the result of an ancient clusterisation
if arguments["in:AllBlocks"] != None:
    file = open(arguments["in:AllBlocks"],'r')
    allblocks = pickle.load(file)
    file.close()
else:
    #Gather genes that can be clusterized (ancGenesInModernGenome)
    #1
    ancGenesInModernGenome = set([])
    for chrom in newModernGenome:
        ancGenesInModernGenome.update(set(newModernGenome[chrom]))
    ancGenesInModernGenome.difference_update(set([None]))
    #2
    ancGenesInModernGenome2 = set(newModernGenomeByGeneNames.keys())
    assert ancGenesInModernGenome == ancGenesInModernGenome2

    print >> sys.stderr, "Gene clusterization, this task can take some time ...",
    allblocks = clusterator.mkClusterByScore(ancGenesInModernGenome, allindices,
                                             forceClustering, clusteredGenesRatio, verbose=True)
    # allblocks est du type : [Bloc1, Bloc2, ...]
    # avec Bloc = [(None,geneName1), (coevolvingScoreG1G2,geneName2), ..., (coevolvingScoreGNGN-1,geneNameN)]

    #Write the result of the clusterisation for avoiding to re-compute this step at the next execution of the script
    if arguments["out:AllBlocks"] != None:
        file = open(arguments["out:AllBlocks"],'w')
        pickle.dump(allblocks, file)
        file.close()

#Extract all pairs of clustered genes that have a coevolScore >= minCoevolScore
pairsOfCoevolvingGenes = []
for block in allblocks:
    for (i,(m,g)) in enumerate(block):
        if m >= arguments['minCoevolScore']:
            #Rq, the 1st gene of a block has a coevolScore equal to None, and None is always < arguments['minCoevolScore'] whatever is the value of arguments['minCoevolScore']
            pairsOfCoevolvingGenes.append((block[i-1][1],g))
print >> sys.stderr, "Nb of coevolving genes = %s" % len(pairsOfCoevolvingGenes)

#Gather all clustered genes that are syntenic (on the same chromosome on the modern genome)
nbSyntenicPairsOfCoevolvingGenes = 0
pairsOfCoevolvingGenesByChrom = collections.defaultdict(list)
for (g1,g2) in pairsOfCoevolvingGenes:
    #Only genes that have no paralogs are considered
    if len(newModernGenomeByGeneNames[g1]) == 1 and len(newModernGenomeByGeneNames[g2]) == 1:
        chrom1 = newModernGenomeByGeneNames[g1][0][0]
        chrom2 = newModernGenomeByGeneNames[g2][0][0]
        if chrom1 == chrom2:
            nbSyntenicPairsOfCoevolvingGenes += 1
            pairsOfCoevolvingGenesByChrom[chrom1].append((g1,g2))
print >> sys.stderr, "Nb of syntenic coevolving genes = %s" % sum([len(pairsOfCoevolvingGenesByChrom[chrom]) for chrom in pairsOfCoevolvingGenesByChrom])

#Rq, the same gene can be involved into two pairs of genes
if len(pairsOfCoevolvingGenes) > 0:
    print >> sys.stderr, "Proportion of syntenic pairs of genes = %s" % ( float(nbSyntenicPairsOfCoevolvingGenes) / len(pairsOfCoevolvingGenes))
else:
    print >> sys.stderr, "There are no pairs of coevolving genes with coevolution scores over %s" % arguments['minCoevolScore']

#Compare the proportion of syntenic gene pairs to the expected proportion of syntenic gene pairs from random selection
nbGenes = 0
nbGenesPerChrom = []
for chrom in newModernGenome:
    nbGenesPerChrom.append(len(newModernGenome[chrom]))
    nbGenes += len(newModernGenome[chrom])
probability2GenesOnSameChrom = float(sum(n*n for n in nbGenesPerChrom)) / float(nbGenes*nbGenes)
print >> sys.stderr, "Expected proportion of syntenic pairs of genes by random selection = %s" % probability2GenesOnSameChrom

#Gather normalised distances between genes of each pair of syntenic coevolving genes
normalisedDistances = []
for chrom in pairsOfCoevolvingGenesByChrom:
    for (g1,g2) in pairsOfCoevolvingGenesByChrom[chrom]:
        x1 = newModernGenomeByGeneNames[g1][0][1] #index of g1 on chrom
        x2 = newModernGenomeByGeneNames[g2][0][1] #index of g2 on chrom
        distance = abs(x2 - x1)
        normalisedDistance = float(distance) / (len(newModernGenome[chrom])-1)
        normalisedDistances.append(normalisedDistance)
l = normalisedDistances
averageNormalisedDistance = float(sum(l))/len(l) if len(l) > 0 else None
print >> sys.stderr, "Average normalised distance = %s" % averageNormalisedDistance
#Compare average normalised distance to the expected average normalised distance from random selection
print >> sys.stderr, "Expected average normalised distance by random selection = %s" % str(1.0/3.0)

#Verify that coevolScore are >= minCoevolScore
for (g1,g2) in pairsOfCoevolvingGenes:
    v1 = clusterator.byname[g1]
    v2 = clusterator.byname[g2]
    (coevolScore,n,s) = clusterator.coevoluTionScore(v1, v2, allindices)
    #print >> sys.stderr, "coevolScore1 = %s" % coevolScore
    assert coevolScore >= arguments['minCoevolScore'], "%s >= %s" % (coevolScore, arguments['minCoevolScore'])
#A fortiori these assertions should be true too
countN = collections.defaultdict(int)
countNS = collections.defaultdict(int)
for chrom in pairsOfCoevolvingGenesByChrom:
    for (g1,g2) in pairsOfCoevolvingGenesByChrom[chrom]:
        v1 = clusterator.byname[g1]
        v2 = clusterator.byname[g2]
        (coevolScore,n,s) = clusterator.coevoluTionScore(v1, v2, allindices)
        #print >> sys.stderr, "coevolScore2 = %s" % coevolScore
        countN[n] += 1
        countNS[(n,s)] += 1
        assert coevolScore >= arguments['minCoevolScore'], "%s >= %s" % (coevolScore, arguments['minCoevolScore'])
print >> sys.stderr, "Nb of each value of n for each syntenic pair of coevolving genes =", dict(countN)
print >> sys.stderr, "Nb of each value of (n,s) for each syntenic pair of coevolving genes =", dict(countNS)

#pour une branche entre node et fils
#indices = frozenset(phylTree.indNames[e] for e in phylTree.allNames if phylTree.isChildOf(e, fils) and (fils != e))
#scoreSameEvents.reinit_stats() # Reinitialize la cache du decorateur memoize
#myClustering.scoreSameEvents(v1, v2, lstindices) #FIXME prendre en compte les Nones et les scores >1
#plt.plot(range(1,len(chrLengthsRG)+1), chrLengthsRG, marker='o', linestyle='', color='k', label='Real value', markeredgewidth=1, markersize=3)
#plt.title('Lengths of the %s longest chromosomes' % arguments["limitOfDrawnChrs"])
#plt.xlabel('chromosomes sorted by length')
#plt.ylabel('lengths of chromosomes (in nb of genes)')
#plt.ylim(ymin=0)
#plt.legend()
##plt.show()
#plt.savefig(sys.stdout, format='svg')
