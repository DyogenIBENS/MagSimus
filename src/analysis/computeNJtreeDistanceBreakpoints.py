#!/usr/bin/python
# -*- coding: utf-8 -*-
import utils.myDiags as myDiags
import utils.myTools as myTools
import utils.myMaths as myMaths
import utils.myLightGenomes as myLightGenomes
import utils.myPhylTree as myPhylTree
#import os
import sys
import collections
import itertools

__doc__ = """ The user should specify either the speciesTree or the allPComps value that are exclusive"""

# Arguments
#TODO implement a filter depending on filterType
modesOrthos = list(myDiags.FilterType._keys)
arguments = myTools.checkArgs(
    # allPairwiseComparisons and for each comparison the ancGene file for gene
    # family definition.
    [()],
    [('genomes', str, 'genes.%s.list.bz2'),
     ('ancGenes', str, 'ancGenes.%s.list.bz2'),
     # The user may define a species tree to automatically define gene families,
     # using the gene content of the LCA of the two compared species.
     ('speciesTree', str, 'None'),
     # Alternatively, the user may specify his own gene families
     ('allPComps', str, 'None'),
     ("minChromLength", int, 1),
     ("filterType", str, modesOrthos),
     ("tandemGapMax", int, 0),
     ("distanceMetric", str, 'CD'),
     ("gapMax", int, 0),
     ('truncationMax', str, 'None'),
     ('identifyMicroRearrangements', bool, False),
     ("pThreshold", float, 0.001),
     ("out:distanceMatrix", str, 'None'),
     ('out:treeWithDistances', str, "treeWithDistances")
     ],
    __doc__
    )
if arguments['truncationMax'] == 'None':
    arguments['truncationMax'] = None
filterType = myDiags.FilterType[modesOrthos.index(arguments["filterType"])]

genomes = arguments['genomes']
ancGenes = arguments['ancGenes']

assert arguments['speciesTree'] == 'None' and arguments['allPComps'] != 'None' or \
    arguments['speciesTree'] != 'None' and arguments['allPComps'] == 'None'

# Define all the pairwise comparisons
allPComps = []
setOfExtantSpecies = set([])
if arguments['speciesTree'] != 'None':
    phylTree = myPhylTree.PhylogeneticTree(arguments["speciesTree"])
    for anc in phylTree.listAncestr:
        for (f1, f2) in itertools.combinations([f for (f, _) in phylTree.items[anc]], 2):
            l1 = phylTree.species[f1]
            l2 = phylTree.species[f2]
            for (sp1, sp2) in itertools.product(l1, l2):
                setOfExtantSpecies = setOfExtantSpecies | {sp1, sp2}
                allPComps.append((sp1, sp2, anc))
    nbExtantSpecies = len(setOfExtantSpecies)
elif arguments['allPComps'] != 'None':
    # check the integrity of the allPComps file
    f = file(arguments['allPComps'], 'r')
    nbCombinations = 0
    for line in f:
        [sp1, sp2, anc] = line.split()
        setOfExtantSpecies = setOfExtantSpecies | {sp1, sp2}
        nbCombinations += 1
    nbExtantSpecies = len(setOfExtantSpecies)
    nbExpectedCombinations = int(myMaths.combinations(nbExtantSpecies, 2))
    print >> sys.stderr, "Expected number of combinations = %s" % nbExpectedCombinations
    if nbCombinations != nbExpectedCombinations:
        raise ValueError("Not enough comparisons of extant species in %s" % allPComps)
    f.close()
    # Load allPComps file
    f = file(allPComps, 'r')
    for line in f.readlines():
        [sp1, sp2, anc] = line.split()
        allPComps.append((sp1, sp2, anc))
    f.close()

nbSbsPComps = myTools.Dict2d(int)
for (sp1, sp2, anc) in allPComps:
    # paiwrise comparison
    nbSbsSp1Sp2 =\
        len(
            myDiags.extractSbsInPairCompGenomes(
                myLightGenomes.LightGenome(genomes % sp1),
                myLightGenomes.LightGenome(genomes % sp2),
                myLightGenomes.Families(ancGenes % anc),
                minChromLength=arguments['minChromLength'],
                filterType=filterType,
                tandemGapMax=arguments['tandemGapMax'],
                gapMax=arguments['gapMax'],
                distanceMetric=arguments['distanceMetric'],
                truncationMax=arguments['truncationMax'],
                identifyMicroRearrangements=arguments['identifyMicroRearrangements'],
                pThreshold=arguments['pThreshold'],
                verbose=True).intoList()
            )
    print >> sys.stderr, "nb sbs in the pairwise comp (%s, %s) = %s" % (sp1, sp2, nbSbsSp1Sp2)
    nbSbsPComps[sp1][sp2] = nbSbsSp1Sp2
    nbSbsPComps[sp2][sp1] = nbSbsSp1Sp2
    nbSbsPComps[sp1][sp1] = 0
    nbSbsPComps[sp2][sp2] = 0

assert len(nbSbsPComps.keys()) == nbExtantSpecies
assert all([len(nbSbsPComps[sp].keys()) == nbExtantSpecies for sp in nbSbsPComps.keys()])
assert all([nbSbsPComps[sp][sp] == 0 for sp in nbSbsPComps.keys()])

orderedNamesOfExtantSpecies = list(setOfExtantSpecies)
lowTriangularMatrix = []
# Build a low triangular matrix
for (i1, sp1) in enumerate(orderedNamesOfExtantSpecies):
    lowTriangularMatrix.append([])
    for (i2, sp2) in enumerate(orderedNamesOfExtantSpecies):
        if i2 <= i1:
            if i1 == i2:
                assert sp1 == sp2
            lowTriangularMatrix[-1].append(nbSbsPComps[sp1][sp2])
        else:
            break

# Assert that pairwise comparisons are organised in a low triangular manner
for i, _ in enumerate(lowTriangularMatrix):
    #print >> sys.stderr, i+1, len(lowTriangularMatrix[i])
    assert len(lowTriangularMatrix[i]) == i+1

#orderedNamesOfExtantSpecies = list(nbSbsPComps.values()[-1]) + [nbSbsPComps.keys()[-1]]
##print >> sys.stderr, nbSbsPComps
##print >> sys.stderr, orderedNamesOfExtantSpecies
##print >> sys.stderr, nbExtantSpecies
#assert nbExtantSpecies == len(orderedNamesOfExtantSpecies)
#
## Assert that pairwise comparisons are organised in a low triangular manner
#for i, sp1 in enumerate(nbSbsPComps.keys()):
#    print >> sys.stderr, i+1, len(nbSbsPComps[sp1].keys())
#    assert len(nbSbsPComps[sp1].keys()) == i+1

## Matrix of distances between extant species
#lowTriangularMatrix = []
## the first distance is equal to 0.0
#lowTriangularMatrix.append([0.0])
#for sp1 in nbSbsPComps.keys():
#    lowTriangularMatrix.append([])
#    for sp2 in nbSbsPComps[sp1].keys():
#        lowTriangularMatrix[-1].append(nbSbsPComps[sp1][sp2])
#    # Values on the diagonal are equal to 0.0
#    lowTriangularMatrix[-1].append(0.0)

## Now from the synteny blocks, build the matrix of distances into the PHYLIP
## format, cf http://evolution.genetics.washington.edu/phylip/doc/neighbor.html
if arguments['out:distanceMatrix'] != 'None':
    f = open(arguments['out:distanceMatrix'], 'w')
    print >> f, "%6.d" % nbExtantSpecies
    #print >> f, orderedNamesOfExtantSpecies[0][:10], "0.0000"
    for (i1, sp1) in enumerate(orderedNamesOfExtantSpecies):
            print >> f, sp1[:10],
            print >> sys.stderr, sp1[:10],
            for (i2, sp2) in enumerate(orderedNamesOfExtantSpecies):
                if i2 <= i1:
                    if i1 == i2:
                        assert sp1 == sp2
                    print >> f, "%6.0f" % lowTriangularMatrix[i1][i2],
                    print >> sys.stderr, "%6.0f" % lowTriangularMatrix[i1][i2],
                else:
                    break
            # new line
            print >> f
            print >> sys.stderr
    f.close()


# Connection to biopython to reconstruct the tree using the NJ method and
# plot the tree.
import imp
# import Bio.Phylo.TreeConstruction as tc
# TODO fix this
# Even if the PYTHONPATH is uploaded with the local version of Biopython, it is
# always the version of the server: /usr/lib/pymodules/python2.7/Bio that is
# loaded...
tc = imp.load_source('TreeConstruction',
                      '/kingdoms/dyogen/workspace4/workspace4/jlucas/Libs/biopython-1.64/build/lib.linux-x86_64-2.7/Bio/Phylo/TreeConstruction.py')
import Bio.Phylo as Phylo
import matplotlib.pyplot as plt
dm = tc._DistanceMatrix(names=orderedNamesOfExtantSpecies, matrix=lowTriangularMatrix)
constructor = tc.DistanceTreeConstructor()
tree = constructor.nj(dm)
print >> sys.stderr, tree
Phylo.draw_ascii(tree, file=sys.stderr)
Phylo.draw(tree, label_func=lambda c: str(c.branch_length), do_show=False, show_confidence=False)
plt.savefig(arguments['out:treeWithDistances'] + '.svg', format='svg')

# Reproduce the data with the NJ implementation in PHYLIP manually, using the
# distance Matrix.

# TODO wrapper for the 'neighbor wrapper' of PHYLIP
## creating the input file
## cf http://evolution.genetics.washington.edu/phylip/doc/main.html
##   section: An example (Unix, Linux or Mac OS X))
#
## remove previous files
#tmpInputFile = os.system("rm outfile")
#tmpInputFile = os.system("rm outtree")
#
#tmpInputFile = os.path.dirname(arguments['out:distanceMatrix']) + "/" + "tmpInputFile.txt"
#f = open('w', tmpInputFile)
#print >> sys.stderr, arguments['out:distanceMatrix']
## Lower-triangular data matrix
#print >> sys.stderr, "L"
## Answer yes to all other default options
#print >> sys.stderr, "Y"
## Replace output file if it already exists
#print >> sys.stderr, "R"
#f.close()
#
#treeWithDistances = arguments['out:treeWithDistances']
## launch the PHYLIP tool neighbors
#if myTools.which('neighbor') is not None:
#    os.system("neighbor  < %s > %s" % (tmpInputFile, treeWithDistances))
#else:
#    print >> sys.stderr, "Error, 'neighbor' 'is not installed properly. Verify that you compiled it and that the folder of the PHYLIP package containing the binary 'neighbor'' is in your PATH env variable."
#    # This will raise an error
#    os.system("neighbor  < %s > %s" % (tmpInputFile, treeWithDistances))
#f.close()
#
#os.system("rm %s" % (tmpInputFile))
