#!/usr/bin/python
# -*- coding: utf-8 -*-
import utils.myTools as myTools
import utils.myLightGenomes as myLightGenomes
import utils.myPhylTree as myPhylTree
import utils.myMaths as myMaths
import sys

arguments = myTools.checkArgs(
                [],
                [('speciesTree', str, 'data/speciesTree.phylTree'),
                 ("genomeFiles", str, 'data/genesST.%s.list.bz2')],
                __doc__,
                loadOnlyDefaultOptions=True,
                )

# LOAD DATA
# the species tree specifies the extant (leaf) genomes that are considered for computing statistics
speciesTree = myPhylTree.PhylogeneticTree(arguments['speciesTree'])
extantGenomes = {}
for extantSpecies in speciesTree.listSpecies:
    extantGenomes[extantSpecies] = myLightGenomes.LightGenome(arguments['genomeFiles'] % str(extantSpecies))

# COMPUTE RES
listOfWas = []
for (genomeName, genome) in extantGenomes.iteritems():
    wa = myMaths.myStats.getWeightedAverage([len(chrom) for chrom in genome.values()])
    print >> sys.stderr, "weighted average of chrom lengths in " + genomeName + " = %.2f genes" % wa
    listOfWas.append(wa)
meanWa = float(sum(listOfWas)) / len(listOfWas)
print >> sys.stderr, "the mean weighted average is = %.2f" % meanWa


