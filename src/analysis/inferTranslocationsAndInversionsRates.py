#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import utils.myTools
import utils.myDiags
import utils.myLightGenomes

__doc__ = """
    Starting from a file.sbs, return the number of sbs in a pairwise comparison of chromosomes of two extant genomes.
    !!! (genome1, genome2) should be the same tuple as used before by phylDiag to generate file.sbs !!!
"""
modesOrthos = list(utils.myDiags.FilterType._keys)
arguments = utils.myTools.checkArgs(
    # file.sbs
    [("syntenyBlocks", file),
     ("genome1", file),
     ("genome2", file)],
    [],
    __doc__
)

sbsInPairComp = utils.myDiags.parseSbsFile(arguments['syntenyBlocks'])
genome1 = utils.myLightGenomes.LightGenome(arguments["genome1"])
genome2 = utils.myLightGenomes.LightGenome(arguments["genome2"])
# This loop is used to ensure that every pair of chromosome (c1, c2) has at
# least an empty list of sbs in sbsInPairComp
for c1 in genome1:
    for c2 in genome2:
        # since sbsInPairComp inherit from collections.defaultdict(lambda
        # collections.defaultdict(list)), if a 2dKey did not exist, it is
        # created and the value is initialased to an empty list []
        c1 = str(c1)
        c2 = str(c2)
        if len(sbsInPairComp[c1][c2]) == 0:
            # this is tautologic (already done by defaultdict)
            sbsInPairComp[c1][c2] = []

def nbSyntenyBlocksInMH(c1, c2, sbsInPairComp):
    if c1 not in sbsInPairComp:
        raise ValueError("c2=%s is not a chromosome of genome2" % c2)
    if c2 not in set(c for c in sbsInPairComp[c1]):
        raise ValueError("c2=%s is not a chromosome of genome2" % c2)
    return len(sbsInPairComp[c1][c2])

for c1 in genome1:
    for c2 in genome2:
        c1 = str(c1)
        c2 = str(c2)
        print >> sys.stderr, "%s\t%s\t%s" % (c1, c2, nbSyntenyBlocksInMH(c1, c2, sbsInPairComp))
