#!/usr/bin/python
# -*- coding: utf-8 -*-

__doc__ = """
        Input:
            Species Tree
        OutPut:
            Print the pair of extant species and their last common ancestor
"""

import itertools

import utils.myTools
import utils.myPhylTree

arguments = utils.myTools.checkArgs( [("speciesTree",file)], [], __doc__ )

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["speciesTree"])

for a in phylTree.listAncestr:
    for (f1,f2) in itertools.combinations([f for (f,_) in phylTree.items[a]], 2):
        l1 = phylTree.species[f1]
        l2 = phylTree.species[f2]
        for (e1,e2) in itertools.product(l1, l2):
            print "%s\t%s\t%s" % (e1,e2,a)
