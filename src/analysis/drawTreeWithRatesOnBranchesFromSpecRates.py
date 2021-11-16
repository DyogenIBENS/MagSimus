#! /usr/bin/python
# -*- coding: utf-8 -*-

import utils.myTools
import utils.myPhylTree
import libs.myMagSimusTools as myMagSimusTools
from Bio import Phylo
import sys

__doc__ = """
    Draw the species tree with the nb/rates of events contained in specRates
"""

arguments = utils.myTools.checkArgs([('speciesTree.nwk', file), ('specRates', file)], [('type', str, 'all'), ('format', str, 'absolute')], __doc__)
assert arguments['type'] in ['Inv-Trans', 'Fiss-Fus', 'geneEvts', 'all'], arguments['type']
assert arguments['format'] in ['absolute', 'rate'], arguments['absOrRate']
speciesTree = utils.myPhylTree.PhylogeneticTree(arguments['speciesTree.nwk'])
branchSpecificParameters = myMagSimusTools.BranchParameters(speciesTree,
                                                            paramBranchSpeParamFile=arguments['specRates'])
tree = Phylo.read(arguments['speciesTree.nwk'], 'newick')
blabels = {}
for c in tree.find_clades():
    if c == tree.root:
        # do not label the branch leading to the root
        continue
    print >> sys.stderr, 'Branch =', c.name, 'branch_length=', c.branch_length
    assert isinstance(c, Phylo.BaseTree.Clade)
    print branchSpecificParameters.paramVal[c.name]
    nbInv = branchSpecificParameters.paramVal[c.name]['chrInvert']
    nbTrans = branchSpecificParameters.paramVal[c.name]['chrTransloc']
    nbFus = branchSpecificParameters.paramVal[c.name]['chrFusion']
    nbFiss = branchSpecificParameters.paramVal[c.name]['chrFission']
    nbGeneDup = branchSpecificParameters.paramVal[c.name]['geneDup']
    nbGeneLoss = branchSpecificParameters.paramVal[c.name]['geneLoss']
    nbGeneBirth = branchSpecificParameters.paramVal[c.name]['geneBirth']
    print tree.find_elements(c.name)
    if arguments['format'] == 'rate':
        nbInv = "%.2f" % (nbInv / float(c.branch_length))
        nbTrans = "%.2f" % (nbTrans / float(c.branch_length))
        nbFus = "%.2f" % (nbFus / float(c.branch_length))
        nbFiss =  "%.2f" % (nbFiss / float(c.branch_length))
        nbGeneDup =  "%.2f" % (nbGeneDup / float(c.branch_length))
        nbGeneLoss =  "%.2f" % (nbGeneLoss / float(c.branch_length))
        nbGeneBirth =  "%.2f" % (nbGeneBirth / float(c.branch_length))
    propTandemDup = "%2.f" % (branchSpecificParameters.paramBranchSpeParamValue[c.name]['geneTandemDupProp'] * 100)

    if arguments['type'] == 'Inv-Trans':
        blabels[c] = "i=%s\nt=%s" % (nbInv, nbTrans)
    elif arguments['type'] == 'Fiss-Fus':
        blabels[c] = "fiss=%s\nfus=%s" % (nbFiss, nbFus)
    elif arguments['type'] == 'geneEvts':
        blabels[c] = "gDup=%s\ngLoss=%s\ngBirth=%s" % (nbGeneDup, nbGeneLoss, nbGeneBirth)
    else:
        assert arguments['type'] == 'all'
        format1 = [['gDel', '=', '%s' % str(nbGeneLoss), 'inv', '=', '%s' % str(nbInv), 'fis', '=', '%s' % nbFiss],
                   ['gBirth', '=', '%s' % str(nbGeneBirth), 'transl', '=', '%s' % str(nbTrans), 'fus', '=', '%s' % nbFus]]
        blabels[c] = utils.myTools.printTable(format1, sys.stderr, spaceBetweenColumns=1)
        blabels[c] = 'gDup=' + str(nbGeneDup) + '  prop. tandem=' +  str(propTandemDup) + '%\n' + blabels[c]
        # "%s,  %s%%, %s, %s\n" % (nbGeneDup, nbGeneLoss, nbGeneBirth, propTandemDup) +\
        #          "%s, %s, %s, %s" % (nbInv, nbTrans, nbFiss, nbFus)

Phylo.draw(tree, branch_labels=blabels)