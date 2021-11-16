#!/usr/bin/python

# MagSimus/src/analysis
# python 2.7.8
# Author: Lucas TITTMANN
# Copyright (c) 2015 IBENS/Dyogen Lucas TITTMANN and Hugues ROEST CROLLIUS
# mail : hrc(at)ens.fr
# Licence: GPL 3.0
# This is free software, you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 (GPL v3) or later and the CeCiLL v2 license in France

# Small script to calculate input for entropy measure assessment

import myDiags

@myTools.tictac
@myTools.verbose
def calcHomologies(g1, g2, ancGenes, output_path, file_prefix, minChromLength=0):
    if isinstance(g1, myGenomes.Genome) and isinstance(g2, myGenomes.Genome):
        g1 = g1.intoDict()
        g2 = g2.intoDict()
    elif isinstance(g1, dict) and isinstance(g2, dict):
        pass
    else:
        raise TypeError('g1 and/or g2 must be either myGenomes.Genome or dict')
    #step 1 :filter genomes and rewrite in tandem blocks if needed
    ##############################################################
    # rewrite genomes by family names (ie ancGene names)
    g1_aID = myMapping.labelWithAncGeneID(g1, ancGenes)
    g2_aID = myMapping.labelWithAncGeneID(g2, ancGenes)
    # genes that are not in ancGene have a aID=None
    # print >> sys.stderr, "genome 2 initially contains %s genes" % sum([len(g2[c2]) for c2 in g2])
    # Must be applied on the two genomes, because of the mode inBothGenomes (InCommonAncestor => not only anchor genes are kept but all genes herited from a gene of the LCA)
    #mfilt2origin1 -> mGf2Go1
    ((g1_aID, mGf2Go1, (nCL1, nGL1)), (g2_aID, mGf2Go2, (nCL2, nGL2))) =\
        filter2D(g1_aID, g2_aID, FilterType[modesOrthos.index('InBothSpecies')], minChromLength)

    def createGeneChrDic(g1_aID):
        from collections import Counter
        g1ong2 = {}
        for k, v in g1_aID.iteritems():
            genes = Counter([gt[0] for gt in v])
            for gene in genes:
                if gene in g1ong2:
                    g1ong2[gene] += [k]*genes[gene]
                else:
                    g1ong2[gene]=[k] * genes[gene]
        return(g1ong2)

    from collections import defaultdict
    g1ong2 = createGeneChrDic(g1_aID)
    g2ong1 = createGeneChrDic(g2_aID)
    a = {}
    for g1chr in g1_aID:
        a[g1chr] = {}
        for g2chr in g2_aID:
            a[g1chr][g2chr] = 0
    # Every homology, wether duplicated or not, counts 1
    matAllDup = copy.deepcopy(a)
    for gene in g1ong2:
        for (c1,c2) in itertools.product(g1ong2[gene], g2ong1[gene]):
            matAllDup[c1][c2] += 1

    # All homologies which are in more than 1 chromosome are ignored
    matNoDup = copy.deepcopy(a)
    for gene in g1ong2:
        c1 = set(g1ong2[gene])
        c2 = set(g2ong1[gene])
        if len(c1) > 1 or len(c2) > 1:
            continue
        matNoDup[c1.pop()][c2.pop()] += 1

    # Homologies are distributed realtive to their duplicate number
    matProb = copy.deepcopy(a)
    for gene in g1ong2:
        nRecHomo = len(g1ong2[gene])*len(g2ong1[gene])
        for (c1, c2) in itertools.product(g1ong2[gene], g2ong1[gene]):
            matProb[c1][c2] += 1.0/nRecHomo

    def writeToCsvFile(csvData, csvFileName, writeRowNames = True):
        text_file = open(csvFileName, 'wb')
        text_file.write('"'+'","'.join([str(i) for i in csvData[csvData.keys()[0]].keys()])+'"\n')
        for recRow in csvData:
            recRowStr = ','.join([str(i) for i in csvData[recRow].values()])+'\n'
            if writeRowNames:
                recRowStr = '"' + str(recRow) + '",' + recRowStr
            text_file.write(recRowStr)
        text_file.close()

    writeToCsvFile(matAllDup, output_path+ file_prefix+'-'+'AllDup.csv')
    writeToCsvFile(matNoDup, output_path+ file_prefix+'-'+'NoDup.csv')
    writeToCsvFile(matProb, output_path+ file_prefix+'-'+'Prob.csv')


    return (matAllDup, matNoDup, matProb)

