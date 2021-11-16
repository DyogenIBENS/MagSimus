#!/usr/bin/python
# -*- coding: utf-8 -*-
import itertools
from numpy import random

import utils.myTools as myTools
import utils.myLightGenomes as myLightGenomes
import utils.myMapping as myMapping
import utils.myPhylTree as myPhylTree
import sys
import collections
import matplotlib.pyplot as plt

GeneDupStat = collections.namedtuple('GeneDupStat', ['d', 'dd', 'td'])

def computeProbaTandemDupOrientationComparedToIni(geneOrientationsInTbs):
    tnGTD = sum([len(tb) - 1 for tb in geneOrientationsInTbs])
    listOfTbSpecificProbas = []
    for geneOrientations in geneOrientationsInTbs:
        counter = collections.Counter(geneOrientations)
        nbO = {}
        nbG = len(geneOrientations)
        # Probability to consider a tandem dup in this tb = nb of tandem dup in this tb / total number of tandem dup
        pTandemDupInThisTb = float(nbG - 1) / tnGTD
        nbO[+1] = counter[+1]
        nbO[-1] = counter[-1]
        assert nbO[+1] + nbO[-1] == nbG
        pIniG = {}
        # Hypothesis: all genes of the tb have the same probability to be the initial gene (before all tandem dups)
        # Thus the probability, in this tb, that the initial gene had a +1 orientation is:
        pIniG[+1] = float(nbO[+1]) / nbG
        # And the probability, in this tb, that the initial gene had a -1 orientation is:
        pIniG[-1] = float(nbO[-1]) / nbG
        assert 0.999999999 <= pIniG[+1] + pIniG[-1] <= 1.0000001

        pTandemDupSame = {}
        # In this tb, probability that the initial gene had an orientation +1 and a tandem dup has the same orientation +1
        pTandemDupSame[+1] = pIniG[+1] * float((nbO[+1] - 1)) / (nbG - 1)
        # In this tb, probability that the initial gene had an orientation -1 and a tandem dup has the same orientation -1
        pTandemDupSame[-1] = pIniG[-1] * float((nbO[-1] - 1)) / (nbG - 1)
        # In this tb, probability that a tandem dup has the same orientation as the orientation of the ini gene
        pTandemDupInTbSameOAsIni = pTandemDupSame[+1] + pTandemDupSame[-1]

        pTandemDupDiff = {}
        # In this tb, probability that the initial gene had an orientation +1 and a tandem dup has the diff orientation -1
        pTandemDupDiff[+1] = pIniG[+1] * float(nbO[-1]) / (nbG - 1)
        # In this tb, probability that the initial gene had an orientation -1 and a tandem dup has the diff orientation +1
        pTandemDupDiff[-1] = pIniG[-1] * float(nbO[+1]) / (nbG - 1)
        # In this tb, probability that a tandem dup has the same orientation as the orientation of the ini gene
        pTandemDupInTbDiffOComparedToIni = pTandemDupDiff[+1] + pTandemDupDiff[-1]

        listOfTbSpecificProbas.append((pTandemDupInThisTb, pTandemDupInTbSameOAsIni, pTandemDupInTbDiffOComparedToIni))

    assert 0.999999999 <= sum(pTandemDupInThisTb for (pTandemDupInThisTb, _, _) in listOfTbSpecificProbas) <= 1.0000001
    # probability that a tandem dup, in the whole genome, has the same orientation of the initial gene of its tb
    pTandemDupSameOAsIni = sum(pTandemDupInThisTb * pTandemDupInTbSameOAsIni for (pTandemDupInThisTb, pTandemDupInTbSameOAsIni, _) in listOfTbSpecificProbas)
    pTandemDupDiffOComparedToIni = sum(pTandemDupInThisTb * pTandemDupInTbDiffOComparedToIni for (pTandemDupInThisTb, _, pTandemDupInTbDiffOComparedToIni) in listOfTbSpecificProbas)
    return (pTandemDupSameOAsIni, pTandemDupDiffOComparedToIni)

def fsimuGeneOrientationsInTbs(realGeneOrientationsInTbs, probaSameOrientation=0.75, duplicatedGeneIsTheIni=False):
    simuGeneOrientationsInTbs = []
    for nbGenesInTbs in [len(tb) for tb in realGeneOrientationsInTbs]:
        currTbOs = []
        #first orientation
        s = random.choice([-1, +1])
        currTbOs.append(s)
        for _ in range(nbGenesInTbs - 1):
            if duplicatedGeneIsTheIni:
                sOfTheDuplicatedGene = currTbOs[0]
            else:
                sOfTheDuplicatedGene = random.choice(currTbOs)
            if random.random() <= probaSameOrientation:
                currTbOs.append(sOfTheDuplicatedGene)
            else:
                currTbOs.append(-sOfTheDuplicatedGene)
        simuGeneOrientationsInTbs.append(currTbOs)
    return simuGeneOrientationsInTbs

# Arguments
#TODO implement a filter depending on filterType
arguments = myTools.checkArgs(
                [('root', str), ('speciesTree', str), ("genomeFiles", str), ("familiesFiles", str)],
                # [('root', str, 'Amniota'),
                #  ('speciesTree', str, 'data/speciesTree.phylTree'),
                #  ("genomeFiles", str, 'data/genesST.%s.list.bz2'),
                #  ("familiesFiles", str, 'data/ancGenes.%s.list.bz2')],
                [],
                __doc__#,
                #loadOnlyDefaultOptions=True,
                )
# import  os
# os.chdir('/home/jlucas/Libs/MagSimus')

speciesTree = myPhylTree.PhylogeneticTree(arguments['speciesTree'])
root = arguments['root']
genomes = {}
for extantSp in speciesTree.species[root]:
    genomes[extantSp] = myLightGenomes.LightGenome(arguments['genomeFiles'] % extantSp)
families = myLightGenomes.Families(arguments['familiesFiles'] % root)

dupsStatsIn = collections.defaultdict(dict)
highestTandemGapMax = 15
tandemGapMaxs = range(highestTandemGapMax + 1)
#tandemGapMaxs = [15]
for (extantSp, genome) in genomes.iteritems():
    genome_fID = myMapping.labelWithFamID(genome, families)
    filtGenome = True
    if filtGenome:
        (genomeFilt, gf2gfID, _) = myMapping.remapFilterGeneContent(genome_fID, {None})
        pass
    else:
        genomeFilt = genome_fID
        gf2gfID = None

    familySizes = collections.Counter(genomeFilt.getGeneNames(asA=list, checkNoDuplicates=False))
    tnGD = sum(nbGn - 1 for nbGn in familySizes.values())
    print >> sys.stderr, "Number of duplicates (original gene exluded) = %s" % tnGD

    for tandemGapMax in tandemGapMaxs:
        # (genome_tb, tb2gfID, nGTD) = myMapping.remapRewriteInTb(genomeFilt, tandemGapMax=0, mOld=gf2gfID)
        # print >> sys.stderr, "Number of tandem duplicates (original gene excluded) = %s (%.2f%%) \t with tandem_gap = %s" % (nGTD, 100 * float(nGTD)/float(nGD), 0)
        (genome_tb, tb2gfID, tnGTD) = myMapping.remapRewriteInTb(genomeFilt, tandemGapMax=tandemGapMax, mOld=gf2gfID)
        print >> sys.stderr, "Number of tandem duplicates (original gene excluded) = %s (%.2f%%) \t with tandem_gap = %s" % (tnGTD, 100 * float(tnGTD)/float(tnGD), tandemGapMax)

        histTbSizes = collections.defaultdict(int)
        geneOrientationsInTbs = []
        nGDD = collections.defaultdict(lambda: -1)
        nGD2 = 0
        tnGTD2 = 0
        for (chr, chrom) in genome_tb.iteritems():
            for (idxTb, tb) in enumerate(chrom):
                nGDD[tb.n] += 1
                sizeTb = len(tb2gfID[chr][idxTb])
                idxGs = tb2gfID[chr][idxTb]
                if len(idxGs) > 1:
                    tnGTD2 += sizeTb - 1
                    geneOrientationsInTbs.append([])
                    for idxG in idxGs:
                        geneOrientationsInTbs[-1].append(genome_fID[chr][idxG].s)
        tnGDD = sum(n for n in nGDD.values())

        assert tnGTD2 == tnGTD
        assert tnGD == tnGDD + tnGTD

        dupsStatsIn[extantSp][tandemGapMax] = (GeneDupStat(tnGD, tnGDD, tnGTD), geneOrientationsInTbs)

for extantSpecies in ['Homo.sapiens', 'Mus.musculus', 'Monodelphis.domestica', 'Gallus.gallus', 'Canis.lupus.familiaris']:
    #print >> sys.stderr, extantSpecies
    # extantSpecies = 'Homo.sapiens'
    # extantSpecies = 'Mus musculus'
    # extantSpecies = 'Monodelphis domestica'
    # extantSpecies = 'Gallus gallus' # plus eleve que 75%
    tandemGapMax = highestTandemGapMax
    geneOrientationsInTb = dupsStatsIn[extantSpecies][tandemGapMax][1]
    (pTandemDupSameOAsIni, pTandemDupDiffOComparedToIni) = computeProbaTandemDupOrientationComparedToIni(geneOrientationsInTb)
    print >> sys.stderr, "%s -> %s tandemGapMax=%s:" % (root, extantSpecies, tandemGapMax)
    print >> sys.stderr, "proba tandem dup same o as the ini gene = %s" % pTandemDupSameOAsIni
    print >> sys.stderr, "proba tandem dup opposite o compared to the ini gene = %s" % pTandemDupDiffOComparedToIni
    # Try to reproduce this statistic with simulations
    probaSameOrientation = 0.75
    for probaSameOrientation in[0.5, 0.60, 0.65, 0.75, 0.80, 0.85]:
        print >> sys.stderr, "\n"
        print >> sys.stderr, probaSameOrientation
        print >> sys.stderr, "\n"
        simuGeneOrientationsInTbs = fsimuGeneOrientationsInTbs(geneOrientationsInTbs, probaSameOrientation=probaSameOrientation, duplicatedGeneIsTheIni=False)
        (pTandemDupSameOAsIni, pTandemDupDiffOComparedToIni) = computeProbaTandemDupOrientationComparedToIni(simuGeneOrientationsInTbs)
        print >> sys.stderr, "Null hypothesis: tandem dup have %s chance to be in the same o as the duplicatedGene gene" % probaSameOrientation
        print >> sys.stderr, "proba tandem dup same o as the ini gene = %s" % pTandemDupSameOAsIni
        print >> sys.stderr, "proba tandem dup opposite o compared to the ini gene = %s" % pTandemDupDiffOComparedToIni

colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
colorcycler = itertools.cycle(colors)

fig = plt.figure()
for extantSp in dupsStatsIn.keys():
    Xs = []
    Ys = []
    for (tandemGapMax, (ds, _)) in dupsStatsIn[extantSp].iteritems():
        assert isinstance(ds, GeneDupStat)
        Xs.append(tandemGapMax)
        Ys.append(float(ds.td) / ds.d)
    plt.plot(Xs, Ys, label=extantSp)
plt.ylim((0.0, 1.0))
plt.xlabel('tandemGapMax in ancGenes')
plt.ylabel('% (tandem duplications/duplications)')
plt.legend(loc=4, fontsize=10)
plt.show()

fig = plt.figure()
for extantSp in dupsStatsIn.keys():
    Xs = []
    Ys = []
    for (tandemGapMax, (ds, geneOrientationsInTbs)) in dupsStatsIn[extantSp].iteritems():
        (pTandemDupSameOAsIni, pTandemDupDiffOComparedToIni) = computeProbaTandemDupOrientationComparedToIni(geneOrientationsInTbs)
        Xs.append(tandemGapMax)
        Ys.append(pTandemDupSameOAsIni)
    plt.plot(Xs, Ys, label=extantSp, color=colorcycler.next())
plt.ylim(0, 1)
plt.xlabel('tandemGapMax in ancGenes')
plt.ylabel('proba. tandem dup. has same o. has the ini. anc. gene')
plt.legend(loc=3)
plt.show()

fig = plt.figure()
for extantSp in dupsStatsIn.keys():
    if extantSp != 'Homo sapiens':
        continue
    Xs = []
    YsAll = []
    YsSize2 = []
    YsSizeOver2 = []
    for (tandemGapMax, (_, geneOrientationsInTbs)) in dupsStatsIn[extantSp].iteritems():
        pTandemDupO = {}
        pTandemDupO['all'] = computeProbaTandemDupOrientationComparedToIni(geneOrientationsInTbs)[0]
        # tbs of size 2
        geneOrientationsInTbsOfSize2 = [geneOrientationsInTb for geneOrientationsInTb in geneOrientationsInTbs if len(geneOrientationsInTb) <= 2]
        assert all([len(geneOrientationsInTb) == 2 for geneOrientationsInTb in geneOrientationsInTbsOfSize2])
        pTandemDupO['size2'] = computeProbaTandemDupOrientationComparedToIni(geneOrientationsInTbsOfSize2)[0]
        # tbs of size over 2
        geneOrientationsInTbsOfSizeOver2 = [geneOrientationsInTb for geneOrientationsInTb in geneOrientationsInTbs if len(geneOrientationsInTb) > 2]
        pTandemDupO['sizeOver2'] = computeProbaTandemDupOrientationComparedToIni(geneOrientationsInTbsOfSizeOver2)[0]
        Xs.append(tandemGapMax)
        YsAll.append(pTandemDupO['all'])
        YsSize2.append(pTandemDupO['size2'])
        YsSizeOver2.append(pTandemDupO['sizeOver2'])
    plt.plot(Xs, YsAll, label=extantSp + ' all', color=colorcycler.next())
    plt.plot(Xs, YsSize2, label=extantSp + ' size2', color=colorcycler.next())
    plt.plot(Xs, YsSizeOver2, label=extantSp + ' sizeOver2', color=colorcycler.next())
plt.ylim(0, 1)
plt.xlabel('tandemGapMax in ancGenes')
plt.ylabel('proba. tandem dup. has same o. has the ini. anc. gene')
plt.legend(loc=4)
plt.show()

fig = plt.figure()
tandemGapMax = highestTandemGapMax
barWidth = (1.0 - 0.2) / len(dupsStatsIn.keys())
for idxSp, extantSp in enumerate(dupsStatsIn.keys()):
    (ds, geneOrientationsInTbs) = dupsStatsIn[extantSp][tandemGapMax]
    histTbSizes = collections.Counter([len(tb) for tb in geneOrientationsInTbs])
    Xs = [x + 0.1 + idxSp * barWidth for x in histTbSizes.keys()]
    assert 0.99999 <= float(sum([len(tb) - 1 for tb in geneOrientationsInTbs])) / ds.td <= 1.00001, "%s/%s" %\
                (float(sum([len(tb) - 1 for tb in geneOrientationsInTbs])), ds.td)
    Ys = [float(sizeTb -1)/ds.td * 100 for sizeTb in histTbSizes.values()]
    plt.bar(Xs, Ys, width=barWidth, label=extantSp, color=colorcycler.next())
plt.title("histogram of tandem blocks for tandemGapMax = %s" % tandemGapMax)
xlim = 10
plt.xlim(2, xlim)
plt.xticks([ind + 0.5 for ind in range(2, xlim)], [str(ind) for ind in range(2, xlim)])
plt.xlabel('size of tandem blocks')
plt.ylabel('% of tandem duplicates')
plt.legend(loc=1)
plt.show()



