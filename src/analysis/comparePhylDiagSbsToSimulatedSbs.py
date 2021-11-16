#! /usr/bin/python
# -*- coding: utf-8 -*-

import sys
import enum
import itertools
from utils import myTools
from utils import myIntervals
from utils import myLightGenomes
from utils import myDiags
from utils import myPhylTree
from utils import myMapping
from utils import myGenomesDrawer

from libs import myBreakpointsAnalyser

__doc__ = """This script compares the sbs of PhylDiag with the reference sbs recorded during simulation of magSimus
          (with the breakpoint analyser feature)
          """
FilterType = enum.Enum('InFamilies', 'InBothGenomes', 'None')

def analyseMistakeOfPhylDiag(geneNameInvolvedInAMistake, sbsInPairComp, simulatedSbsGenome,
                             genome1, genome2, mistake='Fp'):
    gn = geneNameInvolvedInAMistake
    assert mistake in {'Fp', 'Fn'}
    assert isinstance(simulatedSbsGenome, myLightGenomes.LightGenome)
    assert isinstance(genome1, myLightGenomes.LightGenome)
    assert isinstance(genome2, myLightGenomes.LightGenome)
    assert genome1.withDict is True
    assert genome2.withDict is True
    assert simulatedSbsGenome.withDict is True
    positionsOnGenomes = []
    if mistake == 'Fp':
        for ((c1, c2), sb) in sbsInPairComp.items2d():
            if geneNameInvolvedInAMistake == sb.la[0][0]:
                positionsOnGenomes.append((myLightGenomes.GeneP(c1, sb.l1[0][0]), myLightGenomes.GeneP(c2, sb.l2[0][0])))
            if geneNameInvolvedInAMistake == sb.la[-1][0]:
                positionsOnGenomes.append((myLightGenomes.GeneP(c1, sb.l1[-1][0]), myLightGenomes.GeneP(c2, sb.l2[-1][0])))
    elif mistake == 'Fn':
        posSbsGenome = simulatedSbsGenome.getPosition(gn, default=None)
        assert posSbsGenome is not None
        # le gene se situe a l'extremite d'un sb de reference (chromExtremity mode)
        assert posSbsGenome.idx in {0, len(simulatedSbsGenome[posSbsGenome.c]) - 1}
        positionsOnGenomes.append((genome1.getPosition(gn, default=None), genome2.getPosition(gn, default=None)))
    # print >> sys.stderr, "%s on G1:%s and G2:%s" % (mistake, str(possGenome1), str(possGenome2))
    return positionsOnGenomes

    # distribLenSbsWithFn = collections.defaultdict(int)
    # for (geneN, _) in sFn:
    #     geneP = sbsGenome.getPosition(geneN, default=None)
    #     if geneP:
    #         lenSbFn = len(sbsGenome[geneP.c])
    #         print >> sys.stderr, "Fp: gene %s in chrR %s of len=%s at idx %s" % (geneN, geneP.c, lenSbFn, geneP.idx)
    #         distribLenSbsWithFn[lenSbFn] += 1
    # sumLens = sum(distribLenSbsWithFn.values())
    # percentageLen1 = 100 * float(distribLenSbsWithFn[1]) / sumLens if sumLens != 0 else 'None'
    # print >> sys.stderr, "(%s%% of len1), distribLenSbsWithFn=%s" %\
    #                      (percentageLen1,
    #                       sorted([(l, nb) for (l, nb) in distribLenSbsWithFn.iteritems()], key=lambda x: x[0]))
    # return distribLenSbsWithFn

arguments = myTools.checkArgs(
    [
        # 'data/speciesTree.phylTree'
        ('speciesTree', file),
        # 'Mus.musculus'
        ('species1', str),
        # 'Gallus.gallus'
        ('species2', str),
        # 'Amniota'
        ('ancGenes', str)
    ],
    [
        ('pSimGenomes', str, 'res/simu1/genes.%s.list.bz2'),
        ('pAncGenes', str, 'res/simu1/ancGenes.%s.list.bz2'),
        ('pSimulatedSbs', str, 'res/simu1/sbs.genes.%s.%s.list.bz2'),
        # filterType should be in ['None', 'InFamilies', 'InBothGenomes']
        ('filterType', str, 'InBothGenomes'),
        #('gapMaxs', str, '(0,1,2,3,4,5,10,15,20)'),
        ('gapMaxs', str, '(0,1,20)'),
        #('gapMaxs', str, '(0,1)'),
        #('gapMaxs', str, '(0,)'),
        #('tandemGapMaxs', tuple, (0, 10, 100)),
        # Benchmark this!!!!
        #('gapMaxMicroInv', int, 0),
        #('identifyMonoGenicInvs', bool, False),
        ('tandemGapMax', int, 5),
        ('pThreshold', str, 'None'),
        ('distanceMetric', str, 'CD'),
        ('distinguishMonoGenicDiags', bool, True),
        #('identifyMicroRearrangements', bool, False),

        #
        ('preComputePairwiseSbs', bool, True),
        #
        ('removeSbsOfLen1', bool, False),
        #
        ('reduceToSameGeneNames', bool, False),
        #
        ('oriented', bool, True),

        ('analyseMistakesOneByOne', bool, False),
        ('mistake', str, 'Fn')
    ],
    __doc__,
    # load Only Default Options
    loadOnlyDefaultOptions=False)

# if arguments['gapMaxMicroInv'] == 'None':
#     arguments['gapMaxMicroInv'] = None
# else:
#     arguments['gapMaxMicroInv'] = int(arguments['gapMaxMicroInv'])
if arguments['gapMaxs'] == 'None':
    arguments['gapMaxs'] = None
else:
    arguments['gapMaxs'] = tuple(eval(arguments['gapMaxs']))

# import os
# os.chdir('/home/jlucas/Libs/MagSimus')
# arguments['speciesTree'] = 'data/speciesTree.phylTree'
# arguments['species1'] = 'Mus.musculus'
# arguments['species2'] = 'Gallus.gallus'
# arguments['ancGenes'] = 'Amniota'
# arguments['filterType'] = 'InFamilies'
# arguments['preComputePairwiseSbs'] = False
# arguments['tandemGapMax'] = 0
# arguments['removeSbsOfLen1'] = True
# arguments['reduceToSameGeneNames'] = False

speciesTree = myPhylTree.PhylogeneticTree(arguments['speciesTree'])

if arguments['preComputePairwiseSbs']:
    myBreakpointsAnalyser.computePairwiseSyntenyBlocksOfMagSimus(speciesTree, arguments['pSimGenomes'],
                                                                 arguments['pAncGenes'], arguments['pSimulatedSbs'])

genome1 = myLightGenomes.LightGenome(arguments['pSimGenomes'] % arguments['species1'])
genome2 = myLightGenomes.LightGenome(arguments['pSimGenomes'] % arguments['species2'])
ancGenes = myLightGenomes.Families(arguments['pAncGenes'] % arguments['ancGenes'])
filterType = list(myDiags.FilterType._keys)
filterType = myDiags.FilterType[filterType.index(arguments["filterType"])]
# to have the species combination in alphabetical order
speciesNames = sorted([arguments['species1'], arguments['species2']])
referenceSbs = myLightGenomes.LightGenome(arguments['pSimulatedSbs'] % (speciesNames[0], speciesNames[1]))

if arguments['removeSbsOfLen1']:
    referenceSbs.removeChrsStrictlySmallerThan(2)
nbSbsR = len(referenceSbs.keys())
dictStatsComp = {}
pThreshold = arguments['pThreshold'] if arguments['pThreshold'] != 'None' else None


# variableArg = ('ibwg-imgi-gmmi-tm', [(False, False, 0, None),
#                                      (True, False, 0, None),
#                                      (False, True, 0, None),
#                                      (True, True, 0, None),
#                                      (True, True, 1, None),
#                                     (True, True, 1, 4)])

variableArg = ('ibwg-imgi-gmmi-tm',
[(True, True, 0, None),
(True, True, 1, None)])

#variableArg = ('ibwg-imgi-tm', [(False, False, None)])
# tandemGapMaxs = arguments['tandemGapMaxs']
# identifyMonoGenicInversionS = [False, True]
# pThresholdS = [1.0, 0.1, 0.01, 0.001, None]
# gapMaxMicroInvS = [0, 1, 2, 3, 4, 5]
# distanceMetricS = ['CD', 'MD', 'DPD', 'ED']
progressBar = myTools.ProgressBar(len(arguments['gapMaxs']) * len(variableArg[1]))
progress = 0

# For the analysis of Fp and Fn
# ((g1_tb, mtb2g1, (nCL1, nGL1)), (g2_tb, mtb2g2, (nCL2, nGL2))) =\
#     myDiags.editGenomes(genome1, genome2, ancGenes,
#                         filterType=filterType, labelWith='FamName',
#                         tandemGapMax=arguments['tandemGapMax'], keepOriginal=True)
# g1_tb.computeDictG2P()
# g1_tb.computeDictG2Ps()
# g2_tb.computeDictG2P()
# g2_tb.computeDictG2Ps()
#######

mistake = arguments['mistake']

for varArg in variableArg[1]:
    dictStatsComp[varArg] = {}
    #(ibwg, imgi, tm) = varArg
    (ibwg, imgi, gmmi, tm) = varArg
    for gapMax in arguments['gapMaxs']:
        print >> sys.stderr, "gapMax = %s" % gapMax
        progress += 1
        progressBar.printProgressIn(sys.stderr, progress)
        sbsInPairComp = myDiags.extractSbsInPairCompGenomes(genome1, genome2, ancGenes,
                                                            filterType=filterType,
                                                            gapMax=gapMax,
                                                            distanceMetric=arguments['distanceMetric'],
                                                            tandemGapMax=arguments['tandemGapMax'],
                                                            gapMaxMicroInv=gmmi,  #arguments['gapMaxMicroInv'],
                                                            distinguishMonoGenicDiags=arguments['distinguishMonoGenicDiags'],
                                                            #identifyMonoGenicInvs=varArg,  #arguments['identifyMonoGenicInvs'],
                                                            identifyMonoGenicInvs=imgi,  #arguments['identifyMonoGenicInvs'],
                                                            #identifyMicroRearrangements=arguments['identifyMicroRearrangements'],
                                                            truncationMax=tm,
                                                            identifyMicroRearrangements=ibwg,
                                                            pThreshold=arguments['pThreshold'],
                                                            verbose=True)

        phylDiagSbsGenome = myDiags.buildGenomeFromSbs(sbsInPairComp, sbLengthThreshold=None)

        if arguments['removeSbsOfLen1']:
            phylDiagSbsGenome.removeChrsStrictlySmallerThan(2)
        nbSbs = len(phylDiagSbsGenome.keys())
        if arguments['removeSbsOfLen1']:
            assert all(len(sb) > 1 for sb in phylDiagSbsGenome.values())
            assert all(len(sb) > 1 for sb in referenceSbs.values())

        if arguments['reduceToSameGeneNames']:
            sgSbs = phylDiagSbsGenome.getGeneNames(checkNoDuplicates=False)
            sgReferenceSbs = referenceSbs.getGeneNames()
            # 1st reduce to the same set of genes
            setGenesInBoth = sgSbs & sgReferenceSbs
            setGenesToRemove = (sgSbs | sgReferenceSbs) - setGenesInBoth
            (phylDiagSbsGenome, _, (nbRemovedSbs1, nbRemovedGenes1)) =\
                myMapping.remapFilterGeneContent(phylDiagSbsGenome, setGenesToRemove)
            print >> sys.stderr, (nbRemovedSbs1, nbRemovedGenes1)
            (referenceSbs2, _, (nbRemovedSbs2, nbRemovedGenes2)) =\
                myMapping.remapFilterGeneContent(referenceSbs, setGenesToRemove)
            print >> sys.stderr, (nbRemovedSbs2, nbRemovedGenes2)
        else:
            referenceSbs2 = referenceSbs

        (eff, (sTp, sFp, sFn)) = myIntervals.compareGenomes(phylDiagSbsGenome, referenceSbs2,
                                                             mode='chromExtremity', oriented=arguments['oriented'], verbose=True)
        sTpNs = set(gn for (gn, _) in sTp) if arguments['oriented'] else sTp
        if arguments['analyseMistakesOneByOne']:
            phylDiagSbsGenome.computeDictG2Ps()
            referenceSbs2.computeDictG2P()
            genome1.computeDictG2P()
            genome2.computeDictG2P()
            sF = sFp if mistake == 'Fp' else sFn
            for g in sF:
                if arguments['oriented']:
                    (gn, _) = g
                else:
                    gn = g
                positionsOnGenomes = analyseMistakeOfPhylDiag(gn, sbsInPairComp, referenceSbs2,
                                                              genome1, genome2, mistake=mistake)
                print >> sys.stderr, "mistake-gn in sTp=%s" % str(gn in sTpNs)
                for (pos1, pos2) in positionsOnGenomes:
                    print >> sys.stderr, "pos1=%s, pos2=%s" % (pos1, pos2)
                    CDF1 = "%s:%s-%s" % (pos1.c, max(1, pos1.idx - 20), min(pos1.idx + 19, len(genome1[pos1.c])))
                    CDF2 = "%s:%s-%s" % (pos2.c, max(1, pos2.idx - 20), min(pos2.idx + 19, len(genome2[pos2.c])))

                    (ibwg, imgi, gmmi, tm)
                    myGenomesDrawer.homologyMatrixViewer(genome1, genome2, ancGenes, CDF1, CDF2,
                                                         filterType=filterType,
                                                         gapMax=gapMax,
                                                         distanceMetric=arguments['distanceMetric'],
                                                         tandemGapMax=arguments['tandemGapMax'],
                                                         identifyMonoGenicInvs=imgi,
                                                         identifyMicroRearrangements=ibwg,
                                                         pThreshold=arguments['pThreshold'],
                                                         gapMaxMicroInv=gmmi,
                                                         verbose=False,
                                                         inSbsInPairComp=sbsInPairComp,
                                                         truncationMax=tm,
                                                         # convertGenicToTbCoordinates=True,
                                                         # chromosomesRewrittenInTbs=True,
                                                         switchOnDirectView=True)
                    raw_input()
        dictStatsComp[varArg][gapMax] = (eff, (nbSbs, nbSbsR))

import matplotlib.pyplot as plt
def plot(dictStatsComp, arguments, variableArg):
    fig, axes = plt.subplots(nrows=3)
    gapMaxs = arguments['gapMaxs']
    # fig, axes = plt.subplots(nrows=5)
    # Xs
    # lines = ["-", "--", "-.", ":"]
    lines = ["-"]
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
    linecycler = itertools.cycle(lines)
    colorcycler = itertools.cycle(colors)
    supTitle = "Efficiency analysis %s - %s\n" % (arguments['species1'], arguments['species2']) +\
               'filterType=%s\n' % arguments['filterType'] +\
               'tandemGapMax=%s' % arguments['tandemGapMax']
               #'distanceMetric=%s\n' % arguments['distanceMetric'] +\
               #'identifyMicroRearrangements=%s\n' % arguments['identifyMicroRearrangements'] +\
               #'removeSbsOfLen1=%s\n' % arguments['removeSbsOfLen1'] +\
               #'pThreshold=%s\n' % arguments['pThreshold'] +\
               #'gapMaxMicroInv=%s\n' % arguments['gapMaxMicroInv'] +\
               #'identifyMonoGenicInvs=%s\n' % arguments['identifyMonoGenicInvs'] +\

    fig.suptitle(supTitle, fontsize=12, )
    sensitivityTitle = 'sensitivity'
    specificityTitle = 'specificity'
    subTitles = (sensitivityTitle, specificityTitle)
    minY = {sensitivityTitle: sys.maxint, specificityTitle: sys.maxint}
    for varArg in variableArg[1]:
        label = list(varArg)
        for i, l in enumerate(label):
            if l is True:
                label[i] = '+'
            elif l is False:
                label[i] = '-'
            elif l is None:
                label[i] = '$\\varnothing$'
            else:
                label[i] = str(l)
        label = tuple(label)
        color = next(colorcycler)
        line = next(linecycler)
        # Ys
        statsPerGapMax = dictStatsComp[varArg]
        sensitivitys = [eff.sn for (eff, (nbSbs, nbSbsR )) in statsPerGapMax.values()]
        specificitys = [eff.sp for (eff, (nbSbs, nbSbsR )) in statsPerGapMax.values()]
        # fps = [eff.fp for (eff, (nbSbs, nbSbsTrue )) in statsPerGapMax.values()]
        # fns = [eff.fn for (eff, (nbSbs, nbSbsTrue )) in statsPerGapMax.values()]
        # Yss = (sensitivitys, specificitys, fps, fns)
        Yss = (sensitivitys, specificitys)
        # subTitles = ('sensitivity', 'specificity', 'fps', 'fns')
        for subTitle, ax, Ys in zip(subTitles, axes, Yss):
            # if subTitle in ('sensitivity', 'specificity'):
            #     ax.set_ylim(0, 1)
            # elif subTitle in ('fps', 'fns'):
            #     ax.set_ylim(0, max(Ys))
            minY[subTitle] = min(Ys) if min(Ys) < minY[subTitle] else minY[subTitle]
            ax.plot(xrange(len(gapMaxs)), Ys, color=color, linestyle=line, label='(' + ','.join(label) + ')')
            ax.set_ylabel(subTitle)
            ax.set_ylim((minY[subTitle], 1.0))
            #ax.set_ylim((0.0, 1.0))
            # ax.set_title(subTitle, fontsize=9)
            ax.set_xticks(xrange(len(gapMaxs)))
            ax.set_xticklabels([str(gapMax) for gapMax in gapMaxs])
        nbSbss = [nbSbs for (eff, (nbSbs, nbSbsR )) in statsPerGapMax.values()]
        # diffNbSbs = [(nbSbs - nbSbsTrue ) for (nbSbs, nbSbsTrue ) in zip(nbSbss, nbSbsRs)]
        # axes[-1].plot(gapMaxs, diffNbSbs, color=color, linestyle='-', label=pThreshold)
        axes[-1].plot(xrange(len(gapMaxs)), nbSbss, color=color, linestyle=line, label=label)
        # axes[-1].set_title('nbSbs', fontsize=9)
        axes[-1].set_ylabel('#synteny blocks')
        axes[-1].set_xlabel('gapMax')
        #axes[-1].set_xticks(xrange(len(gapMaxs)), (str(gapMax) for gapMax in gapMaxs))
        #ax.set_xticklabels(xrange(len(gapMaxs)))
        axes[-1].set_xticks(xrange(len(gapMaxs)))
        axes[-1].set_xticklabels([str(gapMax) for gapMax in gapMaxs])
        #axes[-1].set_xlim((0, 20))
    nbSbsEs = [nbSbsR for (eff, (nbSbs, nbSbsR)) in statsPerGapMax.values()]
    #line2, = axes[-1].plot(gapMaxs, nbSbsEs, color='red', linestyle='-', label='Sim', linewidth=2)
    line2, = axes[-1].plot(xrange(len(gapMaxs)), nbSbsEs, color='red', linestyle='-', linewidth=1.5)
    legend2 = axes[-1].legend([line2], ('Simulation',), loc="upper right")
    leg = axes[0].legend(loc=4, ncol=2, shadow=True, title=variableArg[0], fancybox=True)
    # set the linewidth of each legend object
    for legobj in leg.legendHandles:
        legobj.set_linewidth(2.0)
    plt.gca().add_artist(legend2)
    fig.tight_layout()
    plt.show(block=True)

plot(dictStatsComp, arguments, variableArg)

