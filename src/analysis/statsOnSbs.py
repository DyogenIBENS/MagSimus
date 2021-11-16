#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
import collections
import pickle
from utils import myTools, myDiags, myLightGenomes, myMaths, myCondor
import os
os.chdir('/home/jlucas/Libs/MagSimus')
if False:
    arguments = myTools.checkArgs(
                    [('genome1', file), ('genome2', file), ('ancGenes', file)],
                    [('tandemGapMaxs', str, '(0,1,2,3,4,5,6,7,8,9,10,15,20)'),
                     ('gapMaxs', str, '(0,1,2,3,4,5,6,7,8,9,10,15,20)'),
                     # ('tandemGapMaxs', str, '(0,1)'),
                     # ('gapMaxs', str, '(0,1)'),
                     ('filterType', str, 'InBothGenomes'),
                     ('gapMaxMicroInv', int, 1),
                     #('identifyMonoGenicInvs', bool, True),
                     ('identifyMonoGenicInvs', bool, False),
                     ('pThreshold', float, 1.0),
                     ('distanceMetric', str, 'CD'),
                     #('distinguishMonoGenicDiags', bool, True),
                     ('distinguishMonoGenicDiags', bool, False),
                     ('identifyMicroRearrangements', bool, True),
                     ('truncationMax', int, None),
                     ],
                    __doc__,
                    loadOnlyDefaultOptions=True,
                    )
    import os
    os.chdir('/home/jlucas/Libs/MagSimus')
    # arguments['genome1'] = 'data/genesST.Mus.musculus.list.bz2'
    # arguments['genome2'] = 'data/genesST.Gallus.gallus.list.bz2'
    # arguments['ancGenes'] = 'data/ancGenes.Amniota.list.bz2'
    # arguments['genome1'] = 'data/genesST.Homo.sapiens.list.bz2'
    # arguments['genome2'] = 'data/genesST.Mus.musculus.list.bz2'
    # arguments['ancGenes'] = 'data/ancGenes.Euarchontoglires.list.bz2'
    # TODO
    arguments['genome1'] = 'res/simu1/genes/genes.Homo.sapiens.list.bz2'
    arguments['genome2'] = 'res/simu1/genes/genes.Mus.musculus.list.bz2'
    arguments['ancGenes'] = 'res/simu1/ancGenes/ancGenes.Euarchontoglires.list.bz2'
    tandemGapMaxs = tuple(eval(arguments['tandemGapMaxs']))
    gapMaxs = tuple(eval(arguments['gapMaxs']))
    filterType = list(myDiags.FilterType._keys)
    filterType = myDiags.FilterType[filterType.index(arguments["filterType"])]
    wa_for = {}
    nbSbs_for = {}
    for tgm in tandemGapMaxs:
        for gm in gapMaxs:
            sbsInPairComp = myDiags.extractSbsInPairCompGenomes(myLightGenomes.LightGenome(arguments['genome1']),
                                                                myLightGenomes.LightGenome(arguments['genome2']),
                                                                myLightGenomes.Families(arguments['ancGenes']),
                                                                filterType=filterType,
                                                                gapMax=gm,
                                                                distanceMetric=arguments['distanceMetric'],
                                                                tandemGapMax=tgm,
                                                                gapMaxMicroInv=arguments['gapMaxMicroInv'],
                                                                distinguishMonoGenicDiags=arguments['distinguishMonoGenicDiags'],
                                                                identifyMonoGenicInvs=arguments['identifyMonoGenicInvs'],
                                                                identifyMicroRearrangements=arguments['identifyMicroRearrangements'],
                                                                truncationMax=arguments['truncationMax'],
                                                                pThreshold=arguments['pThreshold'],
                                                                verbose=True)
            nbOfSbs = sum([len(sbsInPairComp[c1][c2]) for (c1, c2) in sbsInPairComp.keys2d()])
            print >> sys.stderr, "Number of synteny blocks = %s" % nbOfSbs
            nbSbs_for[(tgm, gm)] = nbOfSbs

            if nbOfSbs > 0:
                sbLens = sorted([len(sb.la) for (_, sb) in sbsInPairComp.iteritems2d()])
                nbSbsOfLen = collections.Counter(sbLens)
                distribSbLens = [" %s:%s" % (nbSbsOfLen[sbLen], sbLen) for sbLen in sorted(nbSbsOfLen.keys())]
                distribSbLens = distribSbLens[:5] + ["..."] + distribSbLens[-3:]
                print >> sys.stderr, "Distribution of sb lengths (nbSbs:length) = %s" % ",".join(distribSbLens)
            wa = myMaths.myStats.getWeightedAverage(sbLens)
            wa_for[(tgm, gm)] = wa
            # save this to avoid executing all when it is not necessary
            pickle.dump(nbSbs_for, open('nbSbs_WithoutMonoSbs.txt', 'wb'))
else:
    nbSbs_for = pickle.load(open('nbSbs.txt', 'rb'))

X = []
Y = []
Z2 = []
for ((tgm, gm), nbSbs) in nbSbs_for.iteritems():
    X.append(float(gm))
    Y.append(float(tgm))
    Z2.append(float(nbSbs))
#Z1 = []
# for ((tgm, gm), wa) in wa_for.iteritems():
#     X.append(float(gm))
#     Y.append(float(tgm))
#     # Z1.append(float(wa))


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
# fig1 = plt.figure()
# ax1 = fig1.add_subplot(111, projection='3d')
# ax1.scatter(X, Y, zs=Z1, color='b')

fig2 = plt.figure()
ax2 = fig2.add_subplot(111, projection='3d')
ax2.scatter(X, Y, zs=Z2, color='b')
plt.show()

X3 = []
Z3 = []
for ((tgm, gm), nbSbs) in nbSbs_for.iteritems():
    if gm == 5:
        X3.append(tgm)
        Z3.append(nbSbs)
XZ3 = zip(X3, Z3)
XZ3 = sorted(XZ3, key=lambda x:x[0])
X3 = [tgm for (tgm, nbSbs) in XZ3]
Z3 = [nbSbs for (tgm, nbSbs) in XZ3]
plt.figure()
plt.plot(X3, Z3, color='black', linestyle='-', label='#sbs')
plt.ylabel('#sbs')
plt.xlabel('tandemGapMax')

# X4 = []
# Z4 = []
# for ((tgm, gm), wa) in wa_for.iteritems():
#     if gm == 5:
#         X4.append(tgm)
#         Z4.append(wa)
# XZ4 = zip(X4, Z4)
# XZ4 = sorted(XZ4, key=lambda x:x[0])
# X4 = [tgm for (tgm, wa) in XZ4]
# Z4 = [wa for (tgm, wa) in XZ4]
# plt.figure()
# plt.plot(X4, Z4, color='black', linestyle='-', label='#sbs')
# plt.ylabel('wa')
# plt.xlabel('tandemGapMax')

# plt.figure()
# plt.plot(wa_for_tgm.keys(), wa_for_tgm.values(), color='black', linestyle='-', label='wa')
# plt.ylabel('weighted average of the distribution of sbs lengths')
# plt.xlabel('tandemGapMax')
