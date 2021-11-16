#!/usr/bin/python
# -*- coding: utf-8 -*-

import utils.myTools as myTools
import utils.myDiags as myDiags
import utils.myMaths as myMaths

import enum

import sys
import matplotlib.pyplot as plt
import collections

__doc__ = """ Compute a graph N50/WA of simulated genomes with N50 of real genomes """

mode = enum.Enum('N50','WA')
modeList = list(mode._keys)
arguments = myTools.checkArgs([("speciesCombi", str)],
                              [("realDiags", str, "%s_vs_%s_f%s.sbs"),
                               ("simuDiags", str, "simu%s/%s_vs_%s_f%s.sbs"),
                               ('mode', mode, modeList)],
                              __doc__)

f = open(arguments['speciesCombi'])
nbPairWiseComp = counter = sum(1 for line in f)
f.close()

#@myTools.tictac
def N50andWaOfSbsOfPairwiseComp(sbsFileName, speciesCombi):
    N50s = collections.defaultdict(int)
    WAs = collections.defaultdict(float)
    for line in open(speciesCombi):
        [sp1, sp2, anc] = line.split()
        # paiwrise comparison
        pComp = sp1 + " " + sp2
        #print >> sys.stderr, sbsFileName % (sp1, sp2, anc)
        sbs = myDiags.parseSbsFile(sbsFileName % (sp1, sp2, anc))
        sbsLengths = [len(sb.la) for (_, sb) in sbs.iteritems2d()]
        N50 = myMaths.myStats.getValueNX(sbsLengths, 50)
        weightedAverage = myMaths.myStats.getWeightedAverage(sbsLengths)
        N50s[pComp] = N50
        WAs[pComp] = weightedAverage
    return (N50s, WAs)

(realN50s, realWAs) = N50andWaOfSbsOfPairwiseComp(arguments['realDiags'],
                                                  arguments['speciesCombi'])

simuN50ss = collections.defaultdict(list)
simuWAss = collections.defaultdict(list)
for i in range(100):
    print >> sys.stderr, "%d%%" % (i+1),
    (simuN50s, simuWAs) = N50andWaOfSbsOfPairwiseComp(arguments['simuDiags'] % (i,'%s','%s','%s'),
                                                      arguments['speciesCombi'])
    simuN50ss[i] = simuN50s
    simuWAss[i] = simuWAs
print >> sys.stderr, ''

assert len(simuN50ss.values()) == 100
assert len(simuN50ss[0].values()) == nbPairWiseComp, "%s %s" % (len(simuN50ss[0].values()), nbPairWiseComp)

xNames = []
s_Yss = []
r_Ys = []
for line in open(arguments['speciesCombi']):
    [sp1, sp2, anc] = line.split()
    # paiwrise comparison
    pComp = sp1 + " " + sp2
    if arguments['mode'] == 'N50':
        r_Ys.append(realN50s[pComp])
        s_Yss.append([N50s[pComp] for (simui, N50s) in simuN50ss.iteritems()])
    elif arguments['mode'] == 'WA':
        r_Ys.append(realWAs[pComp])
        s_Yss.append([WAs[pComp] for (simui, WAs) in simuWAss.iteritems()])
    commonName = {'Homo.sapiens': 'Human',
                  'Mus.musculus': 'Mouse',
                  'Anolis.carolinensis': 'Lizard',
                  'Gallus.gallus': 'Chicken',
                  'Monodelphis.domestica': 'Opossum'}
    sp1 = commonName[sp1]
    sp2 = commonName[sp2]
    pComp = sp1 + "_vs_" + sp2
    xNames.append(pComp)

len(r_Ys) == nbPairWiseComp, "len(r_Ys)=%s" % len(r_Ys)
len(s_Yss) == nbPairWiseComp, "len(s_Yss)=%s" % len(s_Yss)

y_max = 0
for ys in s_Yss:
    y_max = max(ys) if max(ys) > y_max else y_max
for y in r_Ys:
    y_max = y if y > y_max else y_max

(fig, ax1) = plt.subplots(figsize=(10, 6))
bp = plt.boxplot(s_Yss)
ax1.set_xlabel('pairwise comparisons')
if arguments['mode'] == 'N50':
    ax1.set_title('N50 of synteny blocks')
    ax1.set_ylabel('N50')
elif arguments['mode'] == 'WA':
    ax1.set_title('WA of synteny blocks')
    ax1.set_ylabel('WA')


#little trick to obtain the x values of the middle of the boxes
X = []
for i in range(len(bp['medians'])):
    xs = bp['medians'][i].get_xdata()
    xaverage = sum(xs)/len(xs)
    X.append(xaverage)

#Set the axes ranges and axes labels
numBoxes = len(s_Yss)
ax1.set_xlim(0.5, numBoxes+0.5)
top = y_max + (10.0 * y_max)/100.0
bottom = 0
ax1.set_ylim(bottom, top)
xtickNames = plt.setp(ax1, xticklabels=xNames)
# 'ha', horizontal alignment either : left, center or right
plt.setp(xtickNames, rotation=45, fontsize=8, ha='right')

#overplot real values
plt.plot(X, r_Ys, color='w', marker='*', markeredgecolor='k', linestyle='None')

# This crops perfectly the figure
plt.tight_layout()
# Show image
# plt.show()
# Save image
plt.savefig(sys.stdout, format='svg')
