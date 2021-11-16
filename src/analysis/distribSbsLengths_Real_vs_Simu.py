#!/usr/bin/python
# -*- coding: utf-8 -*-

import utils.myTools as myTools
import utils.myDiags as myDiags

# import matplotlib
from matplotlib import pyplot as plt
import collections
# import re #pour parser
# import os, sys
import sys
# import numpy as np

arguments = myTools.checkArgs([("Real_S1_vs_S2.sbs", file)],
                              [('in:Sim%s_S1_vs_S2.sbs', str, 'in:Sim%s_S1_vs_S2.sbs')],
                               __doc__)
# ("ICfactorOfSigma", float, 1.96) au cas ou on veut l'IC de confiance grise

#def distribSbsLengthsFromSbsFile(fileName):
#    f = open(fileName, 'r')
#    idxSb_old = 0
#    nbHps = 0
#    realDistribSbsLengths = collections.defaultdict(int)
#    for line in f.readlines():
#        idxSb = int(line.split('\t')[0])
#        if idxSb == idxSb_old:
#            nbHps += 1
#            continue
#        else:
#            realDistribSbsLengths[nbHps] += 1
#            nbHps = 0
#            idxSb_old = idxSb
#    f.close()
#    X = []
#    Y = []
#    for (nbHps, nbSbs) in realDistribSbsLengths.iteritems():
#        X.append(nbHps)
#        Y.append(nbSbs)
#    return (X, Y)
#


def distribSbsLengthsFromSbsFile(fileName):
    X = []
    Y = []
    distribSbsLengths = collections.defaultdict(int)
    sbs = myDiags.parseSbsFile(fileName)
    for ((c1, c2), sb) in sbs.iteritems2d():
        nbHps = len(sb.la)
        distribSbsLengths[nbHps] += 1
    for (nbHps, nbSbs) in distribSbsLengths.iteritems():
        X.append(nbHps)
        Y.append(nbSbs)
    return (X, Y)

(X1, Y1) = distribSbsLengthsFromSbsFile(arguments['Real_S1_vs_S2.sbs'])

X2s = []
Y2s = []
for i in range(100):
    (X2, Y2) = distribSbsLengthsFromSbsFile(arguments['in:Sim%s_S1_vs_S2.sbs'] % i)
    X2s.append(X2)
    Y2s.append(Y2)

# l= {sb length: [nb of sbs that have this length in simu 1, ...]}
l = collections.defaultdict(list)
assert len(X2s) == len(Y2s) == 100
for simui in range(100):
    X2 = X2s[simui]
    Y2 = Y2s[simui]
    assert len(X2) == len(Y2)
    for i in range(len(X2)):
        l[X2[i]].append(Y2[i])

#Calcul de l'intervalle de confiance pour les blocks des dists simules

#import numpy
#X2=[]
#Y2mean=[]
#Y2low=[]
#Y2high=[]
#Y2min=[]
#Y2max=[]
#for (i, (x,ys)) in enumerate(l.items()):
#       n=len(ys)
#       ymean = numpy.mean(ys)
#       if n == 1:
#               ysigma = 0
#       else:
#               ysigma = numpy.sqrt((1.0/(n-1)) * sum([(y-ymean)*(y-ymean) for y in ys]))
#       ylow = ymean - arguments["ICfactorOfSigma"]*ysigma
#       yhigh = ymean + arguments["ICfactorOfSigma"]*ysigma
#       X2.append(x)
#       Y2mean.append(ymean)
#       Y2low.append(ylow)
#       Y2high.append(yhigh)
#       Y2min.append(min(ys))
#       Y2max.append(max(ys))

#########################################################
#       An other way to compute mu and sigma            #
#########################################################
#
#from scipy.stats import norm
#import matplotlib.mlab as mlab
# best fit of data
#(mu, sigma) = norm.fit(l[5])
#print >> sys.stderr, "mu =", mu
#print >> sys.stderr, "sigma =",sigma
# # the histogram of the data
#n, bins, patches = plt.hist(l[5], 30, normed=1, facecolor='green', alpha=0.75)
#
# # add a 'best fit' line
#y = mlab.normpdf( bins, mu, sigma)
#plt.plot(bins, y, 'r--', linewidth=2)
#
#plt.xlabel('Smarts')
#plt.ylabel('Probability')
#plt.title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=%.3f,\ \sigma=%.3f$' %(mu, sigma))
#plt.grid(True)
#plt.show()
##fig1.savefig("/tmp/tmp.svg", format='svg')
#sys.exit(0)

# preparation du boxplot
data = []
for (i, (x, ys)) in enumerate(l.items()):
    data.append(ys)

#print >> sys.stderr, re.split('\[|\]|,', arguments["yBounds"])[1:-1]a

#upper.fill_between(X2, Y2low, Y2high, color=(0.45, 0.45, 0.45), linewidth=0.5)
#lower.fill_between(X2, Y2low, Y2high, color=(0.45, 0.45, 0.45), linewidth=0.5)
#upper.plot(X2, Y2mean, marker='o', linestyle='', color='k', label='Means of Simulations', markersize=5)
#lower.plot(X2, Y2mean, marker='o', linestyle='', color='k', label='Means of Simulations', markersize=5)
bp = plt.boxplot(data, positions=[x for (x, _) in l.items()], patch_artist=True)
plt.setp(bp['fliers'], marker='None')
# Hide outlier points
plt.setp(bp['boxes'], facecolor='grey', alpha=0.5)
plt.setp(bp['whiskers'], linestyle='-', color='k')
plt.setp(bp['caps'], color='k')
#ymax=bp['caps'][0].get_data()[1][0]
#plt.ylim(0,ymax+ymax*0.075)

#xtickNames = plt.setp(ax1, xticklabels=np.repeat(randomDists, 2))
ymax = bp['caps'][0].get_data()[1][0]

#print >> sys.stderr, [2] + np.arange(5, max(X2)+1, 5.0)
#print >> sys.stderr, [2] + range(5, max(X2)+1, 5)
plt.xticks([2] + range(5, int(max(X2))+1, 5))
#ax.set_xticks([x for x in X2 if x%5 == 0], [x for x in X2 if x%5 == 0])
#xtickNames = plt.setp(ax, xticklabels=[x for x in X2 if x%5 == 0])
#plt.setp(xtickNames)

#plt.setp(lower.get_xticklabels(), rotation='vertical', fontsize=14) # Works

#plt.xticks([i for (i,_) in enumerate(X2)], [x for (_,x) in enumerate(X2)])
#facecolors='none' allow to keep empty circles for the plotting
plt.plot(X1, Y1, marker='o', linestyle='', color='w', label='Real value', markersize=4, markeredgewidth=1.5, markeredgecolor='k')

plt.ylabel('Nb of synteny blocks')
plt.xlabel('Lengths of synteny blocks')
#plt.tight_layout()
plt.xlim(0, 100)
#plt.title("OSB length distribution in Human-Mouse comparison \n Confidence Interval : %s*sigma around mean" % arguments["ICfactorOfSigma"] )
plt.title("Distribution of the lengths of synteny blocks\nHuman-Mouse comparison")
plt.legend()
#plt.show()
plt.savefig(sys.stdout, format='svg')
