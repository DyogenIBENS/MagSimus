#!/usr/bin/python
# -*- coding: utf-8 -*-
# requires python >2.6
import utils.myTools as myTools
import utils.myLightGenomes as myLightGenomes

import sys
from matplotlib import pyplot as plt

__doc__ = """ draw chromosome lengths of real and simulated genomes """

# test these
arguments = myTools.checkArgs( \
[("realGenome", file)], [("in:simulatedGenomes", str, "MagSimus/Data69/%s/Genes/genes.Homo.sapiens.list.bz2"), \
#("ICfactorOfSigma", float, 1.96) # Au cas ou on veut l'image avec l'IC grise
("limitOfDrawnChrs", int, 23)],
 __doc__)

genome1 = myLightGenomes.LightGenome(arguments["realGenome"])
genome1.removeUnofficialChromosomes()
chrLengthsRG = [len(chrom) for chrom in genome1.values()]
chrLengthsRG.sort(reverse=True)
if len(chrLengthsRG) > arguments["limitOfDrawnChrs"]:
    chrLengthsRG = chrLengthsRG[:arguments["limitOfDrawnChrs"]]

chrLengthsSGs = []
for i in range(100):
    genome2 = myLightGenomes.LightGenome(arguments["in:simulatedGenomes"] % i)
    chrLengthsSG = [len(chrom) for chrom in genome2.values()]
    chrLengthsSG.sort(reverse=True)
    if len(chrLengthsSG) > arguments["limitOfDrawnChrs"]:
        chrLengthsSG = chrLengthsSG[:arguments["limitOfDrawnChrs"]]
    if len(chrLengthsSG) < arguments["limitOfDrawnChrs"]:
        print >> sys.stderr, "Le genome numero %s a uniquement %s chromosomes" % (i, len(chrLengthsSG))
    chrLengthsSGs.append(chrLengthsSG)
# chrLengthsSGs = [ chrLengthsSG1, chrLengthsSG2, ...]
# chrLengthsSG = [tailleEnGenesDu1erChr, tailleEnGenesDu2emeChr, ... ] avec les chr par taille decroissante

# Renvoit une liste des chromosomes d'indice i dans un karyotype trie par taille decroissante de taille arguments["limitOfDrawnChrs"]
def chrLengthsOf_I_th_chrs(i, chrLengthsSGs):
    assert i < arguments["limitOfDrawnChrs"] # par construction
    listOfChrsOfLengthi=[]
    for chrLengthsSG in chrLengthsSGs:
        # Liste des tailles des chromosomes du karyotype classe par ordre decroissant
        assert len(chrLengthsSG) <= arguments["limitOfDrawnChrs"]
        if len(chrLengthsSG) > i: # Au cas ou il n'y a pas suffisament de chromosomes deans le karyotype
            listOfChrsOfLengthi.append(chrLengthsSG[i])
    return listOfChrsOfLengthi

data=[]
maxChrLength=0
for i in range(arguments["limitOfDrawnChrs"]):
    ys = chrLengthsOf_I_th_chrs(i, chrLengthsSGs)
    maxChrLength = max(max(ys), maxChrLength)
    data.append(ys)
bp = plt.boxplot(data, patch_artist=True) # patch_artist permet de colorier la box apres
plt.setp(bp['boxes'],facecolor='grey', alpha=0.5)
plt.setp(bp['fliers'], marker='None')
# Hide outlier points
plt.setp(bp['whiskers'], linestyle='-', color='k')
plt.setp(bp['caps'], color='k')
ymax=bp['caps'][0].get_data()[1][0]
plt.ylim(0,ymax+ymax*0.075)

#import numpy
#YSGmeans=[]
#YSGmins=[]
#YSGmaxs=[]
#YSGlows=[]
#YSGhighs=[]
#for i in range(arguments["limitOfDrawnChrs"]):
#       ymax= 0
#       ymin= sys.maxint
#       ys = chrLengthsOf_I_th_chrs(i, chrLengthsSGs)
#       n = len(ys)
#       ymean = numpy.mean(ys)
#       if n == 1:
#               ysigma = 0
#       else:
#               ysigma = numpy.sqrt((1.0/(n-1)) * sum([(y-ymean)*(y-ymean) for y in ys]))
#       for y in ys:
#               ymax = ymax if ymax > y else y
#               ymin = ymin if ymin < y else y
#       ylow =ymean-arguments["ICfactorOfSigma"]*ysigma
#       yhigh=ymean+arguments["ICfactorOfSigma"]*ysigma
#
#       YSGmeans.append(ymean)
#       YSGmins.append(ymin)
#       YSGmaxs.append(ymax)
#       YSGlows.append(ylow)
#       YSGhighs.append(yhigh)
#plt.fill_between(range(1,arguments["limitOfDrawnChrs"]+1), YSGlows, YSGhighs, facecolor=(0.45, 0.45,0.45), linewidth=0.5, label="Range of Simulations")
#plt.plot(range(1,arguments["limitOfDrawnChrs"]+1), YSGmeans, marker='o', linestyle='', color='k', label='Means of Simulations')

# show all values of simulations in green
# for i in range(arguments["limitOfDrawnChrs"]):
#       chrLengthsOf_I_th_chrs_ = chrLengthsOf_I_th_chrs(i  , chrLengthsSGs)
#       plt.plot([i+1 for _ in chrLengthsOf_I_th_chrs_], chrLengthsOf_I_th_chrs_, marker='o', linestyle='', color='g', markersize=5)

plt.plot(range(1,len(chrLengthsRG)+1), chrLengthsRG, marker='o', linestyle='', color='k', label='Real value', markeredgewidth=1, markersize=3)
plt.title('Lengths of the %s longest chromosomes' % arguments["limitOfDrawnChrs"])
plt.xlabel('chromosomes sorted by length')
plt.ylabel('lengths of chromosomes (in nb of genes)')
plt.ylim(ymin=0)
plt.legend()
#plt.show()
plt.savefig(sys.stdout, format='svg')
