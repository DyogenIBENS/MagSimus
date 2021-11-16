#!/usr/bin/python
# -*- coding: utf-8 -*-
# requires python >2.6
import utils.myTools as myTools
import utils.myLightGenomes as myLightGenomes

import sys
from matplotlib import pyplot as plt

__doc__ = """ compare (distributions of) chromosomes lengths of two genomes """

# test these
arguments = myTools.checkArgs( \
[("genome1", file), ("genome2", file)], [("limitOfDrawnChrs", str, 'None')],
 __doc__)

limitOfDrawnChrs = arguments['limitOfDrawnChrs']
if limitOfDrawnChrs == 'None':
    limitOfDrawnChrs = None
else:
    limitOfDrawnChrs = int(limitOfDrawnChrs)
genome1 = myLightGenomes.LightGenome(arguments["genome1"])
genome1.removeUnofficialChromosomes()
chrLengths1 = [len(chrom) for chrom in genome1.values()]
chrLengths1.sort(reverse=True)
if len(chrLengths1) > limitOfDrawnChrs:
    chrLengths1 = chrLengths1[:limitOfDrawnChrs]
if len(chrLengths1) < limitOfDrawnChrs:
    print >> sys.stderr, "genome1 %s has only %s chromosomes" % len(chrLengths1)

genome2 = myLightGenomes.LightGenome(arguments["genome2"])
chrLengths2 = [len(chrom) for chrom in genome2.values()]
chrLengths2.sort(reverse=True)
if len(chrLengths2) > limitOfDrawnChrs:
    chrLengths2 = chrLengths2[:limitOfDrawnChrs]
if len(chrLengths2) < limitOfDrawnChrs:
    print >> sys.stderr, "genome2 %s has only %s chromosomes" % len(chrLengths2)

plt.plot(range(1,len(chrLengths1)+1), chrLengths1, marker='o', linestyle='', color='k', label='genome1', markeredgewidth=1, markersize=6)
plt.plot(range(1,len(chrLengths2)+1), chrLengths2, marker='o', linestyle='', color='r', label='genome2', markeredgewidth=1, markersize=6)
plt.title('Lengths of the %s longest chromosomes' % max([limitOfDrawnChrs, len(chrLengths1), len(chrLengths1)]))
plt.xlabel('chromosomes sorted by length')
plt.ylabel('lengths of chromosomes (#genes)')
plt.ylim(ymin=0)
plt.legend()
plt.show(block=True)
plt.savefig(sys.stdout, format='svg')
