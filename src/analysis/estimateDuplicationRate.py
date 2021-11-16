#!/usr/bin/python

# MagSimus/src/analysis
# python 2.7.8
# Author: Lucas TITTMANN, except for the function "nans". There, the source is written within the function.
# Copyright (c) 2015 IBENS/Dyogen Lucas TITTMANN and Hugues ROEST CROLLIUS
# mail : hrc(at)ens.fr
# Licence: GPL 3.0, except for function "nans"
# This is free software, you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 (GPL v3) or later and the CeCiLL v2 license in France

# This script simulates gene events and tries to asess numerically the underlying gene event rates
# from the gene tree estimates

from collections import Counter
from collections import defaultdict
import numpy as np
import scipy.optimize as sco

def simulateRandomGeneEvents(dList, branchLength, nStartGenes, dt = 0.001):
    [dd, de, dn] = dList

    nDup = 0
    t = 0
    geneBirths = [t + np.random.exponential(dn)]
    geneErased = []
    geneRobust = nStartGenes*[True]

    recGeneName = nStartGenes
    genes = []
    for gene in xrange(nStartGenes):
        genes.append({'name': gene,
                      'tDel': t + np.random.exponential(de),
                      'tDup': t + np.random.exponential(dd)
        })

    while t <= branchLength:
        while t > geneBirths[-1]:
            geneBirths.append(geneBirths[-1]+np.random.exponential(dn))
            genes.append({'name': recGeneName,
                          'tDel': geneBirths[-1] + np.random.exponential(de),
                          'tDup': geneBirths[-1] + np.random.exponential(dd)
            })
            recGeneName += 1

        for (iGene, gene) in enumerate(genes):
            if t > gene['tDup']:
                nDup += 1
                genes.append({'name': gene['name'],
                              'tDel': gene['tDup'] + np.random.exponential(de),
                              'tDup': gene['tDup'] + np.random.exponential(dd)
                })
                genes[iGene]['tDup'] = gene['tDup'] + np.random.exponential(dd)
                if iGene < nStartGenes:
                    geneRobust[iGene] = False
            elif t > gene['tDel']:
                geneErased.append(gene)
                del genes[iGene]
                if iGene < nStartGenes:
                    geneRobust[iGene] = False

        t += dt

    geneFamilyCopyNumbers = Counter([gene['name'] for gene in genes])
    statistics = {'nErased': len(geneErased),
                  'nErasedObs': len(set([gene['name'] for gene in geneErased if gene['name'] < nStartGenes])),
                  'nRobust': sum(geneRobust),
                  'nRobustObs': Counter(geneFamilyCopyNumbers.itervalues())[1],
                  'nDuplication': nDup,
                  'nDuplicationObs': sum([recCopyNumber-1 for (recName, recCopyNumber) in geneFamilyCopyNumbers.iteritems() if recName < nStartGenes]),
                  'nNewGeneNeverSeen': len([recName for recName in xrange(nStartGenes, recGeneName) if recName not in geneFamilyCopyNumbers.keys()])
                  }
    return(statistics)

branchLength = 1.0
nStartGenes = 20000
reference = {'nRobustObs': 14000, 'nDuplicationObs': 5000, 'nErasedObs': 2000}
dt = 0.001

dd = 4
de = 10
dn = 0.001

def calcScore(dList, branchLength=branchLength, nStartGenes=nStartGenes, dt=dt, reference = reference, nRep = 100):
    meanStat = defaultdict(float)
    for i in xrange(nRep):
        newStat = simulateRandomGeneEvents(dList, branchLength=branchLength, nStartGenes=nStartGenes, dt = dt)
        for (statName, statValue) in newStat.iteritems():
            meanStat[statName] += float(statValue) / nRep

    sumRob = meanStat['nRobustObs'] - reference['nRobustObs']
    sumDup = meanStat['nDuplicationObs'] - reference['nDuplicationObs']
    sumEra = newStat['nErasedObs'] - reference['nErasedObs']
    print(meanStat)
    return(sumRob+sumDup+sumEra)

sco.minimize(calcScore, [dd, de, dn], args=(branchLength, nStartGenes, dt, reference),
             bounds=[(0.1, nStartGenes), (0.1, nStartGenes), (0.0000001, 10)], options={'maxiter': 1000, 'disp': True})

from sympy import *
from sympy.assumptions.assume import global_assumptions
x = Symbol('x', real=True, nonnegative=True)
k = Symbol('k', integer = True, nonnegative=True)
n = Symbol('n', integer = True, nonnegative=True)
de = Symbol('de', real=True, positive=True)
dd = Symbol('dd', real=True, positive=True)
l = Symbol('l', real=True, positive=True)
global_assumptions.add(Q.nonnegative(x))
global_assumptions.add(Q.nonnegative(k))
global_assumptions.add(Q.integer(k))
global_assumptions.add(Q.nonnegative(n))
global_assumptions.add(Q.integer(n))
global_assumptions.add(Q.positive(l))
global_assumptions.add(Q.positive(de))
global_assumptions.add(Q.positive(dd))
integrate(de/factorial(k) * x**k * exp(-x*(1+de)), x)

integrate(x**k * exp(x*(dd-de-1) * (l-x)**(n-1) * ((dd-de)*(l-x)-n)), x)