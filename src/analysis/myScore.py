#!/usr/bin/python

# MagSimus/src/analysis
# python 2.7.8
# Author: Lucas TITTMANN
# Copyright (c) 2015 IBENS/Dyogen Lucas TITTMANN, Joseph LUCAS and Hugues ROEST CROLLIUS
# mail : hrc(at)ens.fr
# Licence: GPL 3.0
# This is free software, you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 (GPL v3) or later and the CeCiLL v2 license in France

# DESCRIPTION
#
# This library contains some classes and functions
# They are used in example code in the: if __name__ == '__main__':
# # Concerning the 3 classes: 1) myScore, 2) myScoreComparison, 3) myScoreComparisonArchive
# 1) For the construction of a score, one needs at least the modernGenomes and ancGeneFiles of a speciesTree.
#   This will create a minimal version of a score, and it is usually used for real data (as we don't have
#   any other data there). For simulated data, one can also specify the simulated ancestral genomes and
#   the stderr of MagSimus, from which the simulated number of chromosomal rearrangements is read.
#   This helps evaluating how far our observation (with PhylDiag and other tools) is away from the
#   'real' (in this case, simulated reality) values.
# 2) The class myScoreComparison takes 2 scores (usually one from simulated data and one from real data)
#   and calculates a couple of comparison statistics, like differences (sim - real), ratios (sim / real),
#   in case of distributions e.g. Kolmogorov-Smirnoff and other distribution comparisons. The names of the
#   class variables in myScoreComparison are the same as in myScore, i.e. gMeasure remains gMeasure,
#   but in the myScoreComparison instead of a number it contains a dictionary with 'Diff' (for difference)
#   and 'Ratio' (for simulated_gMeasure / real_gMeasure). In some cases, as with the distributions, where
#   the difference has not that much information, new variables are created (e.g. chrLengths contains the
#   difference between the two chromosomeSize distributions, but also additional information about
#   intersections, slope of the difference, p-value of the difference, etc. Another way to look at distri-
#   bution differences is the Kolomogorov-Smirnoff statistic, which is calculated in chrLengths_KS)
# 3) Finally, myScoreComparisons can be stored in the myScoreComparisonArchive. It fulfills two objectives:
#   a) It stores most of the information of myScoreComparison in self.archive. Usually, all simulations
#       based on the same specRates/parameterFile are called repetitions, and all repetitions have the same ID.
#       Where in the myScore/myScoreComparison there are class variables, e.g. chrLengths, there is a key in the
#       self.archive dictionary, e.g. value = self.archive['chrLengths'][mySpecies][myComparisonID][NumberOfRepetition]
#   b) With the help of myScoreComparison.calcStatistics() summary statistics for all the ScoreComparisons
#       are calculated. Afterwards, graphics can be drawn with plotScoreDistribution.
#   In case of many ScoreComparisons to be stored, one can break up the archive into serveral archives and calculate
#   metaStatistics with myScoreComparisonArchive.calcMetaScore().
import collections

import scipy.optimize as sco
import scipy.stats as scs
import scipy.interpolate as sci
import numpy as np
import math
import sys
import os
import time
import copy
import itertools
import cPickle as pickle
from collections import defaultdict
from collections import Counter
from collections import OrderedDict

#import myOptimizerGraphics
from utils import myMaths
import utils.myDiags as myDiags
import utils.myLightGenomes as myLightGenomes
import utils.myMapping as myMapping
import utils.myPhylTree as myPhylTree
import utils.myFile as myFile
import matplotlib.pyplot as plt
import libs.myBreakpointsAnalyser as myBA
import utils.myIntervals as myIntervals
import utils.myTools as myTools

def diffAndRatio(value, valueRef):
    diff = value - valueRef
    ratio = myScore.calcRatio(value, valueRef, 1)
    return {'Diff': diff, 'DiffRatio': ratio}

def extrap1d(interpolator):
    """
    A very simple interpolator that use a linear extrapolation from border values
    from: http://stackoverflow.com/questions/2745329/how-to-make-scipy-interpolate-give-an-extrapolated-result-beyond-the-input-range
    Example:
        from scipy.interpolate import interp1d
        from scipy import arange, array, exp
        x = arange(0,10)
        y = exp(-x/3.0)
        f_i = interp1d(x, y)     //this interpolates the exp func on the range of x
        f_x = extrap1d(f_i)
    :param interpolator:
    :return: a function that can be applied on any new_x array
    """
    xs = interpolator.x
    ys = interpolator.y

    def pointwise(x):
        if x < xs[0]:
            return ys[0]+(x-xs[0])*(ys[1]-ys[0])/(xs[1]-xs[0])
        elif x > xs[-1]:
            return ys[-1]+(x-xs[-1])*(ys[-1]-ys[-2])/(xs[-1]-xs[-2])
        else:
            return interpolator(x)

    def ufunclike(xs):
        return np.array(map(pointwise, np.array(xs)))

    return ufunclike

def calcRankedValuesEmpFuncFromList(L, nbBins=None, normed=True, normedXaxis=True, increasingOrder=True):
    """
    Return the ranked list (possibly normed) of an array of equidistant sorted values (typically a karyotype).
    The x-axis (possibly normed) corresponds to the relative ranks of the values by specified order
    :param L: [..., value, ...]
    :return: (new_X, rankedValues)
    """
    assert isinstance(L, list)
    if not myTools.isSorted(L, increasingOrder=increasingOrder, stricly=False):
        L.sort(reverse=(not increasingOrder))
    Y = np.array(L)
    if nbBins is None:
        nbBins = len(Y)
    X = np.arange(0, len(Y), 1)
    if normedXaxis:
        new_X = np.arange(0, 1, 1.0 / float(nbBins))
        new_Y = extrapolate(zip(X, Y), new_X * len(Y))
    else:
        new_X = np.arange(0, len(Y), len(Y) / float(nbBins))
        new_Y = extrapolate(zip(X, Y), new_X)

    if normed:
        new_Y /= float(sum(new_Y))

    rankedValues = new_Y
    #cdf = np.cumsum(ranked)
    if nbBins is None:
        assert len(new_X) == len(L)
        if normedXaxis is False:
            all([x[i+1] - x[i] for (i, x) in enumerate(new_X)[:-1]])
    return (new_X, rankedValues)

def plotRankedValuesFromList(L, nbBins=None, normed=False, normedXaxis=False, increasingOrder=False):
    (x, rv) = calcRankedValuesEmpFuncFromList(L, nbBins=nbBins, normed=normed, normedXaxis=normedXaxis, increasingOrder=increasingOrder)
    plt.figure()
    plt.plot(x, rv)
    plt.show(block=False)

def plotRankedValuesFromTwoLists(L, Lref, nbBins=None, normed=False, normedXaxis=False, increasingOrder=False):
    (xref, rvref) = calcRankedValuesEmpFuncFromList(Lref, nbBins=nbBins, normed=normed, normedXaxis=normedXaxis, increasingOrder=increasingOrder)
    (x, rv) = calcRankedValuesEmpFuncFromList(L, nbBins=len(xref), normed=normed, normedXaxis=normedXaxis, increasingOrder=increasingOrder)
    plt.figure()
    plt.plot(xref, rvref, color='k', marker='o')
    plt.plot(xref, rv, color='b', marker='o')
    plt.show(block=False)

def extrapolate(empfuncNeedingExtrapolation, Xref):
    """
    interpolates empfuncNeedingExtrapolation to have values at the x of distribRef
    """
    Xi = []
    Yi = []
    for (xi, yi) in empfuncNeedingExtrapolation:
        Xi.append(xi)
        Yi.append(yi)
    # interpolate
    f_i = sci.interp1d(Xi, Yi, kind='linear')
    # extrapolate
    f_x = extrap1d(f_i)
    newYi = f_x(Xref)
    # normalize
    newYi /= sum(newYi)
    return newYi

class myScoreComparisonArchive:

    def calcBestScore(self, score, species):
        return self.rankList[score][species][0]

    def calcRankOfCompID(self, scoreName, spOrCombo, compID):
        for (iEntry, entry) in enumerate(self.rankList[scoreName][spOrCombo]):
            if compID == entry['ID']:
                return (iEntry, entry)

    @staticmethod
    def createFromDir(path, onlySim='4', metaDir=False):
        """
        :param path: either a str with the path, or a list of paths
        :param filterBy:
        :return:
        """

        def loadFromFolderContainingOpts(folderContainingOpts, onlySim=onlySim, compPrefix=''):
            for simIdx in os.listdir(folderContainingOpts):
                if simIdx == onlySim:
                    simPath = folderContainingOpts + '/' + simIdx
                    for repetIdx in os.listdir(simPath):
                        if repetIdx.isdigit():
                            repetPath = simPath + '/' + repetIdx
                            for fileName in os.listdir(repetPath):
                                if fileName.endswith((".comp")):
                                    pathToComp = repetPath + '/' + fileName
                                    print >> sys.stderr, 'Loading comp %s' %  pathToComp
                                    # comp.ID = iSim-iRepet-real
                                    comp = myScoreComparison.load(pathToComp)
                                    # DEBUG
                                    comp.ID = str(simIdx) + '-' + str(repetIdx) + '-' + 'real'
                                    print >> sys.stderr, 'Adding comp %s (compID=%s) to the archive with prefix %s' %  (pathToComp, comp.ID, compPrefix)
                                    mca.addComp(comp, compPrefix)

        mca = myScoreComparisonArchive()

        # recursively iterate through a directory and all its subdirectories
        #for pathToFolder, dirsInFolder, filesInFolder in os.walk(path):
        if not metaDir:
            loadFromFolderContainingOpts(path, onlySim=onlySim, compPrefix='')
        else:
            for fileOrFolder in os.listdir(path):
                folderPathContainingSims = path + '/' + fileOrFolder
                loadFromFolderContainingOpts(folderPathContainingSims, onlySim=onlySim, compPrefix=fileOrFolder)

            #
            #
            # if filterBy in pathToFolder:
            #     for fileName in filesInFolder:
            #         if fileName.endswith((".comp")):
            #             pathToComp = pathToFolder + '/' + fileName
            #             lastFolderNameInPath = path.split('/')[-1] if path.split('/')[-1] != '' else path.split('/')[-2]
            #             # index of the first component not in path
            #             idxExtOfPath = pathToComp.split('/').index(lastFolderNameInPath) + 1
            #             compNamePrefix = '-'.join(pathToComp.split('/')[idxExtOfPath:-1])
            #             print >> sys.stderr, 'Loading comp %s' %  pathToComp
            #             comp = myScoreComparison.load(pathToComp)
            #             print >> sys.stderr, 'Adding comp %s to the archive with prefix %s' %  (pathToComp, compNamePrefix)
            #             mca.addComp(comp, compNamePrefix)
            # #else:
            #print >> sys.stderr, '%s ignored since the filter=%s is not in' % (pathToFolder, filterBy)
        mca.calcStatistics()
        print >> sys.stderr, 'MCA created and statistics calculated !'
        return mca

    def cumulateScoreForRanklist(self, scoreCompsBySpOrCombo, relative=False):
        # This method cumulates scores within the rankList over all species or combos of species
        # If relative == False, then the statistics are added (for absolute fitting) arithmetic mean
        # If relative == True, then the RatioToBestStatistics are multiplied (for relative fitting) geometrical mean

        # scoreCompsBySpOrCombo[spOrCombos] = [..., {'ID': 'simuIdx-real', 'Stat': mean, 'RatioToBest': ratio}, ...]
        if relative:
            cumDict = {compScoreID: 1 for compScoreID in self.compIDs}
        else:
            cumDict = {compScoreID: 0 for compScoreID in self.compIDs}
        numberOfCumulations = 0

        for spOrCombos, listOfCompScores in scoreCompsBySpOrCombo.iteritems():
            numberOfCumulations += 1
            for compDict in listOfCompScores:
                # 'ID' is the scoreComp Name : 'simuIdx-real'
                if relative:
                    cumDict[compDict['ID']] *= compDict['RatioToBest']
                else:
                    cumDict[compDict['ID']] += compDict['Stat']

        if relative:
            # If 10 values were multiplied, hence take 10th root of cumulated value
            cumList = [(compID, np.power(val, 1.0/numberOfCumulations)) for (compID, val) in cumDict.iteritems()]
        else:
            # If 10 values were added, hence cumulated value/10
            cumList = [(compID, myScore.calcRatio(val, numberOfCumulations)) for (compID, val) in cumDict.iteritems()]
        cumList.sort(key=lambda x: x[1])

        return cumList

    def calcRankList(self, scoreIDsSelection=[], calcCumulatedScore=True, valueIfDivisionByZero=1000):
        # scoreCompID si 'simuIdx-real'
        # self.scores[scoreName][spOrCombo][scoreCompID] = {'min': X, 'max': X, 'mean': X, 'sd': X, 'n': X}
        self.rankList = {}
        if len(scoreIDsSelection) == 0:
            scoreIDsSelection = self.scores.iterkeys()
        for scoreName in scoreIDsSelection:
            self.rankList[scoreName] = {}
            # self.scores[scoreName][spOrCombo].update(self.calcStatsForlist(vals))
            for (spOrCombos, scoresByIDs) in self.scores[scoreName].iteritems():
                rankedScoreComps = []
                for (scoreCompID, stats) in scoresByIDs.iteritems():
                    if not isinstance(stats, dict):
                        continue
                    if 'KS' in scoreName:
                        # KS is either LogStat or LogpValue, with Stat and pValue >0 and sometimes <1,
                        # hence building abs-value before sorting doesn't make sense
                        rankedScoreComps.append((scoreCompID, stats['mean']))
                    else:
                        rankedScoreComps.append((scoreCompID, np.abs(stats['mean'])))
                rankedScoreComps.sort(key=lambda x: x[1])
                # compute RatioToBest
                compScoresBySpOrCombo = []
                bestScore = rankedScoreComps[0][1]
                for (scoreCompID, mean) in rankedScoreComps:
                    if 'Log' in scoreName:
                        ratio = np.exp(mean - bestScore)
                    else:
                        ratio = myScore.calcRatio(mean, bestScore, valueIfDivisionByZero)
                    compScoresBySpOrCombo.append({'ID': scoreCompID, 'Stat': mean, 'RatioToBest': ratio})

                # self.rankList[scoreName][spOrCombos] = [..., {'ID': 'simuIdx-real', 'Stat': mean, 'RatioToBest': ratio}, ...]
                self.rankList[scoreName][spOrCombos] = compScoresBySpOrCombo

        if calcCumulatedScore:
            for (scoreName, compScoresBySpOrCombo) in self.rankList.iteritems():
                compScoresBySpOrCombo['Cumulated'] = {'Abs': self.cumulateScoreForRanklist(self.rankList[scoreName], relative=False),
                                                      'Rel': self.cumulateScoreForRanklist(self.rankList[scoreName], relative=True)}

    @staticmethod
    def _addObsToSD(nOld, meanOld, sdOld, newObs):
        out = float(nOld-1)/nOld * sdOld*sdOld + float(newObs - meanOld)*(newObs - meanOld)/(nOld+1)
        out = np.sqrt(out)
        return out

    @staticmethod
    def _combineTwoSDs(nX, meanX, sdX, nY, meanY, sdY):
        out = (nX-1)*sdX**2 + nX*meanX**2 + (nY-1)*sdY**2 + nY*meanY**2 - float(nX*meanX + nY*meanY)**2/(nX + nY)
        out = float(out)/(nX + nY - 1)
        out = np.sqrt(out)
        return out

    @staticmethod
    def combineStats(dict1, dict2):
        out = {}
        if dict1['n'] == 1:
            if dict2['n'] == 1:
                out['sd'] = float(dict1['mean'] - dict2['mean'])**2 / 2
            else:
                out['sd'] = myScoreComparisonArchive._addObsToSD(dict2['n'], dict2['mean'], dict2['sd'], dict1['mean'])
        else:
            if dict2['n'] == 1:
                out['sd'] = myScoreComparisonArchive._addObsToSD(dict1['n'], dict1['mean'], dict1['sd'], dict2['mean'])
            else:
                out['sd'] = myScoreComparisonArchive._combineTwoSDs(dict1['n'], dict1['mean'], dict1['sd'], dict2['n'], dict2['mean'], dict2['sd'])

        out['min'] = min(dict1['min'], dict2['min'])
        out['max'] = max(dict1['max'], dict2['max'])
        out['mean'] = (dict1['n'] * dict1['mean'] + dict2['n'] * dict2['mean']) / (dict1['n'] + dict2['n'])
        out['n'] = dict1['n'] + dict2['n']
        return out

    @staticmethod
    def calcMetaScore(listOfMCAs):
        mca = listOfMCAs.pop()
        if isinstance(mca, basestring):
            mca = myScoreComparisonArchive.load(mca)
        metaScores = copy.deepcopy(mca)

        stats = ['n', 'mean', 'max', 'min', 'sd']

        while len(listOfMCAs) > 0:
            mca = listOfMCAs.pop()
            if isinstance(mca, basestring):
                mca = myScoreComparisonArchive.load(mca)

            for (scoreName, specs) in metaScores.scores.iteritems():
                for (spec, comps) in specs.iteritems():

                    dictMeta = {stat: comps[stat] for stat in stats}
                    dictNew = {stat: mca.scores[scoreName][spec][stat] for stat in stats}
                    metaScores.scores[scoreName][spec].update(myScoreComparisonArchive.combineStats(dictMeta, dictNew))

                    compsBoth = defaultdict(int)
                    for comp in comps.iterkeys():
                        if comp not in stats:
                            compsBoth[comp] += 1

                    for comp in mca.scores[scoreName][spec].iterkeys():
                        if comp not in stats:
                            compsBoth[comp] += 2

                    for (comp, case) in compsBoth.iteritems():
                        # In case == 1, nothing has to be done
                        if case == 2:
                            # Case if comp exist only in new MCA
                            metaScores.scores[scoreName][spec][comp] = mca.scores[scoreName][spec][comp]
                        elif case == 3:
                            # Case if comp exists in both MCA and hence stats must be combined:
                            metaScores.scores[scoreName][spec][comp] = myScoreComparisonArchive.combineStats(metaScores.scores[scoreName][spec][comp],
                                                                                                      mca.scores[scoreName][spec][comp])
        return metaScores

    @staticmethod
    def plotScoreDistribution(scores, scoreName, markComps=['0-real'], consideredSpsOrCombos=[],
                              sortScoresByMeans=False, saveImageToPath='', helpLines='',
                              mainTitle=None, format='png'):
        """
        :param scores:
        :param scoreName:
        :param markComps: highlight some comps
        :param consideredSpsOrCombos: only consider consideredSpsOrCombos in the list
        :param sortScoresByMeans:
        :param saveImageToPath:
        :param helpLines:
        :param mainTitle:
        :param format:
        :return:
        """
        # score should be something like mca.scores
        # scoreName should be something like 'dist-nInv-DiffRatio'
        # to not markComps, simply put an empty list in there
        # helpLines can either be '' for nothing, 'meanCI' for confidence intervall of the mean, or 'minmax'

        from itertools import cycle
        lines = ["-", "--", "-."]
        if helpLines != '':
            lines = ["-"]
        colors = ["r", "g", "b", 'black', 'lightblue']
        linecycler = cycle(lines)
        colorcycler = cycle(colors)
        labels = []
        legendEntries = []
        fig = plt.figure(figsize = (16.5, 10.5))
        ax = fig.add_subplot(111)
        ax.get_yaxis().get_major_formatter().set_useOffset(False)


        for (spOrCombo, vals) in scores[scoreName].iteritems():
            if len(consideredSpsOrCombos) > 0 and spOrCombo not in consideredSpsOrCombos:
                continue
            xyStats = [(d['mean'], d['sd'], d['n'], d['min'], d['max'], simRepetID) for (simRepetID, d) in vals.iteritems() if isinstance(d, dict)]
            if sortScoresByMeans:
                xyStats.sort()
            else:
                # rank by simRepetID
                xyStats.sort(key=lambda x: x[-1])
                # FIXME: do this crapy code is usefull ?
                # # score by simRepetID
                # try:
                #     xMeanI = [(int(xyStat[-1].split('-')[-2]), xyStat[0], i) for (i, xyStat) in enumerate(xyStats)]
                # except ValueError:
                #     # xyStat[-1] = simRepetID
                #     # xyStats[0] = mean
                #     xMeanI = [(xyStat[-1], xyStat[0], i) for (i, xyStat) in enumerate(xyStats)]
                # xMeanI.sort()
                # xyStats = [xyStats[xyStats[2]] for xyStats in xMeanI]

            color = next(colorcycler)
            ys = [stats[0] for stats in xyStats]
            xs = [stats[-1] for stats in xyStats]
            for compName in markComps:
                #print >> sys.stderr, compName
                #print >> sys.stderr, vals
                compNames = [s for s in vals.iterkeys() if compName in s]
                #print >> sys.stderr, compNames
                if len(compNames) > 0:
                    ax.plot(ys.index(vals[compNames[0]]['mean']), vals[compNames[0]]['mean'], color = color, marker = 'o', markersize=4)

            thisPlot, = ax.plot(range(len(ys)), ys, next(linecycler), color = color, markersize=2)
            legendEntries.append(thisPlot)

            if helpLines.lower() != '':
                if helpLines.lower() == 'meanci':
                    upperLine = [(x, stats[0] + 1.96*stats[1]/np.sqrt(stats[2])) for (x, stats) in enumerate(xyStats) if stats[2]>1]
                    lowerLine = [(x, stats[0] - 1.96 * stats[1] / np.sqrt(stats[2])) for (x, stats) in enumerate(xyStats) if stats[2]>1]
                elif helpLines.lower() == 'minmax':
                    upperLine = [(x, stats[4]) for (x, stats) in enumerate(xyStats)]
                    lowerLine = [(x, stats[3]) for (x, stats) in enumerate(xyStats)]
                else:
                    raise ValueError('Wrong input in helpLines !')
                ax.plot([x[0] for x in upperLine], [x[1] for x in upperLine], linestyle = ':', color = color)
                ax.plot([x[0] for x in lowerLine], [x[1] for x in lowerLine], linestyle = ':', color = color)

            if isinstance(spOrCombo, tuple):
                spOrCombo = ' - '.join(spOrCombo)
            shortName = spOrCombo.split(' ')
            shortName = ''.join([word[0] for word in shortName]).upper()
            labels.append(shortName)

        if mainTitle is None:
            mainTitle = 'Distribution of ' + scoreName
            if helpLines != '':
                mainTitle += '\n(helplines indicate ' + helpLines + ')'
            if sortScoresByMeans:
                mainTitle += '\n(sorted by y-value)'
        ax.set_title(mainTitle, fontsize=20)
        ax.set_xlabel('Simulation', fontsize=16)
        ax.set_ylabel('Score', fontsize=16)
        ax.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.1,
                         box.width, box.height * 0.9])

        # Put a legend below current axis
        ax.legend(legendEntries, labels, loc='upper center', bbox_to_anchor=(0.5, -0.05),
                  fancybox = True, shadow=True, ncol=5)
        # ax.legend(legendEntries, labels, loc = 2, ncol = 2, borderaxespad = 0.)
        if saveImageToPath != '':
            fig.savefig(saveImageToPath, dpi=150, format=format)

        return fig

    @staticmethod
    def calcStatsForlist(values, withGeoMean=False, absValues=False):
        out = {'min': np.min(values),
               'max': np.max(values),
               'n': len(values)
               }
        if withGeoMean:
            if absValues:
                tmpValues = [myMaths.ratioAbs(v) for v in values]
            else:
                tmpValues = values
            out['geoMean'] = myMaths.geoMean(tmpValues)
            out['geoSd'] = myMaths.geoSd(tmpValues)
        else:
            if absValues:
                tmpValues = [abs(v) for v in values]
            else:
                tmpValues = values
            out['mean'] = np.mean(tmpValues)
            out['sd'] = np.std(tmpValues, ddof = 1)
        return out

    def calcStatistics(self):
        """
        calc self.scores and calc ranks for score names in self.calcRanksFor
        :return:
        """

        if len(self.compIDs) == 0:
            print >> sys.stderr, 'Cannot calculate statistics if no ScoreComparisons are added !'
            return
        # self.archive['chrLengths'][spOrCombos][myComparison.ID][NumberOfRepetition]
        # myComparison.ID = scoreSim.ID + '-' +  scoreReal.ID
        for (scoreName, scoreCompPerSpOrCombo) in self.archive.iteritems():
            # DEBUG
            # print >> sys.stderr, '%s !!!' % scoreName
            self.scores[scoreName] = {}
            assert 'Diff' in scoreName
            withGeoMean = True if 'Ratio' in scoreName else False
            allMeans = []
            allGeoMeans = []
            for (spOrCombo, scoreRepetsByCompID) in scoreCompPerSpOrCombo.iteritems():
                self.scores[scoreName][spOrCombo] = {}
                allVals = collections.defaultdict(list)
                # FIXME, compute score by sim (for all repets) iSim-iRepet-real
                for (scoreCompID, listValuesRepets) in scoreRepetsByCompID.iteritems():
                    # self.calcStatsForlist(vals) = {'min': X, 'max': X, 'mean': X, 'sd': X, 'n': X}
                    # scoreCompID = iSim-iRepet-real = scoreRepetID
                    # assert len(listValuesRepets) == 1
                    # print >> sys.stderr,  scoreName, spOrCombo, scoreCompID,  listValuesRepets[0]
                    # raw_input()
                    # if withGeoMean and listValuesRepets[0] > 5:
                    #     print >> sys.stderr,  scoreName, spOrCombo, scoreCompID,  listValuesRepets[0]
                    #     raw_input()
                    self.scores[scoreName][spOrCombo][scoreCompID] = listValuesRepets[0]
                    if len(scoreCompID.split('-')) == 3:
                        iSim = scoreCompID.split('-')[-3]
                        # iRepet = scoreCompID.split('-')[-2]
                        # sum (over allSimID) of all simulated values compared to real real
                        # allVals['Sum-all-real'] += listValuesRepets
                        allVals[iSim] += listValuesRepets
                    #else:
                        # DEBUG
                        #print sys.stderr, 'Hey!!!!', scoreCompID
                # self.scores[scoreName][spOrCombo].update(self.calcStatsForlist(allVals['Sum-all-real'], withGeoMean=withGeoMean))
                # self.calcStatsForlist(allVals) = {'min': X, 'max': X, 'mean': X, 'sd': X, 'n': X}

                assert len(set([scoreCompID.split('-')[-3] for scoreCompID in scoreRepetsByCompID.keys() if len(scoreCompID.split('-')) == 3])) == 1
                for iSim in set([scoreCompID.split('-')[-3] for scoreCompID in scoreRepetsByCompID.keys() if len(scoreCompID.split('-')) == 3]):
                    # calc over all repets
                    stats = self.calcStatsForlist(allVals[iSim], withGeoMean=withGeoMean)
                    # self.scores[scoreName][spOrCombo][str(iSim) + '-' + 'real'] = stats
                    self.scores[scoreName][spOrCombo].update(stats)
                    if not withGeoMean:
                        allMeans.append(stats['mean'])
                        print >> sys.stderr, '%s(%s) = %s' %  (scoreName, spOrCombo, stats['mean'])
                    else:
                        allGeoMeans.append(stats['geoMean'])
                        print >> sys.stderr, '%s(%s) = %s' %  (scoreName, spOrCombo, stats['geoMean'])


            # calc stats over all spOrCombo
            if not withGeoMean:
                # no interesting meaning
                #print >> sys.stderr, 'allMeans = %s' % allMeans
                self.scores[scoreName]['mean'] = self.calcStatsForlist(allMeans, withGeoMean=withGeoMean, absValues=True)
                #print >> sys.stderr, 'self.scores[%s][\'mean\'] = %s' % (scoreName, self.scores[scoreName]['mean']['mean'])
                # DEBUG
                print >> sys.stderr, '%s(%s) = %s' %  (scoreName, 'mean', self.scores[scoreName]['mean'])
            else:
                # this is meaning full
                #print >> sys.stderr, 'allGeoMeans = %s' % allGeoMeans
                self.scores[scoreName]['geoMean'] = self.calcStatsForlist(allGeoMeans, withGeoMean=withGeoMean, absValues=True)
                #print >> sys.stderr, 'self.scores[%s][\'geoMean\'] = %s' % (scoreName, self.scores[scoreName]['geoMean']['geoMean'])
                # DEBUG
                print >> sys.stderr, '%s(%s) = %s' %  (scoreName, 'geoMean', self.scores[scoreName]['geoMean'])

        self.scores['mean'] = self.calcStatsForlist([stats['mean']['mean'] for (scoreName, stats) in self.scores.iteritems() if 'Ratio' not in scoreName], withGeoMean=False, absValues=True)
        self.scores['geoMean'] = self.calcStatsForlist([stats['geoMean']['geoMean'] for (scoreName, stats) in self.scores.iteritems() if 'Ratio' in scoreName], withGeoMean=True, absValues=True)
        #self.calcRankList(scoreIDsSelection=self.calcRanksFor)
        # DEBUG
        print >> sys.stderr, '%s(%s) = %s' %  ('general diff', 'mean', self.scores['mean'])
        print >> sys.stderr, '%s(%s) = %s' %  ('general diff', 'geoMean', self.scores['geoMean'])


    def calcNumberOfDataPointsPerComp(self):
        return len([val3 for val1 in self.archive.values() for val2 in val1.values() for val3 in val2.values()])

    def addComp(self, scoreComp, namePrefix=''):
        """
        :param scoreComp:
        :param namePrefix:
        load self.archive[scoreName][spOrCombo][scoreComp.ID].append(val)
        """
        assert isinstance(scoreComp, myScoreComparison)

        addingFirstComp = False
        # Check if new comparison was over same species
        if len(self.namesSpecies) == 0:
            print sys.stderr, 'adding the first comparison'
            addingFirstComp = True
            self.namesSpecies = sorted(scoreComp.nChr.keys())
            self.namesBranches = sorted(scoreComp.branches.keys())
            self.namesCombos = sorted(scoreComp.dist.keys())

        elif self.namesSpecies != sorted(scoreComp.nChr.keys()) or \
            self.namesBranches != sorted(scoreComp.branches.keys()) or \
            self.namesCombos != sorted(scoreComp.dist.keys()):
            raise ValueError('New comparison was over other species than first comparison, yet in an archive all comparisons must be over same species !')

        # scoreComp.ID = scoreSim.ID + '-' + scoreReal.ID
        # !!!!!!
        # 0-1-real : iRepet-iSim-real
        # !!!!!!
        scoreComp.ID = namePrefix + scoreComp.ID
        newCompID = scoreComp.ID
        if newCompID not in self.compIDs:
            self.compIDs.append(newCompID)

        for scoreName in self.scoreNames:
            if addingFirstComp:
                self.archive[scoreName] = {}
            splitedScoreName = scoreName.split('-')
            useLog = False
            if splitedScoreName[-1][:3] == 'Log':
                # FIXME useLog ???
                useLog = True
                splitedScoreName[-1] = splitedScoreName[-1][3:]

            paramNameOfCompScore = splitedScoreName[0]
            # access the parameter of the myScoreComparison object 'scoreComp'
            paramOfCompScore = eval('scoreComp.' + paramNameOfCompScore)
            if len(splitedScoreName) == 1:
                raise ValueError('scoreNames have all at least two components separated by a "-"')
            else:
                for (spOrCombo, subScoreDict) in paramOfCompScore.iteritems():
                    if addingFirstComp:
                        self.archive[scoreName][spOrCombo] = defaultdict(list)
                    if len(splitedScoreName) == 2:
                        k1 = splitedScoreName[1]
                        # DEBUG
                        # print >> sys.stderr, scoreName, k1, subScoreDict
                        val = subScoreDict[k1]
                    elif len(splitedScoreName) == 3:
                        k1 = splitedScoreName[1]
                        k2 = splitedScoreName[2]
                        val = subScoreDict[k1][k2]
                    else:
                        raise ValueError('scoreName %s too long ( >3 components) !' % splitedScoreName)
                    if useLog:
                        # natural logarithm of the value
                        val = np.log(val)

                    self.archive[scoreName][spOrCombo][scoreComp.ID].append(val)
                # DEBUG
                #valsToShow = self.archive[scoreName][spOrCombo][scoreComp.ID] if len(self.archive[scoreName][spOrCombo][scoreComp.ID]) > 3 else self.archive[scoreName][spOrCombo][scoreComp.ID]
                #print >> sys.stderr, 'self.archive[%s][%s][%s] = %s' % (scoreName, spOrCombo, scoreComp.ID, valsToShow)

    def save(self, fileName=None):
        if fileName is None:
            fileName = self.ID + '.archive'
        with open(fileName, 'wb') as outFile:
            pickle.dump(self, outFile, -1)

    @staticmethod
    def load(fileName):
        with open(fileName, 'rb') as inFile:
            score = pickle.load(inFile)
        return score

    def getGammaParameters(self, outPath = 'gammaParameters.csv',
                           parameterFileName = 'parametersG_ABC.v83-7',
                           scoreName = 'sbsDistrib_KS-LogStatistic',
                           meanType = 'arithmetic'):
        # Needs access to parameter files
        # Mean type specifies how cumulation should work: arithmetic or geometric
        def writeOut(out, outPath):
            outStr = ''
            for line in out:
                for entry in line:
                    outStr += str(entry) + '\t'
                outStr = outStr[:-1] + '\n'

            with open(outPath, "wb") as text_file:
                text_file.write(outStr)
                def __init__(self, ID = '1'):
                    self.ID = str(ID)

        # At first, create summary statistics csv
        firstLine = ['sim', 'gammaShape', 'gammaScale', 'cumAbsMean', 'cumRelMean', 'HS-MM_mean', 'HS-MM_min', 'HS-MM_max']
        nCompleteLine = len(firstLine)
        out = [firstLine]
        for simName in self.compIDs:
            sim = '/'.join(simName.split('-')[:-2]) + '/'
            para = myScore.readParameterFile(sim+parameterFileName)
            line = [sim, para['chr:invDistMean'], para['chr:invDistShape']]
            if 'Cumulated' in self.rankList[scoreName].keys():
                for tempSim in self.rankList[scoreName]['Cumulated']['Abs']:
                    if simName == tempSim[0]:
                        line.append(tempSim[1])
                        break
                for tempSim in self.rankList[scoreName]['Cumulated']['Rel']:
                    if simName == tempSim[0]:
                        line.append(tempSim[1])
                        break
            else:
                line.append(float(sum([vals[simName]['mean'] for vals in self.scores[scoreName].itervalues()]))/len(self.namesCombos))
                line.append('NA')

            line.append(self.scores[scoreName][('Homo sapiens', 'Mus musculus')][simName]['mean'])
            line.append(self.scores[scoreName][('Homo sapiens', 'Mus musculus')][simName]['min'])
            line.append(self.scores[scoreName][('Homo sapiens', 'Mus musculus')][simName]['max'])
            if len(line) == nCompleteLine:
                out.append(line)

        writeOut(out, outPath)

        # Secondly, create all data csv
        firstLine = ['sim', 'gammaShape', 'gammaScale', 'cumAbs']
        nCompleteLine = len(firstLine)
        out = [firstLine]
        for (simName, vals) in self.calcMean(scoreName, meanType = meanType).iteritems():
            sim = '/'.join(simName.split('-')[:-2]) + '/'
            para = myScore.readParameterFile(sim+parameterFileName)
            line = [sim, para['chr:invDistMean'], para['chr:invDistShape']]
            for val in vals:
                lineCopy = copy.deepcopy(line)
                lineCopy.append(val)

                if len(lineCopy) == nCompleteLine:
                    out.append(lineCopy)

        modifiedOutPath = outPath.split('.')
        modifiedOutPath[0] += '_allData'
        modifiedOutPath = '.'.join(modifiedOutPath)
        writeOut(out, modifiedOutPath)

    def calcMean(self, scoreName, meanType = 'arithmetic'):
        # Average score over all species or combinations
        # meanType can be 'arithmetic' or 'geometric'
        cumulated = {}
        if meanType == 'arithmetic':
            for combo in self.namesCombos:
                for (sim, vals) in self.archive[scoreName][combo].iteritems():
                    try:
                        old = cumulated[sim]
                    except KeyError:
                        old = [0 for i in xrange(len(vals))]
                    cumulated[sim] = [old[i] + vals[i] for i in xrange(len(vals))]

            for (sim, vals) in cumulated.iteritems():
                cumulated[sim] = [val/float(len(self.namesCombos)) for val in vals]

        elif meanType == 'geometric':
            for combo in self.namesCombos:
                for (sim, vals) in self.archive[scoreName][combo].iteritems():
                    try:
                        old = cumulated[sim]
                    except KeyError:
                        old = [1 for i in xrange(len(vals))]
                    cumulated[sim] = [old[i] * vals[i] for i in xrange(len(vals))]

            for (sim, vals) in cumulated.iteritems():
                cumulated[sim] = [np.power(val, 1.0/(len(self.namesCombos))) for val in vals]

        else:
            raise ValueError('Wrong mean type !')
        return cumulated

    def __init__(self, ID='1', withBA=False):
        self.ID = str(ID)

        # data for checking that all compScores are comparable, have the same species
        # scoreComps of the same simulation (repetitions) have the same compID
        self.compIDs = []
        self.namesSpecies = []
        self.namesBranches = []
        self.namesCombos = []
        self.withBA = withBA
        # Archives for all the comparisons loaded scoreComps they are dictionaries
        # key = ID of the scoreComparison, e.g. 0-real, and the value is a list of all scoreComps of the repetitions
        # (comparisons with the same ID)

        # most basic scoreNames for myOptimizer
        self.scoreNames = { # gene events
                            #'geneEventRatesOnBranches-geneDup-LogRatio', 'geneEventRatesOnBranches-geneLoss-LogRatio',
                            # rearrangements
                            #   fissions and fusions
                            'nChr-Diff', 'nChr-DiffRatio',
                            #   inversions and translocations on branches
                            'branches-nInv-Diff', 'branches-nTransl-Diff', 'branches-nInv-DiffRatio', 'branches-nTransl-DiffRatio',
                            #   inversions and translocations in pairwise comparisons
                            'dist-nInv-Diff', 'dist-nTransl-Diff', 'dist-nInv-DiffRatio', 'dist-nTransl-DiffRatio',
                            # others
                            'dist-nSB-Diff', 'dist-nSB-DiffRatio',
                            'ks_chrLengths-Diff', 'ks_chrLengths-DiffRatio',
                            'ks_sbsLengths-Diff', 'ks_sbsLengths-DiffRatio'
                            }
        # this is useless. Comparisons are between sim and real. And BA is used to compare intra-sim scores
        # if self.withBA:
        #     self.scoreNames |= {# inversions and translocations
        #                         'dist-nInvBA-Diff', 'dist-nInvBA-DiffRatio', 'dist-nTranslBA-Diff', 'dist-nTanslBA-DiffRatio'
        #                         # inversions and translocations on branches
        #                         'branches-nInvBA-Diff', 'branches-nInvBA-DiffRatio', 'branches-nTranslBA-Diff', 'branches-nTanslBA-DiffRatio'
        #                         # others
        #                         'dist-nSBBA-Diff', 'dist-nSBBA-DiffRatio'
        #                         }
        #self.scoreNames = {'nChr-Diff', 'dist-nSB-Diff'} # (ks_D, ks_pVal),   , 'dist-nSBBA-Diff'

        # BA: breakpointAnalyser in magSimus1
        # if withBA:
        #     self.scoreNames = {'chrLengths-AbsoluteMeanFitting', 'chrLengths-MeanFitting', 'chrLengths-RelativeSteepness', 'chrLengths_KS-LogStatistic', 'chrLengths_KS-pValue',
        #                       'sbsDistrib-AbsoluteMeanFitting', 'sbsDistrib-MeanFitting', 'sbsDistrib-RelativeSteepness', 'sbsDistrib_KS-LogStatistic', 'sbsDistrib_KS-LogStatisticBlock', 'sbsDistrib_KS-LogpValue',
        #                       'gMeasure-LogRatio',
        #                       'nChr-Diff',
        #                       'dist-nSB-Diff', 'dist-nSB-DiffRatio', 'dist-nTransl-Diff', 'dist-nTransl-DiffRatio', 'dist-nInv-Diff', 'dist-nInv-DiffRatio',
        #                       'branches-nTransl-Diff', 'branches-nTransl-DiffRatio', 'branches-nInv-Diff', 'branches-nInv-DiffRatio',
        #                       'error-M2006_Inv-Sim', 'error-M2006_Inv-LSE', 'error-M2006_Inv-LSE+M2006', 'error-M2006_Inv-LSE+M2006+PhylDiag',
        #                       'error-M2006_Transl-Sim', 'error-M2006_Transl-LSE', 'error-M2006_Transl-LSE+M2006', 'error-M2006_Transl-LSE+M2006+PhylDiag',
        #                       'geneEventRatesOnBranches-geneBirth-LogRatio', 'geneEventRatesOnBranches-geneDup-LogRatio', 'geneEventRatesOnBranches-geneLoss-LogRatio'#, 'geneEventsToMRCA-geneBirth-LogRatio', 'geneEventsToMRCA-geneDup-LogRatio', 'geneEventsToMRCA-geneLoss-LogRatio', 'geneEventsToMRCA-geneTandemDup-LogRatio', 'geneEventsToMRCA-geneTandemDupProp-LogRatio',
        #                       }
        # else:
        #     self.scoreNames = {'chrLengths-AbsoluteMeanFitting', 'chrLengths-MeanFitting', 'chrLengths-RelativeSteepness', 'chrLengths_KS-LogStatistic', 'chrLengths_KS-pValue',
        #                       'sbsDistrib-AbsoluteMeanFitting', 'sbsDistrib-MeanFitting', 'sbsDistrib-RelativeSteepness', 'sbsDistrib_KS-LogStatistic', 'sbsDistrib_KS-LogStatisticBlock', 'sbsDistrib_KS-LogpValue',
        #                       'gMeasure-LogRatio',
        #                       'nChr-Diff',
        #                       'dist-nSB-Diff', 'dist-nSB-DiffRatio', 'dist-nTransl-Diff', 'dist-nTransl-DiffRatio', 'dist-nInv-Diff', 'dist-nInv-DiffRatio',
        #                       'branches-nTransl-Diff', 'branches-nTransl-DiffRatio', 'branches-nInv-Diff', 'branches-nInv-DiffRatio',
        #                       'error-M2006_Inv-Sim', 'error-M2006_Inv-LSE+M2006+PhylDiag',
        #                       'error-M2006_Transl-Sim', 'error-M2006_Transl-LSE+M2006+PhylDiag',
        #                       'geneEventRatesOnBranches-geneBirth-LogRatio', 'geneEventRatesOnBranches-geneDup-LogRatio', 'geneEventRatesOnBranches-geneLoss-LogRatio'#, 'geneEventsToMRCA-geneBirth-LogRatio', 'geneEventsToMRCA-geneDup-LogRatio', 'geneEventsToMRCA-geneLoss-LogRatio', 'geneEventsToMRCA-geneTandemDup-LogRatio', 'geneEventsToMRCA-geneTandemDupProp-LogRatio',
        #                       }

        if not all([scoreName.replace('-', '').replace('_', '').replace('+', '').isalnum() for scoreName in self.scoreNames]):
            raise ValueError('ScoreNames are only allowed to contain letters, numbers, _, + and - !')

        # in this list are the score names for which ranks are computed
        self.calcRanksFor = []
        # self.calcRanksFor = ['chrLengths-AbsoluteMeanFitting', 'chrLengths_KS-LogStatistic',
        #                      'sbsDistrib-AbsoluteMeanFitting', 'sbsDistrib_KS-LogStatistic', 'sbsDistrib_KS-LogStatisticBlock',
        #                      'gMeasure-LogRatio',
        #                      'nChr-Diff',
        #                      'dist-nSB-Diff', 'dist-nInv-Diff', 'dist-nTransl-Diff',
        #                      'branches-nInv-Diff', 'branches-nTransl-Diff',
        #                      'geneEventRatesOnBranches-geneBirth-LogRatio', 'geneEventRatesOnBranches-geneDup-LogRatio',
        #                      'geneEventRatesOnBranches-geneLoss-LogRatio']

        self.archive = {}
        self.scores = {}

# DEBUG
# lref = [1] * 150 + [2] * 50 + [3] * 40 + [4] *32
# l = [1] * 140 + [2] * 40 + [3] * 30 + [4] *24

def compareListsOfLengths(L, Lref):
    """
    Either compare a list of chromosome lengths or a list of sbs lengths
    :return:
    """
    # FIXME
    cumulated = True
    normalizedIn = None
    #print >> sys.stderr, type(Lref[0])
    #print >> sys.stderr, type(L[0])
    refHist = collections.Counter(Lref)
    hist = collections.Counter(L)
    cumHist = 0
    cumRefHist = 0
    maxDiff = 0
    maxDiffRatio = 1

    for length in sorted(set(refHist.keys() + hist.keys())):
        cumHist += hist[length] if length in hist else 0
        cumRefHist += refHist[length] if length in refHist else 0

        Diff = cumHist - cumRefHist
        if abs(Diff) > abs(maxDiff):
            maxDiff = Diff
            xMaxDiff = length

        num = float(cumHist) if cumHist > 0 else 1.0
        den = float(cumRefHist) if cumRefHist > 0 else 1.0
        ratioDiff =  num / den
        if myMaths.ratioAbs(ratioDiff) > myMaths.ratioAbs(maxDiffRatio):
            maxDiffRatio = ratioDiff
            xMaxRatio = length

    # print >> sys.stderr, 'xMaxDiff=', xMaxDiff
    # print >> sys.stderr, 'xMaxRatio=', xMaxRatio
    assert maxDiffRatio is not None
    #maxDiffs.append((max([np.abs(refHist[idx][1]-hist[idx][1]) for idx in xrange(min(len(refHist), len(hist)))]), iHist))

    # nbBins = len(Lref)
    # #DEBUG
    # # plt.figure()
    # # plt.plot(range(len(Lref)), Lref, color='k', marker='o')
    # # plt.plot(range(len(L)), L, color='b', marker='o')
    # # plt.show()
    # # ranked values normed: RVN
    # (X, rvn_ref) = calcRankedValuesEmpFuncFromList(Lref, nbBins=nbBins, normed=True, normedXaxis=True, increasingOrder=False)
    # (X2, rvn) = calcRankedValuesEmpFuncFromList(L, nbBins=nbBins, normed=True, normedXaxis=True, increasingOrder=False)
    # assert len(X2) == len(X)
    # assert np.allclose(X2, X), "%s == %s" (X2, X)
    # # DEBUG
    # # plt.figure()
    # # plt.plot(new_X, pdf_ref, color='k', marker='o')
    # # plt.plot(new_X, pdf, color='b', marker='o')
    # # plt.show()
    # (diff, flag) = self.calcDiffOfTwoFunctions(zip(X, rvn), zip(X, rvn_ref))
    # normedL = np.array(L) / float(sum(L))
    # normedLref = np.array(Lref) / float(sum(Lref))

    #(ks_D, ks_pVal) = scs.ks_2samp(normedL, normedLref)
    #return ((ks_D, ks_pVal), diff)
    #print >> sys.stderr, maxDiff, maxDiffRatio
    return {'Diff':maxDiff, 'DiffRatio':maxDiffRatio}

class myScoreComparison:

    @staticmethod
    def calcDiffOfTwoFunctions(empFunc, empFuncRef):
        # 'empFunc' stands for distribution
        # The difference between empFunc and empFuncRef is calculated at all x values
        # (empFunc[x] - empFuncRef[x])
        # In case of different x values, 'empFunc' gets interpolated to fit to the
        # same x values as empFuncRef. In this case, a flag is added, as this
        # might lead to distorted differences
        # Both distributions must be sorted by x value
        # OUTPUT: returns tuple of (differences of distributions, flag)

        flag = ''
        if len(empFunc) != len(empFuncRef) or any([empFunc[i][0] != empFuncRef[i][0] for i in xrange(len(empFunc))]):
            # If empFunc have different lengths or different x
            flag = 'different x values'
            # As distributions can't be directly compared, empFunc gets interpolated and extrapolated
            newX = [x for (x,_) in empFuncRef]
            empFunc = zip(newX, extrapolate(empFunc, newX))

        # Get difference between distributions
        # Due to float number magic, it is wisely to round at some point and to add a 0, to avoid having -0.0
        diff = [(empFuncRef[i][0], round(empFunc[i][1] - empFuncRef[i][1], 6) + 0) for i in xrange(len(empFuncRef))]

        return (diff, flag)

    @staticmethod
    def calcProbabilityDensityFromSBLengths(SBLengths, propToGenes=True):
        # SBLength stands for synteny block lengths
        # Takes as input a dict of observed synBlock lengths and
        # returns an orderedTupleList [(synBlockLength: probabilityForGeneToBeInDiagonalOfThatLength), ...]
        # if propToGenes == True: a SB of length 100 counts 100 times; if False, it counts only once
        if propToGenes:
            synBlock_PDF = [(synBlockLength, synBlockLength * count) for (synBlockLength, count) in SBLengths.iteritems()]
        else:
            synBlock_PDF = [(synBlockLength, count) for (synBlockLength, count) in SBLengths.iteritems()]
        synBlock_sum = sum([genes for (synBlockLength, genes) in synBlock_PDF])
        synBlock_PDF = [(synBlockLength, float(genes) / synBlock_sum) for (synBlockLength, genes) in synBlock_PDF]
        synBlock_PDF.sort()
        return synBlock_PDF

    @staticmethod
    def fillHistogrammGapsWithNumber(histogramm, number=0):
        """
        :param histogramm: [..., (length, count), ...]
        :param number: the number that fills the gaps in the distribution
        """
        assert isinstance(histogramm, list)
        maxLength = max([x[0] for x in histogramm])
        index = 0
        filledHistogramm = []
        for length in xrange(1, maxLength+1):
            if length < histogramm[index][0]:
                filledHistogramm.append((length, number))
            elif length == histogramm[index][0]:
                filledHistogramm.append((length, histogramm[index][1]))
                index += 1
            else:
                raise ValueError("You shouldn't be landing here")
        return filledHistogramm

    @staticmethod
    def compareCDFs(distrib, distribRef):
        # 'distrib' stands for distribution
        # Both distributions should be an ordered (by x, unique x) list with [(x, y)]
        # It calculates several statistics for the comparison, namely the difference
        # between both for each x (distrib[x] - distribRef[x]), intersection locations,
        # slope of difference distribution, the mean fitting (mean over differences),
        # and the absolute mean fitting (AMF).
        # OUTPUT: Comparison statistics as a dictionary
        out = {'Flag': []}

        (diff, flag) = myScoreComparison.calcDiffOfTwoFunctions(distrib, distribRef)
        out['Diff'] = diff
        out['Flag'].append(flag)

        # Get intersection locations
        out['Intersections'] = []
        lastSign = 0
        for (i, (x, y)) in enumerate(out['Diff']):
            if i > 0:
                if y != 0 and math.copysign(1, y) != lastSign:
                    if lastSign != 0:
                        # check if there was a sign change
                        out['Intersections'].append(i)
                    lastSign = math.copysign(1, y)

        if len(out['Intersections']) > 1:
            out['Flag'].append('multiple intersections')
        x = [entry[0] for entry in out['Diff']]
        y = [entry[1] for entry in out['Diff']]
        slope, intercept, r_value, pValue, std_err = scs.linregress(x, y)
        out['RelativeSteepness'] = slope
        out['RelativeSteepness_significance'] = pValue
        out['AbsoluteMeanFitting'] = float(sum([abs(recy) for recy in y[1:(len(y)-1)]]))/len(y)
        out['MeanFitting'] = float(sum(y)) / len(y)
        return out

    @staticmethod
    def calcNormalisedCDF_FromOrderedTupleList(orderedTupleList):
        # INPUT: list with [(smallest_x, y1), (second_smallest_x, y2), ...]
        # OUTPUT: cumulated density also as orderedTupleList
        out = [(0, 0)]
        totalSum = sum([y for (x, y) in orderedTupleList])
        cumSum = 0
        for (x, y) in orderedTupleList:
            y = float(y) / totalSum
            cumSum += y
            out.append((x, cumSum))
        return (out)

    @staticmethod
    def compareDistsOrBranchesEstOnSimAndEstOnReal(simEst, realEst):
        # 'simEst' stands for estimations based on a simulation,
        # 'realEst' stands for estimations based on real genomes
        # Compares two dictionaries, usually with {species: float}
        # This function calculates the Diff-erence between the float in simEst and in realEst

        res = {}
        for spOrCombo in realEst:
            #res[spOrCombo] = {}
            if spOrCombo not in simEst:
                continue
            res[spOrCombo] = diffAndRatio(simEst[spOrCombo], realEst[spOrCombo])
            #res[spOrCombo]['Diff'] = diffAndratio['Diff']
            #res[spOrCombo]['DiffRatio'] = diffAndratio['Ratio']
            # FIXME: I do not understand
            # reliability was 'diffEst'
            # Furthermore, a DiffRatio is also calculated, using the provided reliability, a dict which contains for each species a ratio
            # This ratio usually is: Value observed from simulated genome / Simulated Value (e.g. 0.6)
            # res[spOrCombo]['Diff'] = simEst[spOrCombo] - realEst[spOrCombo]
            # if spOrCombo in reliability and reliability[spOrCombo]['Ratio'] != 0 and reliability[spOrCombo]['Ratio'] != -1:
            #     # Case if reliabilityRatio was well calculated
            #     res[spOrCombo]['DiffRatio'] = float(res[spOrCombo]['Diff']) / reliability[spOrCombo]['Ratio']
            # else:
            #     # Case if, e.g., the simulatedValue = 0.5 and estimatedValue = 0 -> reliabilityRatio = 0
            #     # Or if the simulatedValue = 0 and estimatedValue = 0.5 -> reliabilityRatio = -1
            #     res[spOrCombo]['DiffRatio'] = -realEst[spOrCombo]
        return res

    def save(self, fileName = None):
        if fileName is None:
            fileName = self.ID + '.comp'

        with open(fileName, 'wb') as outFile:
            pickle.dump(self, outFile, -1)

    @staticmethod
    def load(fileName):
        with open(fileName, 'rb') as inFile:
            score = pickle.load(inFile)
        return score

    def calcChrLengthsKS(self, simChrLengths, scoreRealChrLengths, scoreRealChrLengths_CDF, chrLengthsDifferences):
        # Calculate Kolmogorov-Smirnoff 2-sided test from chrLengths
        obsReal = []
        for (chr, nGenes) in enumerate(scoreRealChrLengths):
            obsReal = obsReal + list((chr + 1) * np.ones(nGenes, dtype=np.int))

        interpolatedSimCDF = [chrLengthsDifferences[i][1] + scoreRealChrLengths_CDF[i][1] for i in xrange(len(scoreRealChrLengths_CDF))]
        interpolatedSimPDF = [interpolatedSimCDF[i] - interpolatedSimCDF[i - 1] for i in xrange(1, len(interpolatedSimCDF))]
        interpolatedGeneNumbers = [int(round(sum(simChrLengths) * interpolatedSimPDF[i], 0)) for i in xrange(len(interpolatedSimPDF))]
        obs = []
        for (chr, nGenes) in enumerate(interpolatedGeneNumbers):
            obs = obs + list((chr + 1) * np.ones(nGenes, dtype=np.int))

        chrLengthsKS = scs.ks_2samp(obs, obsReal)
        out = {'Statistic': chrLengthsKS[0], 'pValue': chrLengthsKS[1]}
        return out

    def __init__(self, scoreSim, scoreReal):
        self.ID = str(scoreSim.ID) + '-' + str(scoreReal.ID)

        # for each extant species
        self.nChr = {}
        self.geneEventsToMRCA = {}

        # for each comparison
        self.ks_chrLengths = {}
        self.ks_sbsLengths = {}
        self.gMeasure = {}
        self.dist = {}
        # self.sbsDistrib = {}

        # for each branch
        self.geneEventRatesOnBranches = {}
        self.branches = {}

        # FIXME: meanning ?
        # self.error = scoreSim.error
        # for (branch, values) in scoreReal.branchesEst.iteritems():
        #     for stat in self.error[branch]:
        #         statShort = stat.split('_')[1]
        #         self.error[branch][stat]['real_LSE+M2006+PhylDiag'] = values['n'+statShort]
        #         self.error[branch][stat]['real_LSE+M2006Unif+PhylDiag'] = values['n'+statShort+'M2006Unif']

        if len(set(scoreSim.species).difference(set(scoreReal.species))) > 0:
            raise ValueError('myScores are not based over same species!')

        # for all extant species
        for sp in scoreSim.species:
            # compare chromosome numbers
            self.nChr[sp] = diffAndRatio(scoreSim.nChr[sp], scoreReal.nChr[sp])

            # compare chrLengths distributions (each chromosome has the same weight)
            # DEBUG
            # print >> sys.stderr, "type(scoreSim) =", type(scoreSim)
            # print >> sys.stderr, "dir(scoreSim) =", dir(scoreSim)
            # print >> sys.stderr, "type(scoreReal) =", type(scoreReal)
            # print >> sys.stderr, "dir(scoreReal) =", dir(scoreReal)
            self.ks_chrLengths[sp] = compareListsOfLengths(scoreSim.chrLengths[sp],
                                                                scoreReal.chrLengths[sp])

            # compare gene events from MRCA to extant species
            # self.tandemDup[spec] = (scoreSim.tandemDup[spec][0] - scoreReal.tandemDup[spec][0], scoreSim.tandemDup[spec][1] - scoreReal.tandemDup[spec][1])
            self.geneEventsToMRCA[sp] = {}
            for event in scoreReal.geneEventsToMRCA[sp].keys():
                self.geneEventsToMRCA[sp][event] = diffAndRatio(scoreSim.geneEventsToMRCA[sp][event], scoreReal.geneEventsToMRCA[sp][event])

        # for all combinations (combos) of extant species
        for combo in scoreSim.combos:
            # compare all distances
            # scoreX.scoreSim[combo]['nSB'] gives the number of synteny blocks between two extant sepcies
            self.dist[combo] = myScoreComparison.compareDistsOrBranchesEstOnSimAndEstOnReal(scoreSim.distEst[combo],
                                                                                            scoreReal.distEst[combo])

            # compare gMeasures
            self.gMeasure[combo] = diffAndRatio(scoreSim.gMeasure[combo]['Statistic'], scoreReal.gMeasure[combo]['Statistic'])

            # compare synteny blocks distribution
            # counting each sb once
            self.ks_sbsLengths[combo] = compareListsOfLengths(scoreSim.sbsLengths[combo],
                                                                   scoreReal.sbsLengths[combo])
            # the weight of a sb is proportional to the sb length
            # FIXME: this may take some time...
            # sbsLengthsSimWeighted = [lenSb for lenSb in scoreSim.sbsLengths[sp] for _ in range(lenSb)]
            # sbsLengthsRealWeighted = [lenSb for lenSb in scoreReal.sbsLengths[sp] for _ in range(lenSb)]
            # (self.ks_sbsLengths[sp], diff_sbs) = self.compareListsOfLengths(scoreSim.sbsLengths[sp], scoreReal.sbsLengths[sp])

            # Distribution comparison between syntenyBlock distributions
            # First, create PDFs out of synBlock-lengths lists
            sbLengthsCounterSim = Counter(scoreSim.sbsLengths[combo])
            sbLengthsCounterReal = Counter(scoreReal.sbsLengths[combo])
            # TODO, understand this
            # sbsLengthsSim_PDF = myScoreComparison.calcProbabilityDensityFromSBLengths(sbLengthsCounterSim)
            # sbsLengthsReal_PDF = myScoreComparison.calcProbabilityDensityFromSBLengths(sbLengthsCounterReal)
            # # Secondly, get CDFs from PDFs (take care of x-values!)
            # sbsLengthsSim_CDF = myScoreComparison.calcNormalisedCDF_FromOrderedTupleList(sbsLengthsSim_PDF)
            # sbsLengthsReal_CDF = myScoreComparison.calcNormalisedCDF_FromOrderedTupleList(sbsLengthsReal_PDF)
            # # Finally, compare both distributions
            # self.sbsDistrib[combo] = myScoreComparison.compareCDFs(sbsLengthsSim_CDF, sbsLengthsReal_CDF)

        # for all branches
        for branch in scoreSim.branchesEst:
            self.branches[branch] = myScoreComparison.compareDistsOrBranchesEstOnSimAndEstOnReal(scoreSim.branchesEst[branch],
                                                                                                 scoreReal.branchesEst[branch])
            self.geneEventRatesOnBranches[branch] = {}
            for event in scoreReal.geneEventRatesOnBranches[branch].keys():
                self.geneEventRatesOnBranches[branch][event] = {
                    'Ratio': myScore.calcRatio(scoreSim.geneEventRatesOnBranches[branch][event],
                                               scoreReal.geneEventRatesOnBranches[branch][event])
                }

class myScore:
    # Calculates different statistics of simulated (or real) genomes and can compare these statistics objects

    # Calculates chrom x chrom number of homologues genes, ignoring unclear cases
    def calcHomologies(self, g1, g2, ancGenes, minChromLength=0):
        assert isinstance(g1, myLightGenomes.LightGenome) and isinstance(g2, myLightGenomes.LightGenome)
        # step 1 :filter genomes and rewrite in tandem blocks if needed
        ##############################################################
        # rewrite genomes by family names (ie ancGene names)
        g1_aID = myMapping.labelWithFamID(g1, ancGenes)
        g2_aID = myMapping.labelWithFamID(g2, ancGenes)
        # genes that are not in ancGene have a aID=None
        ((g1_aID, mGf2Go1, (nCL1, nGL1)), (g2_aID, mGf2Go2, (nCL2, nGL2))) = \
            myDiags.filter2D(g1_aID, g2_aID, myDiags.FilterType[1], minChromLength)

        def createGeneChrDic(g1_aID):

            g1ong2 = {}
            for k, v in g1_aID.iteritems():
                genes = Counter([gt[0] for gt in v])
                for gene in genes:
                    if gene in g1ong2:
                        g1ong2[gene] += [k] * genes[gene]
                    else:
                        g1ong2[gene] = [k] * genes[gene]
            return (g1ong2)

        g1ong2 = createGeneChrDic(g1_aID)
        g2ong1 = createGeneChrDic(g2_aID)
        matNoDup = {}
        for g1chr in g1_aID:
            matNoDup[g1chr] = {}
            for g2chr in g2_aID:
                # Pseudo-count
                matNoDup[g1chr][g2chr] = 1

        # All homologies which are in more than 1 chromosome are ignored
        for gene in g1ong2:
            c1 = set(g1ong2[gene])
            c2 = set(g2ong1[gene])
            if len(c1) > 1 or len(c2) > 1:
                continue
            matNoDup[c1.pop()][c2.pop()] += 1

        return matNoDup

    @staticmethod
    def calcGMeasureForOneCountMatrix(homolMat):
        # Returns a dict of 3: the global G-measure, log pValue and the homologyMatrix
        spec1sum = {key: sum(cont.values()) for (key, cont) in homolMat.iteritems()}
        spec2sum = {}
        for cont1 in homolMat.values():
            for (key2, cont2) in cont1.items():
                if key2 in spec2sum.keys():
                    spec2sum[key2] += cont2
                else:
                    spec2sum[key2] = cont2

        n = sum(spec1sum.values())
        gTotal = 0
        # gMat is the chromosome by chromosome gMeasure
        # gMat = {key1: {key2: 0 for key2 in cont1.keys()} for (key1, cont1) in homolMat.items()}
        for (key1, cont1) in homolMat.items():
            for (key2, cont2) in cont1.items():
                expected = spec1sum[key1] * spec2sum[key2] * 1.0 / n
                gRec = 2 * homolMat[key1][key2] * (
                    np.log(homolMat[key1][key2]) - np.log(expected))
                # gMat[key1][key2] = gRec
                gTotal += gRec
        degreesOfFreedom = (len(homolMat.keys())-1) * (len(homolMat.itervalues().next())-1)
        pValue_log = scs.chi2.logpdf(gTotal, degreesOfFreedom)
        return (gTotal, pValue_log)

    # Deprecated as  not used any more (use myMapping.calcDupDist)
    @staticmethod
    def calcTandemDup(genome, family, allowedGap=0, calcCDF = True):
        genome_fID = myMapping.labelWithFamID(genome, family)

        (genomeFilt, gf2gfID, _) = myMapping.remapFilterGeneContent(genome_fID, {None})

        familySizes = defaultdict(int)
        for chrom in genomeFilt:
            for gene in genomeFilt[chrom]:
                familySizes[gene.n] += 1

        nGeneDupl = 0
        for familySize in familySizes.values():
            nGeneDupl += familySize - 1

        dupsInDist = 0
        for genes in genomeFilt.itervalues():
            for iGene in xrange(len(genes) - 1 - allowedGap):
                correctedDist = [genes[iGene].n != genes[iGene + gap].n for gap in xrange(1, allowedGap)]
                if all(correctedDist) and genes[iGene].n == genes[iGene + allowedGap + 1].n:
                    dupsInDist += 1

        if calcCDF:
            for gap in xrange(allowedGap):
                for iGene in xrange(len(genes) - 1 - gap):
                    correctedDist = [genes[iGene].n != genes[iGene + gap].n for gap in xrange(1, gap + 1)]
                    if all(correctedDist) and genes[iGene].n == genes[iGene + gap + 1].n:
                        dupsInDist += 1

        return (float(dupsInDist), float(nGeneDupl))

    def M2006(self, t, c, d, p, q, ci):
        # ci either being a dict or an int
        if isinstance(ci, dict):
            ci = sum(ci.values())
        sumsum = 0
        for pcont in p:
            for qcont in q:
                sumsum += (1 - qcont) ** (2 * t * pcont)
        out = (c * d - ci) - float(d - 1) / d * sumsum
        return out

    def M2006Unif(self, t, c, d, ci):
        # ci either being a dict or an int
        if isinstance(ci, dict):
            ci = sum(ci.values())
        # cf manuscript
        out = (c * d - ci) - c * float(d - 1) * (float(d-1)/d) ** (float(2 * t) / c)
        return out

    def calcConnectionMatrix(self, speciesTree, combos):
        Xdic = {}
        XCols = self.branches
        XRows = combos
        for line in XRows:
            Xdic[line] = {}
            spec1 = line[0]
            spec2 = line[1]
            link = speciesTree.dicLinks[spec1][spec2]
            mrca = speciesTree.lastCommonAncestor([spec1, spec2])
            for (i, cont) in enumerate(XCols):
                if cont != mrca and cont in link:
                    Xdic[line][cont] = 1
                else:
                    Xdic[line][cont] = 0

        return Xdic

    def calcBranchEstimatesFromDist(self, speciesTree, dist, distName):
        X = np.zeros([len(self.combos), len(speciesTree.parent)])
        XCols = self.branches
        XRows = self.combos
        Y = np.zeros(len(self.combos))
        for (iline, line) in enumerate(XRows):
            spec1 = line[0]
            spec2 = line[1]
            combo = (spec1, spec2)
            link = speciesTree.dicLinks[spec1][spec2]
            mrca = speciesTree.lastCommonAncestor([spec1, spec2])
            for (i, cont) in enumerate(XCols):
                if cont != mrca and cont in link:
                    X[iline, i] = 1
            Y[iline] = dist[combo][distName]

        for col1 in xrange(X.shape[1] - 1):
                    for col2 in xrange(col1 + 1, X.shape[1]):
                        if all(X[:, col1] == X[:, col2]):
                            duplCol = [col1, col2]
                            break

        Xs = np.delete(X, duplCol[1], 1)
        prop = float(speciesTree.parent[XCols[duplCol[1]]].distance) / (
            speciesTree.parent[XCols[duplCol[0]]].distance +
            speciesTree.parent[
                XCols[duplCol[1]]].distance)
        YEst = sco.nnls(Xs, Y)[0]
        YEst = np.insert(YEst, duplCol[1], YEst[duplCol[0]] * prop)
        YEst[duplCol[0]] = YEst[duplCol[0]] * (1 - prop)

        branchesEst = {}
        for (i, key) in enumerate(list(XCols)):
            branchesEst[key] = list(YEst)[i]

        return branchesEst

    def effectiveValuesFromSimulationLogErr(self, pathToLogErr):
        # Read the logErr to determine the realized number of chromosomal rearrangments
        # in the simulation. Fill the class variable self.branchesSim
        nChrAncestors = {}
        branchesSim = {}

        with open(pathToLogErr) as f:
            content = f.readlines()

        start = 0
        notFoundAmniotaNChr = True
        for (i, line) in enumerate(content):
            if line[:9] == '# Branch ':
                if (start > 0):
                    branchesSim[branchName] = {'ind': (start, i)}
                start = i
                branchName = line.split('-> ')[1].split(' (')[0]
            elif notFoundAmniotaNChr and line[:17] == '| startChrLengths':
                nChrAncestors['Amniota'] = len(line.split('|')[2].strip().split(' '))
                notFoundAmniotaNChr = False

        branchesSim[branchName] = {'ind': (start, i)}

        for (branchName, item) in branchesSim.items():
            statsCurrBranch = {}
            for line in content[item['ind'][0]:item['ind'][1]]:
                if line[:13] == 'propTandemDup':
                    propTandemDup = line.split(' ')[3][:-1]
                    if propTandemDup == 'NTR':
                        statsCurrBranch['propTandemDup'] = 0.5
                    else:
                        statsCurrBranch['propTandemDup'] = float(propTandemDup)/100
                elif branchName not in self.nChr and line[:12] == '#chromosomes':
                    nChrAncestors[branchName] = int(line.split('=')[1].strip())
                elif line[:10] == '#chrInvert':
                    statsCurrBranch['nInv'] = int(line.split('=')[1].strip().split(' ')[0])
                elif line[:12] == '#chrTransloc':
                    statsCurrBranch['nTransl'] = int(line.split('=')[1].strip().split(' ')[0])
                elif line[:10] == '#chrFusion':
                    statsCurrBranch['nFusion'] = int(line.split('=')[1].strip().split(' ')[0])
                elif line[:11] == '#chrFission':
                    statsCurrBranch['nFission'] = int(line.split('=')[1].strip().split(' ')[0])
                elif line[:40] == 'number of bra-standard inner breakpoints':
                    statsCurrBranch['BPDetectableWOGE'] = int(line.split('=')[1].strip())
                elif line[:38] == 'number of bra-inner breakpoints reused':
                    BPReuseA = int(line.split('=')[1].strip())
                elif line[:52] == 'number of bra-breakpoints at chromosomal extremities':
                    BPReuseB = int(line.split('=')[1].strip())
                elif line[:57] == 'number of bra-synteny-blocks removed because of gene loss':
                    statsCurrBranch['BPLossByGE'] = int(line.split('=')[1].strip())
                elif line[:44] == 'number of bra-undetectable inner breakpoints':
                    statsCurrBranch['BPRecovered'] = int(line.split('=')[1].strip())

            if 'BPDetectableWOGE' in statsCurrBranch:
                statsCurrBranch['BPReuse'] = BPReuseA + BPReuseB
                statsCurrBranch['BPDetectable'] = statsCurrBranch['BPDetectableWOGE'] - statsCurrBranch['BPLossByGE']
                statsCurrBranch['BPTheoretical'] = statsCurrBranch['BPDetectableWOGE'] + statsCurrBranch['BPReuse'] + statsCurrBranch['BPRecovered']
            branchesSim[branchName] = statsCurrBranch

        return (nChrAncestors, branchesSim)

    def calcDistBetweenExtantSpecies(self, valuesPerBranch, connectionMatrix, speciesTree):
        dist = {}
        for (irow, combo) in enumerate(connectionMatrix):
            # 'nBreakpoints': min(nChr[combo[0]], nChr[combo[1]]) - 1
            dist[combo] = defaultdict(float)
            dist[combo]['nSB'] += self.nChr[speciesTree.lastCommonAncestor(combo)]
            for (icol, branch) in enumerate(connectionMatrix[combo]):
                if connectionMatrix[combo][branch] == 1:
                    dist[combo]['nSB'] += (valuesPerBranch[branch]['nFission'] + 2*valuesPerBranch[branch]['nInv'] +
                                                2*valuesPerBranch[branch]['nTransl'])
                    dist[combo]['nInv'] += valuesPerBranch[branch]['nInv']
                    dist[combo]['nTransl'] += valuesPerBranch[branch]['nTransl']
                    dist[combo]['nFission'] += valuesPerBranch[branch]['nFission']
                    if 'BPDetectableWOGE' in valuesPerBranch[branch]:
                        dist[combo]['BPTheoretical'] += valuesPerBranch[branch]['BPTheoretical']
                        dist[combo]['BPDetectable'] += valuesPerBranch[branch]['BPDetectable']
                        dist[combo]['BPReuse'] += valuesPerBranch[branch]['BPReuse']
                        dist[combo]['BPLossByGE'] += valuesPerBranch[branch]['BPLossByGE']
        return dist

    def calcSyntenyBlocks(self, genome1, genome2, familyMRCA, verbose=False):
        # use phylDiag
        sbsInPairComp = myDiags.extractSbsInPairCompGenomes(genome1, genome2, familyMRCA,
                                                            tandemGapMax=self.phylDiagTandemGapMax,
                                                            gapMax=self.phylDiagGapMax,
                                                            distanceMetric='CD',
                                                            pThreshold=self.phylDiagPValue,
                                                            distinguishMonoGenicDiags=True,
                                                            identifyMicroRearrangements=True,
                                                            identifyMonoGenicInvs=True,
                                                            gapMaxMicroInv=1,
                                                            truncationMax=self.phylDiagTruncationMax,
                                                            filterType=self.phylDiagFilter,
                                                            verbose=verbose)
        sbLengths = [len(sb.la) for (_, sb) in sbsInPairComp.iteritems2d()]
        nSB = len(sbLengths)
        return (sbsInPairComp, nSB, sbLengths)

    def calcTranslInvEstimation(self, nSB, nChr1, nChr2, ci, chrLengths1=None, chrLengths2=None, fM2006=M2006):
        if fM2006 == self.M2006:
            assert (chrLengths1 is not None) and (chrLengths2 is not None)
            p = [float(x) / sum(chrLengths1) for x in chrLengths1]
            q = [float(x) / sum(chrLengths2) for x in chrLengths2]
        c = nChr1
        d = nChr2
        if fM2006 == self.M2006:
            nTransl = sco.newton(fM2006, 20, args=(c, d, p, q, ci), tol=0.01)
        elif fM2006 == self.M2006Unif:
            nTransl = sco.newton(self.M2006Unif, 20, args=(c, d, ci), tol=0.01)
        else:
            raise ValueError('function fM2006 should be either self.M2006, or self.M2006Unif')
        nInv = float(nSB - max(c, d)) / 2 - nTransl
        return (nTransl, nInv)

    @staticmethod
    def calcRatio(a, b, valueIfDivisionByZero = -1):
        if b != 0:
            ratio = float(a) / b
        else:
            if a == 0:
                # Define 0 / 0 == 1
                ratio = 1
            else:
                # Case if division by 0, hence ratio would be infinite -> default to -1
                ratio = valueIfDivisionByZero
        return ratio

    def compareSimulationAndEstimation(self):
        # TODO, do it for M2006 non-uniform
        suffixes = ['', 'BA'] if self.withBA else ['']

        # compare pairwise comparisons of extant species
        for combo in self.distEst.iterkeys():
            if combo not in self.distDiff:
                self.distDiff[combo] = {}
            for statBaseName in ('nInv', 'nTransl', 'nSB'):
                for suffix in suffixes:
                    statEstName = statBaseName
                    statSimName = statBaseName + suffix
                    statDiffName = statSimName
                    self.distDiff[combo][statDiffName] = diffAndRatio(self.distEst[combo][statEstName], self.distSim[combo][statSimName])

        # compare values on each branch
        for combo in self.branchesSim.iterkeys():
            self.branchesDiff[combo] = {}
            for statBaseName in ('nInv', 'nTransl', 'nSB'):
                for suffix in suffixes:
                    statEstName = statBaseName
                    statSimName = statBaseName + suffix
                    statDiffName = statSimName
                    self.branchesDiff[combo][statDiffName] = diffAndRatio(self.branchesEst[combo][statEstName],
                                                                          self.branchesSim[combo][statSimName])

    def calcErrorSplitting(self):
        for branch in self.branches:
            self.error[branch] = {}
            for key in ('M2006_Transl', 'M2006_Inv'):
                toI = key.split('_')[1]
                self.error[branch][key] = OrderedDict()
                self.error[branch][key]['Sim'] = self.branchesSim[branch]['n'+toI]

                if self.withBA:
                    self.error[branch][key]['LSE'] = self.branchesEst[branch]['nLSE'+toI]
                    self.error[branch][key]['LSE+M2006'] = self.branchesEst[branch]['n'+toI+'BAEstM2006']
                    self.error[branch][key]['LSE+M2006Unif'] = self.branchesEst[branch]['n'+toI+'BAEstM2006Unif']

                self.error[branch][key]['LSE+M2006+PhylDiag'] = self.branchesEst[branch]['n'+toI]
                self.error[branch][key]['LSE+M2006Unif+PhylDiag'] = self.branchesEst[branch]['n'+toI+'M2006Unif']

    @staticmethod
    def readParameterFile(inPath):
        userRatesReader = myFile.myTSV.readTabular(inPath, [str, str])

        inDict = {}
        for x in userRatesReader:
            x = list(x)
            x[1] = x[1].replace('\r', '')
            try:
                if x[0] == "startChrLengths":
                    inDict[x[0]] = tuple(float(i) for i in x[1].strip().replace('(', '').replace(')', '').split(','))
                elif (x[0] == "forceClustering" or x[0] == "b_randomAccel" or x[0] == "printIntermGenomes"):
                    if x[1].strip() == "True":
                        inDict[x[0]] = True
                    elif x[1].strip() == "False":
                        inDict[x[0]] = False
                    else:
                        raise ValueError
                elif ((x[0][:3] == 'chr' and x[0][:22] != 'chr:invDist')
                      or x[0][:5] == 'rate:' or x[0][:5] == 'gene:' or x[0] == 'globalFactor'):
                    inDict[x[0]] = float(x[1])
                else:
                    inDict[x[0]] = x[1]
            except ValueError:
                print >> sys.stderr, x
                print >> sys.stderr, 'Wrong parameter for "' + x[
                    0] + '" (e.g. character where should be a number). Default used.'
            except KeyError:
                print >> sys.stderr, 'Parameter "' + x[0] + '" unknown and hence ignored.'

        return inDict

    def writeToCSV(self, outPath, variable = 'error'):
        # Variable can either be 'error', 'branches', 'dist'
        out = ''
        if variable == 'error':
            var = self.error
            for (key, specDict) in var.iteritems():
                out += key+'\n'
                # Head of table:
                out += '\t' + '\t'.join(specDict.itervalues().next().keys()) + '\n'
                for (spec, statDict) in specDict.iteritems():
                    out += str(spec) + '\t' + '\t'.join([str(x) for x in statDict.itervalues()]) + '\n'
                out += '\n'
        else:
            if variable == 'dist':
                var = self.distEst
            elif variable == 'branches':
                var = self.branchesEst
            else:
                raise('Wrong variable name!')
            out = '\t'.join(var.itervalues().next().keys()) + '\n'
            for (spec, statDict) in var.iteritems():
                out += str(spec) + '\t' + '\t'.join([str(x) for x in var[spec].itervalues()]) + '\n'
            out += '\n'

        with open(outPath, "w") as text_file:
            text_file.write(out)

    def save(self, fileName = None):
        if fileName is None:
            fileName = self.ID + '.score'
        with open(fileName, 'wb') as outFile:
            pickle.dump(self, outFile, -1)

    @staticmethod
    def load(fileName):
        with open(fileName, 'rb') as inFile:
            score = pickle.load(inFile)
        return score

    @staticmethod
    def calcGeneEvents(genomeOrFamDesc, family, tandemDupGap = None):
        geneEvents = {}
        geneEvents['geneLoss'] = myMapping.calcNumberOfGeneDeletions(genomeOrFamDesc, family)
        geneEvents['geneBirth'] = myMapping.calcNumberOfGeneBirths(genomeOrFamDesc, family)
        geneEvents['geneDup'] = myMapping.calcNumberOfGeneDuplications(genomeOrFamDesc, family)
        if tandemDupGap is not None:
            (distDict, genomeFilt) = myMapping.calcDupDist(genomeOrFamDesc, family)
            geneEvents['geneTandemDup'] = myMapping.calcTandemDupFromDist(distDict, gapMax=tandemDupGap)
            geneEvents['geneTandemDupProp'] = myScore.calcRatio(geneEvents['geneTandemDup'], geneEvents['geneDup'])
        return geneEvents

    def calcGeneEventRatesOnBranches(self, speciesTree, loadedGenomes, loadedFamilies):
        for branch in speciesTree.indBranches.keys():
            anc = speciesTree.parent[branch].name
            branchLength = speciesTree.parent[branch].distance
            ancFamily = loadedFamilies[anc]
            if branch in loadedGenomes.keys():
                genomeOrFamDesc = loadedGenomes[branch]
            else:
                genomeOrFamDesc = loadedFamilies[branch]
            geneEvents = myScore.calcGeneEvents(genomeOrFamDesc, ancFamily)
            self.geneEventRatesOnBranches[branch] = {eventName: float(eventValue)/branchLength for (eventName, eventValue) in geneEvents.iteritems()}

    def computeBA_Sbs(self, sp1, sp2, genome1, genome2, filePath):
        sbsGenomeBA = myLightGenomes.LightGenome(filePath + 'sbs.genes.' + '.'.join(sp1.split(' '))
                                                  + '.' + '.'.join(sp2.split(' '))+'.list.bz2')
        def geneGenomeToStrGenome(genome1):
            genome1str ={}
            for (chr, genes) in genome1.iteritems():
                genome1str[chr] = [gene.n for gene in genes]
            return(genome1str)
        genome1str = geneGenomeToStrGenome(genome1)
        genome2str = geneGenomeToStrGenome(genome2)

        sbMapDict = defaultdict(list)
        for key, sb in sbsGenomeBA.iteritems():
            for (chr1, genes1) in genome1str.iteritems():
                if sb[0].n in genes1:
                    break
            for (chr2, genes2) in genome2str.iteritems():
                if sb[0].n in genes2:
                    break
            if chr2 not in sbMapDict[chr1]:
                sbMapDict[chr1].append(chr2)
        ciBA = sum([len(x) for x in sbMapDict.itervalues()])
        nSBBA = len(sbsGenomeBA)
        sbLengthsBA = [len(sb) for sb in sbsGenomeBA.values()]
        return (sbsGenomeBA, nSBBA, ciBA, sbLengthsBA)

    def compareEstSbsToRefSbsFromBA(self, sbsGenomeBA, sbsInPairComp, verbose=True):
        phylDiagSbsGenome = myDiags.buildGenomeFromSbs(sbsInPairComp, sbLengthThreshold=None)
        (eff, (sTp, sFp, sFn)) = myIntervals.compareGenomes(phylDiagSbsGenome, sbsGenomeBA,
                                                     mode='chromExtremity', oriented=True, verbose=verbose)
        res = {'TP': eff.tp, 'FP': eff.fp, 'FN': eff.fn, 'TN': eff.tn}
        return res

    def __init__(self, genomeFiles, ancFiles, treePath,
                 logErrPath=None,
                 firstAncestralGenome='Amniota',
                 tandemDupGap=5,
                 phylDiagGapMax=5,
                 phylDiagPValue=None,
                 phylDiagTruncationMax=5,
                 phylDiagFilter=myDiags.FilterType.InBothGenomes,
                 ID='1',
                 verbose=False,
                 verbose2=False):

        startTime = time.clock()
        phylDiagTime = 0
        # Switch off stderr stream
        if not verbose:
            stderrOld = sys.stderr
            stdoutOld = sys.stdout
            nullDevice = open(os.devnull, 'w')
            sys.stderr = nullDevice
            sys.stdout = nullDevice

        #################################################
        self.ID = str(ID)

        #################################################
        # Parameters for phylDiag when calculating observed breakpoints

        self.phylDiagTandemGapMax = tandemDupGap
        self.phylDiagGapMax = phylDiagGapMax
        self.phylDiagPValue = phylDiagPValue
        self.phylDiagTruncationMax = phylDiagTruncationMax
        self.phylDiagFilter = phylDiagFilter

        #################################################
        # Data related variables

        # As an object as class variable is difficult to serialize,
        # the speciesTree itself is not stored as an object Variable
        speciesTree = myPhylTree.PhylogeneticTree(treePath)
        self.species = list(speciesTree.listSpecies)
        self.species.sort()
        self.branches = speciesTree.parent.keys()
        self.combos = []
        self.MRCAOfAll = ancFiles.replace('%s', '{0}').format(firstAncestralGenome)

        #################################################
        # myScore relevant statistics
        # {Species: Number of chromosomes}
        self.nChr = {}
        # {Species: [..., chrLength, ...] by decreasing lengths (in genes)}
        self.chrLengths = {}
        # {SpeciesCombo: {SBLength: count}}
        self.sbsLengths = {}

        # {SpeciesCombo: {distance type (either nBreakpoints, nInv or nTransl):
        #   estimated distance based on PhylDiag and Mazowita2006}}
        self.distEst = {}
        # {BranchName: {distance type (either nInv or nTransl):
        #   estimated distance based on PhylDiag, Mazowita2006 and Non-linear least squares}
        self.branchesEst = {}

        # {geneEvent: {Species: value}}
        # self.tandemDup = {}
        self.geneEventsToMRCA = {}
        self.geneEventRatesOnBranches = {}
        # {SpeciesCombo: (totalGMeasure, log pValue)}
        self.gMeasure = {}

        # Realized estimation error
        # FIXME understand its meaning
        self.error = {}

        # {SpeciesCombo: Input for formula (9), Mazowita 2006}
        self.M2006Input = {}

        # If the breakpoint analyser was used
        filePath = '/'.join(genomeFiles.split('/')[:-1]) + '/'
        if any(['sbs.genes' in fileName for fileName in os.listdir(filePath)]):
            self.withBA = True
        else:
            # in the case of real genomes for instance, or if magSimus was not launched with the breakpoint analyzer to save time
            self.withBA = False

        if self.withBA:
            self.phylDiagEffComparedToBA = {}
            self.nbSbsFromBA = {}
            sbsFiles = filePath + 'sbs.genes.%s.%s.list.bz2'
            # calculate and write reference-sbs of all pairwise comparisons
            myBA.computePairwiseSyntenyBlocksOfMagSimus(speciesTree, genomeFiles, ancFiles, sbsFiles)
            # record efficiencies for each pairwise comparison when the BA is actif
            self.phylDiagEff = {}

        ################################################
        # The next 4 variables stay empty if logErrPath is not given, i.e.
        # if the number of simulated events is not specified

        # {SpeciesCombo: {type of distance (either nBreakpoints, nInv or nTransl):
        #   with MagSimus simulated distance}}
        self.distSim = {}
        # {BranchName: {comparisonFactor (either nInv or nTransl):
        #   with MagSimus simulated distance and spread on branches by Non-linear least squares}}
        self.branchesSim = {}

        # {SpeciesCombo: {type of distance (either nBreakpoints, nInv or nTransl):
        #   {either difference (estimated - simulated) or ratio (estimated / simulated}}}
        # For the ratio, x/0 is defined as 1 for x == 0 and -1 else
        self.distDiff = {}
        # {BranchName: {comparisonFactor (either nInv or nTransl):
        #   {either difference (observed - simulated) or ratio (observed / simulated}}}
        # For the ratio, x/0 is defined as 1 for x == 0 and -1 else
        self.branchesDiff = {}

        ################################################
        # End of class variable and begin of myScore calculation

        # Calculate all species combos
        # sort each combo by alphabetical order, and alls combos by alphabetical order too
        # Thus: [('A', 'B'), ('A', 'C'), ('B', 'C'), ...]
        for (isp1, sp1) in enumerate(self.species):
            for isp2 in xrange(isp1+1, len(self.species)):
                sp2 = self.species[isp2]
                self.combos.append(tuple(sorted((sp1, sp2))))
        self.combos = sorted(self.combos)
        # self.combos = [(sp1, sp2) for (iCombo, (sp1, sp2)) in enumerate(itertools.combinations(self.species, 2))]

        # 1) Read number of realized events in simulation from standard error output
        #    Count the number of chromosomes in simulated ancestral genomes
        if logErrPath is not None:
            connectionMatrix = self.calcConnectionMatrix(speciesTree, self.combos)
            (nChrAncestors, self.branchesSim) = self.effectiveValuesFromSimulationLogErr(logErrPath)
            assert 'Amniota' in nChrAncestors
            self.nChr.update(nChrAncestors)
            # compute self.distSim[combo][statName] for statName in (nInv, nTransl, nSB, ...) for each combo
            self.distSim = self.calcDistBetweenExtantSpecies(self.branchesSim, connectionMatrix, speciesTree)

        # 2) load extant and ancestral genomes (families)
        # 2.1) load ancestral genomes
        loadedFamilies = {}
        for anc in set(speciesTree.allDescendants[firstAncestralGenome]) & set(speciesTree.listAncestr):
            loadedFamilies[anc] = myLightGenomes.Families(ancFiles.replace('%s', '{0}').format(anc.replace(' ', '.')))
        familyOfMRCA = loadedFamilies[firstAncestralGenome]

        ############################
        # single genome statistics
        ############################
        # 2.2) Load extant genomes
        loadedGenomes = {}
        for sp in self.species:
            genome = genomeFiles.replace('%s', '{0}').format(sp.replace(' ', '.'))
            genome = myLightGenomes.LightGenome(genome)
            genome.removeUnofficialChromosomes()
            # "sex chromosomes Y andWwere also not considered. Y is the small, male sex chromosome in mammals, andWis the
            # small, female sex chromosome in birds. Both chromosomes have in common that they do not have many genes on
            # them and have very altered evolutionary constraints compared to the other chromosomes. An important
            # question is if chromosome X and Z, the partners of Y and W, should be included, as they also behave quite
            # differently compared to other chromosomes due to their nature as sex chromosomes. As they nevertheless fit
            # well within the size distribution of autosomes and represent a substantial part of the genome, we decided
            # to follow Mazowita et al. (2006) and to include them in our data." (Lucas Tittman)
            # Delete Y Chromosome
            for chr in ['Y', 'W']:
                genome.pop(chr, None)
            loadedGenomes[sp] = genome

            # ###############################################################################
            # 2.2.1) Chromosome number
            ################################################################################
            # this cancels what we computed above when we compute the score of a simulation
            self.nChr[sp] = len(genome)

            ################################################################################
            # 2.2.2) Chromosome size distribution
            ################################################################################
            if sp not in self.chrLengths:
                self.chrLengths[sp] = [len(chr) for chr in genome.values()]
                # sort by decreasing lengths
                self.chrLengths[sp].sort(reverse=True)

            ################################################################################
            # 2.2.3) Calculate the geneEventsToMRCA (Most recent common ancestor)
            ################################################################################
            self.geneEventsToMRCA[sp] = self.calcGeneEvents(genome, familyOfMRCA, tandemDupGap)

        #################################################################
        # Comparison of two extant genomes
        #################################################################
        # DEBUG
        assert round(myMaths.combinations(len(self.species), 2)) == len(self.combos)
        progressBar = myTools.ProgressBar(len(self.combos))
        for (iCombo, (sp1, sp2)) in enumerate(self.combos):
            genome1 = loadedGenomes[sp1]
            genome2 = loadedGenomes[sp2]
            combo = (sp1, sp2)
            if combo not in self.distEst:
                assert (combo not in self.gMeasure) and (combo not in self.M2006Input)
                self.distEst[combo] = {}
                self.gMeasure[combo] = {}
                self.M2006Input[combo] = {}
            if combo not in self.distSim:
                (sp1, sp2) = combo
                assert (sp2, sp1) not in self.distSim
                # Might be already instanciated because of loading data from logErrPath (cf 1) )
                self.distSim[combo] = {}

            mrca = speciesTree.lastCommonAncestor([sp1, sp2])
            assert mrca in loadedFamilies
            familyMRCA = loadedFamilies[mrca]
            ################################################################################
            # 4) Distribution of syntenyBlockLengths
            ################################################################################
            phylDiagTime -= time.clock()
            (sbsInPairComp, self.distEst[combo]['nSB'], self.sbsLengths[combo]) = self.calcSyntenyBlocks(genome1, genome2, familyMRCA, verbose=verbose2)
            phylDiagTime += time.clock()
            if self.withBA:
                # sbsLengthsBA is not used
                (sbsGenomeBA, self.distSim[combo]['nSBBA'], ciBA, sbsLengthsBA) = self.computeBA_Sbs(sp1, sp2, genome1, genome2, filePath)
                self.phylDiagEff[combo] = self.compareEstSbsToRefSbsFromBA(sbsGenomeBA, sbsInPairComp, verbose=verbose2)
            ################################################################################
            # 5) Calculate G-measure
            ################################################################################
            homolMat = self.calcHomologies(genome1, genome2, familyMRCA)
            (self.gMeasure[combo]['Statistic'], self.gMeasure[combo]['pValue_log']) = \
                myScore.calcGMeasureForOneCountMatrix(homolMat)

            # ###############################################################################
            # 6) Estimate number of translocations and inversions
            ################################################################################
            # TODO : nTranslBAEstM2006 -> nTranslBAM2006 ...
            # rename variables for visibility
            nChr1, chrLengths1 = (self.nChr[sp1], self.chrLengths[sp1])
            nChr2, chrLengths2 = (self.nChr[sp2], self.chrLengths[sp2])
            distE = self.distEst[combo]
            distS = self.distSim[combo]

            # x = {..., ck: [sbs], ....} with ck chromosomes of sp2 that contain sbs of chromosome co
            # ci is thus a kind of number of homologous chromosomes
            ci = sum(len(c2ToSbsSharedWithc1.keys()) for (c1withSbs, c2ToSbsSharedWithc1) in sbsInPairComp.items())
            # TODO it might be interesting to do it also for M2006PropToLength (called M2006)
            (distE['nTransl'], distE['nInv']) = self.calcTranslInvEstimation(distE['nSB'], nChr1, nChr2, ci,
                                                                            fM2006=self.M2006Unif)
            self.M2006Input[combo]['ci'] = ci
            if self.withBA:
                (distS['nTranslBA'], distS['nInvBA']) = self.calcTranslInvEstimation(distS['nSBBA'], nChr1, nChr2, ciBA,
                                                                                    fM2006=self.M2006Unif)
                self.M2006Input[combo]['ciBA'] = ciBA

            # print progress
            progressBar.printProgressIn(sys.stderr, iCombo, prefixLabel='PAIRWISE COMPS ')

        ###################################################################################
        # 7) Estimate branch length from estimated distance matrix between extant species
        ###################################################################################
        self.branchesEst = {branch: {} for branch in self.branches}
        for statName in ['nInv', 'nTransl', 'nSB']:
            statsPerBranches = self.calcBranchEstimatesFromDist(speciesTree, self.distEst, statName)
            for branch in self.branches:
                self.branchesEst[branch][statName] = statsPerBranches[branch]
        # This second part is for calculating the branch estimates for the theoretical dist values to know
        # the error of LSE
        suffixes = ['', 'BA'] if self.withBA else ['']
        if logErrPath is not None:
            for statBaseName in ('nInv', 'nTransl', 'nSB'):
                for suffix in suffixes:
                    statName = statBaseName + suffix
                    # FIXME solve
                    # no 'nInv' in self.distSim[combo]
                    statsPerBranches = self.calcBranchEstimatesFromDist(speciesTree, self.distSim, statName)
                    for branch in self.branches:
                        self.branchesSim[branch][statName] = statsPerBranches[branch]

        ################################################################################
        # 8) getGeneEventRates per Branch
        ################################################################################
        # TODO understand and check that it works
        self.calcGeneEventRatesOnBranches(speciesTree, loadedGenomes, loadedFamilies)

        ################################################################################
        # 9) Compare real values to estimation
        ################################################################################
        if logErrPath is not None:
            self.compareSimulationAndEstimation()
            # FIXME ???? this load self.error, but I do not understand its meaning... Just give a comprehensible name ?
            # self.calcErrorSplitting()

        print >> sys.stderr, 'myScore object creation finished after {0}s (used by PhylDiag: {1}s)!'.format(
            round(time.clock() - startTime, 1), round(phylDiagTime, 1))
        # DEBUG !!!!!!
        #self.save('/home/jlucas/test.score')
        #print >> sys.stderr, 'DEBUG : myScore save in /home/jlucas/test.score'
        # Switch on stderr stream
        if not verbose:
            sys.stderr = stderrOld
            sys.stdout = stdoutOld

def calcOptimalPhylDiagConfig(modernGenomes, geneFamilies, speciesTreePath,
                              score, combos='all', filterType=myDiags.FilterType[1], gapMaxValues=xrange(0, 12, 1),
                              pValueValues=[0.001, 1]):
    if combos == 'all':
        combos = score.combos
    genomesAreSimulated = len(score.distSim) != 0
    simulatedBP = {}
    if genomesAreSimulated:
        simulatedBP = {key: score.distSim[key]['nSB'] for key in score.distSim}
    speciesTree = myPhylTree.PhylogeneticTree(speciesTreePath)
    comboSel = xrange(len(combos))
    estimatedBP = {}
    detectableBP = {}
    for combo in comboSel:
        spec1 = combos[combo][0]
        # recGenome1 = dirWithSimGenomes + prefixGene + spec1 + '.list.bz2'
        pathGenome1 = modernGenomes.replace('%s', '{0}').format(spec1.replace(' ', '.'))
        genome1 = myLightGenomes.LightGenome(pathGenome1)
        # Delete Y Chromosome
        if 'Y' in genome1.keys():
            genome1.pop('Y')

        spec2 = combos[combo][1]
        pathGenome2 = modernGenomes.replace('%s', '{0}').format(spec2.replace(' ', '.'))
        genome2 = myLightGenomes.LightGenome(pathGenome2)
        # Delete Y Chromosome
        if 'Y' in genome2.keys():
            genome2.pop('Y')

        mrca = speciesTree.lastCommonAncestor([spec1, spec2])
        familyMRCA = myLightGenomes.Families(geneFamilies.replace('%s', '{0}').format(mrca))

        # Detectable breakpoints since common ancestor
        if genomesAreSimulated:
            # In this case one can compare to simulated ancestor,
            # and hence calculate the branch lengthes in breakpoints independently
            pathGenomeAnc = modernGenomes.replace('%s', '{0}').format(mrca)
            genomeAnc = myLightGenomes.LightGenome(pathGenomeAnc)

            sbsInPairComp = myDiags.extractSbsInPairCompGenomes(genome1, genomeAnc, familyMRCA, tandemGapMax=10, gapMax=2,
                                                                distanceMetric='CD', pThreshold=1,
                                                                truncationMax=1000, filterType=filterType,
                                                                minChromLength=2, sameStrand=True,
                                                                nbHpsRecommendedGap=2, targetProbaRecommendedGap=0.01,
                                                                validateImpossToCalc_mThreshold=3,
                                                                distinguishMonoGenicDiags=True,
                                                                identifyMonoGenicInvs=True,
                                                                identifyMicroRearrangements=True,
                                                                optimisation=None,
                                                                verbose=True)
            nBP1 = -1
            for (key1, cont1) in sbsInPairComp.items():
                for (key2, cont2) in cont1.items():
                    nBP1 += len(cont2)

            sbsInPairComp = myDiags.extractSbsInPairCompGenomes(genome2, genomeAnc, familyMRCA, tandemGapMax=10, gapMax=2,
                                                                distanceMetric='CD', pThreshold=1,
                                                                truncationMax=1000, filterType=filterType,
                                                                minChromLength=2, sameStrand=True,
                                                                nbHpsRecommendedGap=2, targetProbaRecommendedGap=0.01,
                                                                validateImpossToCalc_mThreshold=3,
                                                                distinguishMonoGenicDiags=True,
                                                                identifyMonoGenicInvs=True,
                                                                identifyMicroRearrangements=True,
                                                                optimisation=None,
                                                                verbose=True)

            nBP2 = -1
            for (key1, cont1) in sbsInPairComp.items():
                for (key2, cont2) in cont1.items():
                    nBP2 += len(cont2)

            detectableBP[(spec1, spec2)] = (nBP1, nBP2)

        # Observed distance between two species for differing gapMax
        estimatedBP[(spec1, spec2)] = {}
        for gapMax in gapMaxValues:
            estimatedBP[(spec1, spec2)][gapMax] = {}
            for pValue in pValueValues:
                sbsInPairComp = myDiags.extractSbsInPairCompGenomes(genome1, genome2, familyMRCA, tandemGapMax=10,
                                                                    gapMax=gapMax, distanceMetric='CD', pThreshold=pValue,
                                                                    identifyMicroRearrangements=True,
                                                                    truncationMax=1000, filterType=filterType, minChromLength=2,
                                                                    sameStrand=True, nbHpsRecommendedGap=2,
                                                                    targetProbaRecommendedGap=0.01,
                                                                    validateImpossToCalc_mThreshold=3,
                                                                    verbose=True)

                synDis = []
                for (key1, cont1) in sbsInPairComp.items():
                    for (key2, cont2) in cont1.items():
                        for syn in cont2:
                            synDis.append(len(syn.la))
                # plt.hist(sbsDistrib, bins=range(min(sbsDistrib), max(sbsDistrib) + 1, 1)); plt.show()

                nBreakPointsEstimated = len(synDis) - 1
                estimatedBP[(spec1, spec2)][gapMax][pValue] = nBreakPointsEstimated

    # Calculate optimal gapmax (in the sense that the observed number of breakpoints is close to the simulated)
    optimalGap = {}
    if genomesAreSimulated:
        for (key, x) in estimatedBP.items():
            optimalGap[key] = {}
            for pValue in pValueValues:
                min = 100000
                minX = 0
                for gapMax in gapMaxValues:
                    a = abs(x[gapMax][pValue] - sum(detectableBP[key]))
                    if a < min:
                        min = a
                        minX = gapMax
                optimalGap[key][pValue] = minX

    out = {'estimatedBP': estimatedBP, 'detectableBP': detectableBP,
           'simulatedBP': simulatedBP, 'optimalGap': optimalGap,
           'gapMaxValues': gapMaxValues, 'genomesAreSimulated': genomesAreSimulated,
           'FilterType': str(filterType), 'gapMaxValues': [i for i in gapMaxValues]}

    print('PhylDiag - Optimal Parameter Finder finished !')
    return (out)

def calcOptimalPhylDiagConfig_LinePlot(optConfig, pValue, saveImageToPath='', ymax=1400):
    from itertools import cycle
    lines = ["-", "--", "-.", ":"]
    colors = ["r", "g", "b", 'black', 'lightblue']
    linecycler = cycle(lines)
    colorcycler = cycle(colors)
    labels = []
    legendEntries = []
    fig = plt.figure(figsize=(16.5, 10.5))
    ax = fig.add_subplot(111)
    for (key, x) in optConfig['estimatedBP'].items():
        color = next(colorcycler)
        if optConfig['genomesAreSimulated']:
            plt.plot(optConfig['optimalGap'][key][pValue], sum(optConfig['detectableBP'][key]), color=color, marker='x')
            plt.plot(optConfig['optimalGap'][key][pValue], optConfig['simulatedBP'][key], color=color, marker='o')
        thisPlot, = plt.plot(optConfig['gapMaxValues'], [vDic[pValue] for (gapMax, vDic) in x.items()], next(linecycler), color=color)
        legendEntries.append(thisPlot)
        label = key[0].split(' ')[0][:1] + key[0].split(' ')[1][:1] + ' - ' + key[1].split(' ')[0][:1] + key[1].split(' ')[1][:1]
        labels.append(label)
    plt.legend(legendEntries, labels, loc=1, ncol=2, borderaxespad=0.)
    myPlotTitle = 'PhylDiag optimal gapMax with pValue = ' + str(pValue) + ', FilterType = ' + optConfig['FilterType']
    if optConfig['genomesAreSimulated']:
        myPlotTitle += '\nsimulated breakpoints (o), detectable breakpoints (x)'
    plt.title(myPlotTitle, fontsize=20)
    plt.xticks(range(min(optConfig['gapMaxValues']), max(optConfig['gapMaxValues']) + 1, 1))
    plt.xlim(min(optConfig['gapMaxValues']), max(optConfig['gapMaxValues']))
    plt.ylim(0, ymax)
    plt.xlabel('Allowed gapmax', fontsize=16)
    plt.ylabel('Number of breakpoints', fontsize=16)
    fig = plt.gcf()
    fig.set_size_inches(16.5, 10.5)
    if saveImageToPath != '':
        fig.savefig(saveImageToPath, dpi=150)
    plt.show()

def calcOptimalPhylDiagConfig_ScatterPlot(optConfig, pValue, saveImageToPath='', ymax=1400):
    labels = []
    for key in optConfig['estimatedBP']:
        label = key[0].split(' ')[0][:1] + key[0].split(' ')[1][:1] + ' - ' + key[1].split(' ')[0][:1] + key[1].split(' ')[1][:1]
        labels.append(label)
    plt.figure(figsize=(8, 7))
    plt.scatter(range(len(labels)), [sum(optConfig['detectableBP'][key]) for key in optConfig['detectableBP'].keys()],
                marker='x', color='r')
    plt.scatter(range(len(labels)),
                [optConfig['simulatedBP'][key] for key in optConfig['detectableBP'].keys()], marker='o',
                color='b')
    plt.scatter(range(len(labels)),
                [optConfig['estimatedBP'][key][1][pValue] for key in optConfig['detectableBP'].keys()], marker='^',
                color='g')
    plt.xticks(range(len(labels)), [x.replace(' - ', '\n') for x in labels], size='small')
    # plt.legend(loc=1, ncol=2, borderaxespad=0.)
    plt.title('Simulated breakpoints (o), detectable BP (x), observed BP (^)\nwith pValue = ' + str(pValue) + ', gapMax = ' +str(1)+ ', Filter = ' +optConfig['FilterType'],
              fontsize=16)
    plt.xlabel('Species comparison', fontsize=16)
    plt.xlim(-1, len(labels))
    plt.ylim(0, ymax)
    plt.ylabel('Number of breakpoints', fontsize=16)
    fig = plt.gcf()
    fig.set_size_inches(8, 7)
    if saveImageToPath != '':
        fig.savefig(saveImageToPath, dpi=150)
    plt.show()

def calcOptimalPhylDiagConfig_SimRealCompPlot(optConfigSim, optConfigReal, pValue, saveImageToPath='', ymax=2):
    from itertools import cycle
    lines = ["-", "--", "-.", ":"]
    colors = ["r", "g", "b", 'black', 'lightblue']
    linecycler = cycle(lines)
    colorcycler = cycle(colors)
    labels = []
    legendEntries = []
    fig = plt.figure(figsize=(16.5, 10.5))
    ax = fig.add_subplot(111)
    for (key, x) in optConfigSim['estimatedBP'].items():
        color = next(colorcycler)
        optimalGap = optConfigSim['optimalGap'][key][pValue]
        # plt.plot(optimalGap, float(sum(optConfigSim['detectableBP'][key]))/optConfigReal['estimatedBP'][key][optimalGap][pValue], color=recColor, marker='x')
        plt.plot(optimalGap, float(optConfigSim['simulatedBP'][key])/optConfigReal['estimatedBP'][key][optimalGap][pValue], color=color, marker='o')
        simRealRatio = [float(vDic[pValue])/optConfigReal['estimatedBP'][key][gapMax][pValue] for (gapMax, vDic) in x.items()]
        thisPlot, = plt.plot(optConfigSim['gapMaxValues'], simRealRatio, next(linecycler), color=color)
        legendEntries.append(thisPlot)
        label = key[0].split(' ')[0][:1] + key[0].split(' ')[1][:1] + ' - ' + key[1].split(' ')[0][:1] + key[1].split(' ')[1][:1]
        labels.append(label)
    plt.legend(legendEntries, labels, loc=1, ncol=2, borderaxespad=0.)
    myPlotTitle = 'BP_Simulation / BP_Real (measured with PhylDiag, pValue = ' + str(pValue) + ', FilterType = ' + optConfigSim['FilterType']
    # myPlotTitle += '\nsimulated breakpoints (o), detectable breakpoints (x)'
    myPlotTitle += ',\nsimulated breakpoints / BP_real (o)'
    plt.title(myPlotTitle, fontsize=20)
    plt.xticks(range(min(optConfigSim['gapMaxValues']), max(optConfigSim['gapMaxValues']) + 1, 1))
    plt.xlim(min(optConfigSim['gapMaxValues']), max(optConfigSim['gapMaxValues']))
    plt.ylim(0, ymax)
    plt.xlabel('Allowed gapmax', fontsize=16)
    plt.ylabel('Number of breakpoints', fontsize=16)
    fig = plt.gcf()
    fig.set_size_inches(16.5, 10.5)
    if saveImageToPath != '':
        fig.savefig(saveImageToPath, dpi=150)
    plt.show()

def calcBreakpoints(spec1, spec2, family, pThreshold = 1, tandemGapMax = 10, gapMax = 0, identifyMicroRearrangements = True, truncationMax = 1000, filterType = myDiags.FilterType[1]):
    stderrOld = sys.stderr
    sys.stderr = open(os.devnull, 'wb')

    spec1 = myLightGenomes.LightGenome(spec1)
    spec2 = myLightGenomes.LightGenome(spec2)
    if 'Y' in spec1.keys():
        spec1.pop('Y')
    if 'Y' in spec2.keys():
        spec2.pop('Y')
    families = myLightGenomes.Families(family)

    sbsInPairComp = myDiags.extractSbsInPairCompGenomes(spec1, spec2, families,
                                                            distanceMetric='CD', pThreshold=pThreshold, tandemGapMax=tandemGapMax, gapMax=gapMax,
                                                            identifyMicroRearrangements=identifyMicroRearrangements,
                                                            truncationMax=truncationMax, filterType=filterType,
                                                            minChromLength=2, sameStrand=True, nbHpsRecommendedGap=2,
                                                            targetProbaRecommendedGap=0.01, validateImpossToCalc_mThreshold=3,
                                                            verbose=False)

    synDis = []
    for (key1, cont1) in sbsInPairComp.items():
        for (key2, cont2) in cont1.items():
            for syn in cont2:
                synDis.append(len(syn.la))
    nBreakPointsEstimated = len(synDis) - 1
    sys.stderr = stderrOld
    return(nBreakPointsEstimated)

def calcTandemDupsFromGenomes(modernGenomes, geneFamilies, species, tandemGapMax=70):
    family = myLightGenomes.Families(geneFamilies)
    genomeList = []
    for spec in species:
        genomeList.append(myLightGenomes.LightGenome(modernGenomes.replace('%s', '{0}').format(spec.replace(' ', '.'))))

    t = np.arange(0, tandemGapMax, 1)

    tandemDupList = {}
    for (speci, genome) in enumerate(genomeList):
        print('Finished ' + str(round(float(speci) / len(genomeList) * 100, 0)) + '%')
        (distDict, genomeFilt) = myMapping.calcDupDist(genome, family)
        tandemDs = [myMapping.calcTandemDupFromDist(distDict, gapMax=i) for i in t]
        totalD = myMapping.calcTotalDupFromDist(distDict)

        # tandemDups = [myMapping.calcTandemDup(genome, family, i) for i in t]
        #
        tandemDupRates = [float(tD) / totalD for tD in tandemDs]
        tandemDupList[species[speci]] = tandemDupRates

    out = {'tandemDupList': tandemDupList, 'geneFamilies': geneFamilies.split('/')[-1].split('.')[1]}
    print('Tandem duplication calculation finished !')
    return (out)

def calcTandemDup_Plot(tandemDup, saveImageToPath=''):
    from itertools import cycle

    lines = ["-", "--", "-.", ":"]
    colors = ["r", "g", "b", 'black', 'lightblue']
    linecycler = cycle(lines)
    colorcycler = cycle(colors)
    plt.figure(figsize=(16.5, 10.5))
    species = []
    for (key, tandemDupRates) in tandemDup['tandemDupList'].iteritems():
        species.append(key)
        plt.plot([tandemDupRate for tandemDupRate in tandemDupRates], next(linecycler), color=next(colorcycler))

    plt.legend(species, loc=4, ncol=2, borderaxespad=0.)
    plt.ylim(0, 1)
    plt.title('Tandem-Duplication rates for gene-families: ' + tandemDup['geneFamilies'], fontsize=20)
    plt.xlabel('Allowed gap size', fontsize=16)
    fig = plt.gcf()
    fig.set_size_inches(16.5, 10.5)
    if saveImageToPath != '':
        fig.savefig(saveImageToPath, dpi=150)
    plt.show()

if __name__ == '__main__':

    pathMagSimus = 'MagSimus/'
    # pathMagSimus = '/home/jlucas/Libs/MagSimus'

    pathRealGenomes = pathMagSimus + 'data/genesST.%s.list.bz2'
    # pathRealGenomes = 'data/genesST.%s.list.bz2'
    pathRealAncGenes = pathMagSimus + 'data/ancGenes.%s.list.bz2'
    # pathRealAncGenes = 'data/ancGenes.%s.list.bz2'
    pathSpeciesTree = pathMagSimus + 'data/speciesTree.phylTree'
    # pathSpeciesTree = 'data/speciesTree.phylTree'

    # real = myScore(pathRealGenomes,
    #                pathRealAncGenes,
    #                pathSpeciesTree,
    #                phylDiagGapMax=5,
    #                verbose=True, ID='real')
    # # record score on hard-drive
    # real.save('./realScore')
    # # load a precomputed score
    real = myScore.load('./scoreReal')
    print >> sys.stderr, 'realScore Loaded'

    pathSimExtantGenomes = pathMagSimus + 'res/simu1/genes.%s.list.bz2'
    pathSimAncGenes = pathMagSimus + 'res/simu1/ancGenes.%s.list.bz2'
    logErrPath = pathMagSimus + 'res/simu1/logErr'
    print >> sys.stderr, 'simScore Loaded'

    # pathSimExtantGenomes = 'Opt9/0/1/genes.%s.list.bz2'
    # pathSimAncGenes = 'Opt9/0/1/ancGenes.%s.list.bz2'
    # logErrPath = 'Opt9/0/1/0-1.stderr'

    # sim = myScore(pathSimExtantGenomes,
    #               pathSimAncGenes,
    #               pathSpeciesTree,
    #               phylDiagGapMax=5,
    #               logErrPath=logErrPath,
    #               verbose=True, ID='1')
    # # record score on hard drive
    # sim.save('./scoreSim')
    # # load a precomputed score
    sim = myScore.load('./scoreSim')

    myComp = myScoreComparison(sim, real)
    print >> sys.stderr, '1st simScoreComp Loaded'
    myComp2 = copy.deepcopy(myComp)
    for (spec, subScores) in myComp2.dist.iteritems():
        for (subScoreName, subSubScores) in subScores.iteritems():
            for (subSubSubScore, val) in subSubScores.iteritems():
                myComp2.dist[spec][subScoreName][subSubSubScore] += np.random.normal(scale=40)
    myComp3 = copy.deepcopy(myComp)
    myComp3.ID = '2-real'

    mca = myScoreComparisonArchive()
    mca.addComp(myComp)
    mca.addComp(myComp2)
    mca.addComp(myComp3)
    mca.calcStatistics()
    #mca.plotScoreDistribution(mca.scores, 'dist-nInv-DiffRatio', markComps=[myComp.ID], saveImageToPath='./test1.svg') #, species = [('Homo sapiens', 'Mus musculus')]
    #mca.plotScoreDistribution(mca.scores, 'dist-nInv-DiffRatio', markComps=[myComp.ID], species=[('Homo sapiens', 'Mus musculus')], saveImageToPath='./test2.svg')
    mca.plotScoreDistribution(mca.scores, 'dist-nSB-Diff', markComps=[], saveImageToPath='./test3.png')
    print >> sys.stderr, 'dist-nSB-Diff distrib wrote'

    if False:
        # Calculate number of breakpoints dependent on gapmax and pValue
        # This takes a while (30min per launch)
        optConfig_sim = calcOptimalPhylDiagConfig(pathSimExtantGenomes,
                                                  pathSimAncGenes,
                                                  pathSpeciesTree,
                                                  sim, pValueValues=[0.001, 1], gapMaxValues=xrange(0, 12, 1))
        optConfig_real = calcOptimalPhylDiagConfig(pathRealGenomes,
                                                   pathRealAncGenes,
                                                   pathSpeciesTree,
                                                   real, pValueValues=[0.001, 1], gapMaxValues=xrange(0, 12, 1))
        exit()

        # Plot the result
        imagePath = pathMagSimus + "PhylDiagOptimalGapMax-%s_pValue%s_InBoth.png"

        calcOptimalPhylDiagConfig_LinePlot(optConfig_sim, 1, ymax=1600, saveImageToPath=(imagePath % ('Simulated', 1)))
        calcOptimalPhylDiagConfig_LinePlot(optConfig_sim, 0.001, ymax=1600, saveImageToPath=(imagePath % ('Simulated', 0.001)))
        calcOptimalPhylDiagConfig_LinePlot(optConfig_real, 1, ymax=1600, saveImageToPath=(imagePath % ('Real', 1)))
        calcOptimalPhylDiagConfig_LinePlot(optConfig_real, 0.001, ymax=1600, saveImageToPath=(imagePath % ('Real', 1)))

        calcOptimalPhylDiagConfig_SimRealCompPlot(optConfig_sim, optConfig_real, 1, ymax=1.3, saveImageToPath=(imagePath % ('SimRealComp', 1)))
        calcOptimalPhylDiagConfig_SimRealCompPlot(optConfig_sim, optConfig_real, 0.001, ymax=1.3, saveImageToPath=(imagePath % ('SimRealComp', 0.001)))

        # Plot to compare simulated, observed and detectable breakpoints
        calcOptimalPhylDiagConfig_ScatterPlot(optConfig_sim, 1)


    # Test effect of identifyMicroRearrangements in HS - MM
    for gapMax in xrange(0, 7):
        iBWG_T = calcBreakpoints(pathRealGenomes % 'Homo.sapiens',
                                 pathRealGenomes % 'Mus.musculus',
                                 pathRealAncGenes % 'Euarchontoglires',
                                 identifyMicroRearrangements=True, gapMax=gapMax)
        iBWG_F = calcBreakpoints(pathRealGenomes % 'Homo.sapiens',
                                 pathRealGenomes % 'Mus.musculus',
                                 pathRealAncGenes % 'Euarchontoglires',
                                 identifyMicroRearrangements=False, gapMax=gapMax)
        print('gapMax: ' + str(gapMax) + ', BP with identifyMicroRearrangements=T: ' + str(iBWG_T) + ' and without: ' + str(iBWG_F))

    ################################################################################
    # Tandem duplication plot
    ################################################################################
    species = ['Homo.sapiens', 'Mus.musculus', 'Canis.lupus.familiaris', 'Monodelphis.domestica', 'Gallus.gallus',
               'Rattus.norvegicus', 'Callithrix.jacchus', 'Oryctolagus.cuniculus', 'Pan.troglodytes', 'Pongo.abelii',
               'Gorilla.gorilla.gorilla', 'Sus.scrofa', 'Macaca.mulatta', 'Papio.anubis', 'Felis.catus',
               'Equus.caballus']
    tandemDup_all = calcTandemDupsFromGenomes(pathMagSimus + 'Resv78/GenesST/genesST.%s.list.bz2',
                                              pathMagSimus + 'Resv78/AncGenes/ancGenes.Amniota.list.bz2',
                          species=species)

    tandemDup_HS_MM = calcTandemDupsFromGenomes(pathMagSimus + 'Resv78/GenesST/genesST.%s.list.bz2',
                                  pathMagSimus + 'Resv78/AncGenes/ancGenes.Euarchontoglires.list.bz2',
                                  species=['Homo.sapiens', 'Mus.musculus'])

    calcTandemDup_Plot(tandemDup_all, saveImageToPath=pathMagSimus + 'TandemDupMultSpec.png')
    calcTandemDup_Plot(tandemDup_HS_MM, saveImageToPath=pathMagSimus + 'TandemDup_HS_MM.png')

    # The tandemDupRatio is higher in modern genes, yet at gap == 0 the effect is marginal. Yet, it increases a lot more in close range
    # FIXME
    [tandemDup_HS_MM['Homo.sapiens'][i] - tandemDup_all['Homo.sapiens'][i] for i in xrange(len(tandemDup_HS_MM['Homo.sapiens']))]
    [tandemDup_HS_MM['Mus.musculus'][i] - tandemDup_all['Mus.musculus'][i] for i in xrange(len(tandemDup_HS_MM['Homo.sapiens']))]

    # Extreme case
    myMapping.calcTandemDup(myLightGenomes.LightGenome(pathMagSimus + 'Resv76/GenesST/genesST.Callithrix.jacchus.list.bz2'), myLightGenomes.Families(pathMagSimus + 'Resv76/AncGenes/ancGenes.Amniota.list.bz2'), 10000)
    myMapping.calcTandemDup(myLightGenomes.LightGenome(pathMagSimus + 'Resv76/GenesST/genesST.Sus.scrofa.list.bz2'), myLightGenomes.Families(pathMagSimus + 'Resv76/AncGenes/ancGenes.Amniota.list.bz2'), 10000)

    # Comp KS and AbsoluteFitting
    # Usually gives the same order
    # FIXME
    [myComp1.chrLengths[spec]['AbsoluteMeanFitting'] for spec in myComp1.chrLengths]
    [myComp1.chrLengths_KS[spec]['Statistic'] for spec in myComp1.chrLengths_KS]

    # python MagSimus/src/analysis/myOptimizer.py Opt MagSimus/data/genesST.%s.list.bz2 MagSimus/data/ancGenes.%s.list.bz2 MagSimus/data/speciesTree.phylTree MagSimus/src/magSimus1.py -startSpecRatesPath="Opt4/4/specRates_MS1.v80-4" -startParameterFilePath="Opt4/4/parameters.v80-4" -scoreCalculationVerbose=False -realScorePath=MagSimus/res/real.score -executionMode=2 -nRepPerInput=10 -specRatesChangeMode="constant"

# from myScore import *
# import utils.myMapping as myMapping
# import utils.myLightGenomes as myLightGenomes
# mSDict = {}
# for i in range(1, 10):
#     mSDict[i] =myScore('Opt6/1/'+str(i)+'/genes.%s.list.bz2', 'Opt6/1/'+str(i)+'/ancGenes.%s.list.bz2',
#                   'Opt6/realData/speciesTree.phylTree', 'Opt6/1/'+str(i)+'/1-'+str(i)+'.stderr', verbose = True)
#
# speciesTree = myPhylTree.PhylogeneticTree('Opt6/realData/speciesTree.phylTree')
# rs = myScore.load('Opt6/realData/real.score')
# for recS in mSDict.itervalues():
#     print(speciesTree.parent['Homo sapiens'].distance * (recS.geneEventRatesOnBranches['Homo sapiens']['geneDup']))
#
# print(speciesTree.parent['Homo sapiens'].distance * rs.geneEventRatesOnBranches['Homo sapiens']['geneDup'])
# eu = myLightGenomes.Families('Opt6/1/4/ancGenes.Euarchontoglires.list.bz2')
# hs = myLightGenomes.LightGenome('Opt6/1/4/genes.Homo.sapiens.list.bz2')
# myScore.calcGeneEvents(hs,eu)
# myMapping.calcNumberOfGeneDuplications(hs, eu)

# mca = myScoreComparisonArchive.createFromDir('Opt6')
# mca.rankList['dist-nBreakpoints-Diff']['Cumulated']
#
# realScore = myScore.load('Sim5204/realData/real.score')
#
# from collections import Counter
# from collections import OrderedDict
# import analysis.getTandemDuplications as getTandemDuplications
# geneDupGapsCount = {}
# for comp in realScore.sbsDistrib.iterkeys():
#     textBest = 'Best (KS: '+str(round(mca5204.scores['sbsDistrib_KS-LogStatistic'][comp]['mean'],2))+', AMF: '+str(round(mca5204.scores['sbsDistrib-AbsoluteMeanFitting'][comp]['mean'],4))+')'
#     textWorst = 'Worst (KS: '+str(round(mca209.scores['sbsDistrib_KS-LogStatistic'][comp]['mean'],2))+', AMF: '+str(round(mca209.scores['sbsDistrib-AbsoluteMeanFitting'][comp]['mean'],4))+')'
#     geneDupGapsCount[comp] = OrderedDict()
#     geneDupGapsCount[comp]['real'] = Counter(realScore.sbsDistrib[comp])
#     geneDupGapsCount[comp][textBest] = Counter(bestScore.sbsDistrib[comp])
#     geneDupGapsCount[comp][textWorst] = Counter(worstScore.sbsDistrib[comp])
#     getTandemDuplications.plotCurves(geneDupGapsCount[comp], xmin=0, xmax = 50, normalized = False, plotCDF=True,
#                                              xTitle='Synteny Block size in genes', yTitle='Number of blocks', plotTitle=comp,
#                                              saveImageToPath='C:/Users/Luc/Desktop/Meeting/synDis_'+'-'.join([compL[0][0] + compL[1][0] for compL in [comp.split(' ') for comp in comp]])+'_changingNInv.png')
#     getTandemDuplications.plotCurves(geneDupGapsCount[comp], xmin = 0, xmax = 50, normalized = True,
#                                              plotCDF = True,
#                                              xTitle = 'Synteny Block size in genes', yTitle = 'Proportion of blocks',
#                                              plotTitle = comp,
#                                              saveImageToPath = 'C:/Users/Luc/Desktop/Meeting/synDisRel_' + '-'.join(
#                                                  [compL[0][0] + compL[1][0] for compL in
#                                                   [comp.split(' ') for comp in comp]]) + '_changingNInv.png')
#
# comp = ('Homo sapiens', 'Mus musculus')
# getTandemDuplications.plotCurves(geneDupGapsCount[comp], xmin = 0, normalized = False,
#                                          plotCDF = True, ymax = 500,
#                                          xTitle = 'Synteny Block size in genes', yTitle = 'Number of blocks',
#                                          plotTitle = comp,
#                                          saveImageToPath = 'C:/Users/Luc/Desktop/Meeting/synDisRel_' + '-'.join(
#                                              [compL[0][0] + compL[1][0] for compL in
#                                               [comp.split(' ') for comp in comp]]) + '_changingNInv_noXMax.png')
