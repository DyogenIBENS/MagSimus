#!/usr/bin/python

# MagSimus/src/analysis
# python 2.7.8
# Author: Lucas TITTMANN
# Copyright (c) 2015 IBENS/Dyogen Lucas TITTMANN, Joseph LUCAS and Hugues ROEST CROLLIUS
# mail : hrc(at)ens.fr
# Licence: GPL 3.0
# This is free software, you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 (GPL v3) or later and the CeCiLL v2 license in France

# Description:
# This script provides functions to analyse and plot simulation results created with myOptimizer
# Among others, there are plots to show the development of specRates during optimization
# Convergence plots for distances and branch lengths
# And plots for the synteny block size distribution

import os
import copy
import numpy as np
# import matplotlib.pyplot as plt
# from pycallgraph import PyCallGraph
# from pycallgraph.output import GraphvizOutput
from collections import defaultdict
from collections import Counter
from collections import OrderedDict
import matplotlib.pyplot as plt
from itertools import cycle
import sys

import analysis.myScore as myScore
import analysis.myOptimizer as myOptimizer
import analysis.getTandemDuplications as getTandemDuplications
import utils.myPhylTree as myPhylTree
import magSimus1


def plotCurves(dictCurvesDictXY, CI=None, saveImageToPath='', normalized=False, plotCumulated=True, xmax=None, ymax=None,
           xmin=0, ymin=0, yLogScale=False, asHists=False,
           plotTitle='Distance of duplication by branch (o: duplications on other chromosome)',
           yTitle='Number of duplications', xTitle='Gap between duplications', showPlot=True):
    """
    :param dictCurvesDictXY: {'name of plot':  {x1: y1, x2: y2, ...}, 'Reality': {x1: y1, ...}}
    :param CI:
    :param saveImageToPath:
    :param normalized:
    :param plotCumulated:
    :param xmax:
    :param ymax:
    :param xmin:
    :param ymin:
    :param yLogScale:
    :param plotTitle:
    :param yTitle:
    :param xTitle:
    :param showPlot:
    :return:
    """
    # Use \ci in KeyName in dictCurvesDictXY to plot CI


    dictCurvesDictXY = copy.deepcopy(dictCurvesDictXY)
    #lines = ["-", "--", "-.", ":"]
    lines = ["-"]
    colors = ["r", "g", "b", 'c', 'm', 'y', 'k']
    linecycler = cycle(lines)
    colorcycler = cycle(colors)
    labels = []
    legendEntries = []

    plt.ioff()
    fig = plt.figure(figsize=(14.5, 8.5))
    ax = fig.add_subplot(111)
    ax.get_yaxis().get_major_formatter().set_useOffset(False)

    if ymax is None:
        ymaxAuto = True
        ymax = 0
    else:
        ymaxAuto = False

    if xmax is None:
        xmaxAuto = True
        xmax = 0
    else:
        xmaxAuto = False

    curveColor = {}
    for (iCurve, curveName) in enumerate(dictCurvesDictXY):
        if normalized:
            sumYs = sum(dictCurvesDictXY[curveName].values())
            dictCurvesDictXY[curveName] = {x: float(y)/sumYs for (x, y) in dictCurvesDictXY[curveName].iteritems()}
        color = next(colorcycler)
        lineType = next(linecycler)
        if 'Gene' in curveName:
            color = 'blue'
            lineType = "--"
        elif 'Stat' in curveName:
            color = 'red'
            lineType = "-."
        elif 'Block' in curveName:
            color = 'green'
            lineType = ".."
        if curveName == 'Reality':
            color = 'black'
            lineType = "-"
        elif '\ci' in curveName:
            if lineType == '-':
                lineType = next(linecycler)
            ax.fill_between(CI[curveName][0], CI[curveName][1], CI[curveName][2], facecolor=color, interpolate=True, linewidth=0.0, alpha=0.1)

        curveColor[curveName] = color
        curveAsListXY = sorted([(x, y) for (x, y) in dictCurvesDictXY[curveName].iteritems() if x != np.finfo(np.float32).max])
        if plotCumulated:
            Xs = [x for (x,y) in curveAsListXY]
            cumCurve = zip(Xs, list(np.cumsum([y for (x, y) in curveAsListXY])))
            curveAsListXY = cumCurve
        if not asHists:
            thisPlot, = ax.plot([x for (x, y) in curveAsListXY], [y for (x, y) in curveAsListXY],
                                lineType, color = color, linewidth=4.0)
        else:
            width = 1.0 / float(len(dictCurvesDictXY) + 1)
            thisPlot = ax.bar([(x + iCurve * width) for (x, y) in curveAsListXY], [y for (x, y) in curveAsListXY],
                                width=width, color=color)

        legendEntries.append(thisPlot)
        labels.append(curveName.replace('\ci',''))
        if ymaxAuto:
            ymax = max(ymax, max([y for (x, y) in curveAsListXY]))
        if xmaxAuto:
            xmax = max(xmax, max([x for (x, y) in curveAsListXY]))

    # FIXME, understand ...
    for (curveName, curveDictXY) in dictCurvesDictXY.iteritems():
        if np.finfo(np.float32).max in curveDictXY.keys():
            ax.plot(xmax, curveDictXY[np.finfo(np.float32).max], color=curveColor[curveName], marker='o', linewidth=2.0)

    ax.set_title(plotTitle, fontsize=20)
    ax.set_xlim((xmin, xmax+1))
    ax.set_ylim((ymin, ymax*1.05))
    ax.set_xlabel(xTitle, fontsize=16)
    ax.set_ylabel(yTitle, fontsize=16)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1,
                     box.width, box.height * 0.9])
    # ind_real = [i for (i, label) in enumerate(labels) if label == 'Reality'][0]
    # Set reality as second entry in legend
    # labels = [labels[1], labels[0], labels[2:]]
    # legendEntries = [legendEntries[1], legendEntries[0], legendEntries[2:]]
    leg = ax.legend(legendEntries, labels, loc = 'upper center', bbox_to_anchor = (0.5, -0.1),
              fancybox = True, shadow = True, ncol = 2)
    for recLeg in leg.legendHandles:
        recLeg.set_linewidth(7.0)
    if yLogScale:
        ax.set_yscale('log')

    fig.set_size_inches(13, 7.5, forward=True)
    if saveImageToPath != '':
        fig.savefig(saveImageToPath, dpi = 150, format=saveImageToPath.split('.')[-1])
    if showPlot:
        plt.show()
    plt.close(fig)

class myResults:

    def __init__(self, optPath):
        pass
        # if isinstance(path, basestring) and path[-1] != '/':
        #     path += '/'
        # self.optPath = optPath

    @staticmethod
    def calcDiffOfSRVals(sr1, sr2):
        # Return the dict of differences between two 2 level dictionaries
        out = {}
        for key1 in sr1:
            out[key1] = {}
            for key2 in sr1[key1]:
                out[key1][key2] = sr1[key1][key2] - sr2[key1][key2]
        return(out)

    # SR means SpecRates
    @staticmethod
    def calcSRDevelopment(path, speciesTree):
        if isinstance(path, basestring) and path[-1] != '/':
            path += '/'

        res = {}
        sim = 0
        while os.path.isdir(path + str(sim)):
            SRName = path + str(sim) + '/' + [file for file in os.listdir(path + str(sim)) if 'specRates' in file][0]
            SR = myOptimizer.myOptimizer.readSpecRates(SRName)
            SRVals = myOptimizer.myOptimizer.changeRatesIntoAbsValues(SR, speciesTree)
            if sim == 0:
                res['Start'] = SRVals
                lastSRVals = SRVals
                res['Diff'] = []
                res['Rates'] = []
            else:
                diff = myResults.calcDiffOfSRVals(SRVals, lastSRVals)
                res['Diff'].append(diff)
                res['Rates'].append(SR)
                lastSRVals = SRVals
                print >> sys.stderr, SR
                # raw_input()
            sim += 1
        return res

    @staticmethod
    def plotSRDevelopment(results, rateName, saveImageToPath='', mainTitle=None, speciesTree=None,
                          plotNbEvents=False):
        if plotNbEvents:
            assert speciesTree is not None
            assert isinstance(speciesTree, myPhylTree.PhylogeneticTree)
        # Criteria may be chrInvert, chrTransloc, geneDup, geneLoss, geneTandemDupProp, ...

        lines = ["-", "--", "-."]
        colors = ["r", "g", "b", 'c', 'm', 'y', 'k']
        linecycler = cycle(lines)
        colorcycler = cycle(colors)
        labels = []
        legendEntries = []
        fig = plt.figure(figsize = (16.5, 10.5))
        ax = fig.add_subplot(111)
        ax.get_yaxis().get_major_formatter().set_useOffset(False)

        for sp in results['Rates'][0].iterkeys():
            ys = [iSimVals[sp][rateName] for iSimVals in results['Rates']]
            if plotNbEvents:
                # multiply by the branch length
                ys = [y * speciesTree.parent[sp].distance for y in ys]

            color = next(colorcycler)
            thisPlot, = ax.plot(range(len(ys)), ys, next(linecycler), color=color, markersize=2)
            legendEntries.append(thisPlot)

            if isinstance(sp, tuple):
                sp = ' - '.join(sp)
            shortName = sp.split('.')
            shortName = ''.join([word[0] for word in shortName]).upper()
            labels.append(shortName)

        if (mainTitle is None):
            mainTitle = 'Rates of ' + rateName
        ax.set_title(mainTitle, fontsize=20)
        ax.set_xlabel('Simulation', fontsize=16)
        ax.set_ylabel('Score', fontsize=16)
        ax.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.1,
                         box.width, box.height * 0.9])

        # Put a legend below current axis
        ax.legend(legendEntries, labels, loc='upper center', bbox_to_anchor=(0.5, -0.05),
                  fancybox=True, shadow=True, ncol = 5)

        # ax.legend(legendEntries, labels, loc = 2, ncol = 2, borderaxespad = 0.)
        if saveImageToPath != '':
            fig.savefig(saveImageToPath, dpi=150)

    @staticmethod
    def extractInvInputCSV(path, outPath='temp.csv'):
        if isinstance(path, basestring) and path[-1] != '/':
            path += '/'

        (args, opts) = magSimus1.getDefaultParameters()
        parasX = ["chr:invDist", "chr:invDistMean", "chr:invDistShape"]
        outList = []
        for file in os.listdir(path):
            if os.path.isdir(path + file):
                max = 0
                for file2 in os.listdir(path + file):
                    if file2.isdigit() and int(file2) > max:
                        max = int(file2)
                f2c = path + file + '/' + str(max) + '/'
                paraFile = myOptimizer.myOptimizer.readParameterFile(f2c + [file3 for file3 in os.listdir(f2c) if 'para' in file3][0])
                outList = []
                for para in parasX:
                    if para in paraFile:
                        outList.append(str(paraFile[para]))
                    else:
                        outList.append(str(opts[para]))
                outLineStart = '\t'.join(outList) + '\t'
                repDict = {}
                reps = 0
                for file3 in os.listdir(f2c):
                    if file3.isdigit():
                        reps += 1
                        f3c = f2c + str(max) + '/'
                        sc = myScore.myScoreComparison.load(f3c + str(max) + '-' + file3 + '.comp')
                        for comp in sc.synDis_KS.iterkeys():
                            compShort = '-'.join([comp.split('.')[0][0]+comp.split('.')[1][0] for comp in comp])
                            if compShort not in repDict:
                                repDict[compShort] = [[],[],[],[]]
                            repDict[compShort][0].append(str(sc.synDis_KS[comp]['Statistic']))
                            repDict[compShort][1].append(str(sc.synDis_KS[comp]['StatisticBlock']))
                            repDict[compShort][2].append(str(sc.synDis[comp]['AbsoluteMeanFitting']))
                            repDict[compShort][3].append(str(sc.synDis[comp]['MeanFitting']))
                for (comp, vals) in repDict.iteritems():
                    outLine = outLineStart + comp + '\t'
                    outLine += '\t'.join(vals[0]) + '\t'
                    outLine += '\t'.join(vals[1]) + '\t'
                    outLine += '\t'.join(vals[2]) + '\t'
                    outLine += '\t'.join(vals[3])
                    outList.append(outLine)

        out = '\n'.join(outList)
        header = 'Dist\tAlpha\tBeta\tCombo\t'
        varNames = ['KS_Genes','KS_Blocks','AbsoluteMeanFitting','MeanFitting']
        # reps here is from the last simulation
        header += '\t'.join([varName + '-' + str(i+1) for varName in varNames for i in xrange(0, reps)]) + '\n'
        out = header + out
                    
        with open(outPath, "wb") as text_file:
            text_file.write(out)

def calcScoreForBestRep(iSim, bestRepCriteria='sbsDistrib_KS-LogStatistic'):
    iRep = str(np.argmin(sca.archive[bestRepCriteria][('Homo sapiens', 'Mus musculus')][iSim]))
    iSim = iSim.split('-')[0]
    print('Repition ' + iRep + ' is the best by chosen criteria.')
    iPath = 'InvTransl_Sim10_Rep1000_SR/'+iSim+'/'+iRep+'/'
    iScore = myScore.myScore(iPath+'genes.%s.list.bz2', iPath+'ancGenes.%s.list.bz2',
                    'MagSimus/data/speciesTree.phylTree',
                    iPath+iSim+'-'+str(iRep)+'.stderr', verbose=True)
    iScore.save(iPath+iSim+'-'+str(iRep)+'.score')
    return(iScore)

def listOfLoadedScores(path, fileEnding):
    """
    :param path to the directory containing scores
    :param fileEnding: extension (XXX.extension) of the score in ['comp',  'score']
    :return:
    """
    sList = []
    checkList = []
    for file in os.listdir(path):
        if os.path.isdir(path+file+'/'):
            for file2 in os.listdir(path+file+'/'):
                if fileEnding in file2:
                    checkList.append(file)
                    sList.append(myScore.myScore.load(path+file+'/'+file2))
    missingList = [i for i in range(1, 101) if str(i) not in checkList]
    if len(missingList) > 0:
        print('Sim where no score was calculated: ' + str(missingList))
    return(sList)

def lengthsToHistoGramm(lengths, normalizedIn, cumulated=True, increasingOrder=True):
    """
    :param lengths: [..., length, ...]
    :param normalizedIn either 'Gene' or None
    :return: histogramm = [..., (value, count), ...]
    """
    counterLengths = Counter(lengths)
    if normalizedIn is None:
        histogramm = sorted([(length, count) for (length, count) in counterLengths.iteritems()], reverse=(not increasingOrder))
        histogramm = myScore.myScoreComparison.fillHistogrammGapsWithNumber(histogramm, number=0)
        if cumulated:
            Y = [count for (length, count) in histogramm]
            X = [length for (length, count) in histogramm]
            cumHistogramm = zip(X, list(np.cumsum(Y)))
            res = cumHistogramm
        else:
            res = histogramm
    elif normalizedIn == 'Genes':
        # each element has a weight prop to its number of genes in the histogram
        histogramm = sorted([(length, count * length) for (length, count) in counterLengths.iteritems()], reverse=(not increasingOrder))
        histogramm = myScore.myScoreComparison.fillHistogrammGapsWithNumber(histogramm, number=0)
        if cumulated:
            Y = [count for (length, count) in histogramm]
            X = [length for (length, count) in histogramm]
            cumHistogramm = zip(X, list(np.cumsum(Y)))
            res = cumHistogramm
        else:
            res = histogramm
    return res

def plotDF(simulations, pathReal, combo, cumulated=True, asHists=True, normalizedIn=None, alpha=0.05, plotOptRep=False,
             savePath='Meeting/', yLogScale=False, xmax=None, showPlot=True, imageFormat='svg'):




    def createPlotInput(pathSim, label, refLengths, combo, normalizedIn, alpha, cumulated):
        """
        :param refLengths: needed to compute the KS, and to select the best repetition
        :param label: label of the ploted distribution
        :param alpha: the alpha for the CI
        :param normalizedIn either 'Gene' or None
        :param cumulated
        :return:
        """
        refHist = lengthsToHistoGramm(refLengths, normalizedIn, cumulated=cumulated, increasingOrder=True)

        if pathSim not in loadedScores:
            loadedScores[pathSim] = (listOfLoadedScores(pathSim, 'score'), listOfLoadedScores(pathSim, 'comp'))

        # scores of all repets
        sList = loadedScores[pathSim][0]
        # score comps
        scList = loadedScores[pathSim][1]
        hists = []
        # for all scores
        # FIXME, for all combos also
        for s in sList:
            hist = lengthsToHistoGramm(s.sbsLengths[combo], normalizedIn, cumulated=cumulated, increasingOrder=True)
            hists.append(hist)
        print >> sys.stderr, 'scores loaded'

        nReps = len(sList)

        # Make all hists the same length (add some values at the end of some hists)
        maxLength = max([max(hist, key=lambda x: x[0]) for hist in hists])[0]
        histsUnilength = copy.deepcopy(hists)
        for hist in histsUnilength:
            cnt = len(hist)
            while cnt <= maxLength:
                if cumulated:
                    hist.append((cnt, hist[-1][1]))
                else:
                    hist.append((cnt, 0))
                cnt += 1
        print >> sys.stderr, 'Make all hists the same length'

        dictHistCumulatedAllRepets = defaultdict(list)
        for hist in histsUnilength:
            for (length, val) in hist:
                dictHistCumulatedAllRepets[length].append(val)

        lowerCurve = []
        upperCurve = []
        # BUG, does not work if there is only one hist
        if len(sList) > 1:
            lowerBound = int(np.ceil(len(sList)*alpha/2.0))
            upperBound = len(sList) - lowerBound - 1
        else:
            lowerBound = 0
            upperBound = 0
        for length, vals in dictHistCumulatedAllRepets.iteritems():
            lowerCurve.append((length, sorted(vals)[lowerBound]))
            upperCurve.append((length, sorted(vals)[upperBound]))
        print >> sys.stderr, 'low and upper curves'

        maxDiffs = []
        kss = []
        for (iHist, hist) in enumerate(hists):
            if normalizedIn is None:
                pass
            elif normalizedIn == 'Genes':
                refHist = [(length, count * length) for (length, count) in refHist]
                hist = [(length, count * length) for (length, count) in hist]
            elif normalizedIn == 'Blocks':
                pass
            maxDiffs.append((max([np.abs(refHist[idx][1]-hist[idx][1]) for idx in xrange(min(len(refHist), len(hist)))]), iHist))

            if not cumulated:
                kss = maxDiffs
            else:
                X = [length for (length, count) in hist]
                cumRefHist = zip(X, np.cumsum([count for (length, count) in refHist]))
                cumHist = zip(X, np.cumsum([count for (length, count) in hist]))
                kss.append((max([np.abs(cumHist[idx][1] - cumRefHist[idx][1]) for idx in xrange(min(len(cumHist), len(cumRefHist)))]), iHist))
                refHist = cumRefHist
                hist[iHist] = cumHist
        #
        # if normalizedIn is None:
        #     absDiffs = []
        #     for (iHist, hist) in enumerate(hists):
        #         absDiffs.append((max([np.abs(refLengths[idx][1]-hist[idx][1]) for idx in xrange(min(len(refLengths), len(hist)))]), iHist))
        # else:
        #     if normalizedIn == 'Genes':
        #         KS_stat = 'Statistic'
        #
        #         # FIXME
        #     else:
        #         # normalizedIn == 'Blocks':
        #         KS_stat = 'StatisticBlock'
        #         # FIXME
        #    absDiffs = [(sc.synDis_KS[combo][KS_stat], iSC) for (iSC, sc) in enumerate(scList)]
        print >> sys.stderr, 'compute max diffs'

        if plotOptRep:
            # sort from lower ks to higher ks
            (smallerKS, maxKSIndex) = sorted(kss)[0]
        else:
            (smallerKS, maxKSIndex) = sorted(kss)[len(kss)/2]

        optHist = {length: count for (length, count) in hists[maxKSIndex]}
        if plotOptRep:
            histText = 'Optimal repetition %s-%s, KS=%s' % (label, str(maxKSIndex), str(np.round(smallerKS, 2)))
            if not cumulated:
                histText += ' maxDiff=%s' % maxDiffs[maxKSIndex][0]
        else:
            # histText = 'Mean Simulation (Sup = '+str(np.round(minKS,2))+', ' + label + ')'
            histText = label + '-' + str(maxKSIndex) + ' KS=' + str(np.round(smallerKS, 2))
        CI = ([vTuple[0] for vTuple in lowerCurve],
              [vTuple[1] for vTuple in lowerCurve],
              [vTuple[1] for vTuple in upperCurve])
        dictOptHist = {'optHist': optHist, 'histText': histText, 'CI': CI}
        print >> sys.stderr, 'return res'
        return (dictOptHist, nReps)

    def imagePathAndPlotProperties(combo, savePath, normalizedIn, yLogScale, plotOptRep, maxRep, imageFormat='svg'):
        """
        :param yLogScale: bool
        :param plotOptRep:  bool
        :param maxRep: max value of repetitions
        :param normalizedIn: either 'Blocks', 'Genes' or None
        :return:
        """

        comboShortened = '-'.join([c.split('.')[0][0]+c.split('.')[1][0] for c in combo])
        imagePath = savePath +str(comboShortened)
        if cumulated:
            imagePath += '_cdf'
        else:
            imagePath += '_pdf'

        ymin = 0
        if normalizedIn is None:
            yTitle = 'Number of synteny blocks'
            imagePath += '_sbs'
            if yLogScale:
                ymin = 1
                imagePath += '_log'
        else:
            yTitle = 'Proportion of ' + normalizedIn
            imagePath += '_propTo' + normalizedIn
            if yLogScale:
                ymin = 1.0/maxRep
                imagePath += '_log'

        if plotOptRep:
            imagePath += '_opt'
        else:
            imagePath += '_mean'
        imagePath += '.' + imageFormat
        out = {'comboShortened': comboShortened,
               'imagePath': imagePath,
               'ymin': ymin,
               'yTitle': yTitle}
        return out

    sr = myScore.myScore.load(pathReal)
    realSBLengths = sr.sbsLengths[combo]
    refHist = lengthsToHistoGramm(realSBLengths, normalizedIn, cumulated=cumulated, increasingOrder=True)
    dictRealHist = {length: count for (length, count) in refHist}
    inputPlot = {'Reality': dictRealHist}

    CI = {}
    maxRep = 0
    for (pathSimsLabel, pathSim) in simulations.iteritems():
        (simDict, nRep) = createPlotInput(pathSim, pathSimsLabel, realSBLengths, combo, normalizedIn, alpha, cumulated)
        inputPlot[simDict['histText']] = simDict['optHist']
        CI[simDict['histText']] = simDict['CI']
        maxRep = max(maxRep, nRep)

    plotPropDict = imagePathAndPlotProperties(combo, savePath, normalizedIn, yLogScale, plotOptRep, maxRep, imageFormat=imageFormat)
    # return(inputCDF)
    plotCurves(inputPlot, CI=CI,
           xmin=1, ymin=plotPropDict['ymin'], normalized=False,
           xmax=xmax,
           plotCumulated=False, yLogScale=yLogScale,
           asHists=asHists,
           xTitle='Synteny block size in genes', yTitle=plotPropDict['yTitle'],
           plotTitle='Synteny block size distribution\nfor '+str(plotPropDict['comboShortened']) + ' (' + str(maxRep) + ' simulations)',
           saveImageToPath=plotPropDict['imagePath'],
           showPlot=showPlot)

if __name__ == '__main__':

    # # Other Meeting
    # outPath = 'C:/Users/Luc/Desktop/Meeting_2105/'
    # sca = myScore.myScoreComparisonArchive.createFromDir('InvTransl_Sim10_Rep1000_SR')
    # sca.plotScoreDistribution(sca.scores, 'dist-nBreakpoints-Diff', saveImageToPath=outPath+'dist-nBreakpoints-Diff.png')
    # sca.plotScoreDistribution(sca.scores, 'dist-nBreakpoints-Diff', helpLines='minmax', saveImageToPath=outPath+'dist-nBreakpoints-Diff_helpLines=minmax.png')
    # sca.plotScoreDistribution(sca.scores, 'dist-nInv-Diff', saveImageToPath=outPath+'dist-nInv-Diff.png')
    # sca.plotScoreDistribution(sca.scores, 'dist-nTransl-Diff', saveImageToPath=outPath+'dist-nTransl-Diff.png')
    #
    # sca.plotScoreDistribution(sca.scores, 'branches-nInv-Diff', saveImageToPath=outPath+'branches-nInv-Diff.png')
    # sca.plotScoreDistribution(sca.scores, 'branches-nInv-Diff', helpLines='minmax', saveImageToPath=outPath+'branches-nInv-Diff_helpLines=minmax.png')
    # sca.plotScoreDistribution(sca.scores, 'branches-nTransl-Diff', saveImageToPath=outPath+'branches-nTransl-Diff.png')
    # sca.plotScoreDistribution(sca.scores, 'branches-nTransl-Diff', helpLines='minmax', saveImageToPath=outPath+'branches-nTransl-Diff_helpLines=minmax.png')
    #
    # sca.rankList['sbsDistrib_KS-LogStatistic']['Cumulated']['Rel'][0]
    # sca.rankList['dist-nBreakpoints-Diff']['Cumulated']['Rel'][0]
    # sca.rankList['branches-nInv-Diff']['Cumulated']['Rel'][0]
    #
    # bestSim = sca.rankList['sbsDistrib_KS-LogStatistic']['Cumulated']['Rel'][0][0]
    #
    # scoreDict = {}
    # scoreDict['real'] = myScore.myScore.load('Sim5204/realData/real.score')
    # scoreDict['0'] = calcScoreForBestRep('0-real')
    # scoreDict['2'] = calcScoreForBestRep('2-real')
    # scoreDict['30'] = calcScoreForBestRep('30-real')
    #
    # geneDupGapsCount = {}
    # for comp in scoreDict['real'].synDis.iterkeys():
    #     geneDupGapsCount[comp] = OrderedDict()
    #     for (scoreName, score) in scoreDict.iteritems():
    #         if scoreName == 'real':
    #             scoreText = 'real'
    #         else:
    #             scoreText = 'Sim ' + scoreName + '(mean KS: '+str(round(sca.scores['sbsDistrib_KS-LogStatistic'][comp][scoreName+'-real']['mean'], 2))+')'
    #         geneDupGapsCount[comp][scoreText] = Counter(score.synDis[comp])
    #     getTandemDuplications.tandemGap_LinePlot(geneDupGapsCount[comp], xmin=0, normalized=False, plotCDF=True,
    #                                              xTitle='Synteny Block size in genes', yTitle='Number of blocks',
    #                                              plotTitle='SB size distribution '+str(comp)+'\n(for simulations, best repetition is plotted with mean KS in legend)',
    #                                              saveImageToPath=outPath+'synDis_'+'-'.join([compL[0][0] + compL[1][0] for compL in [tmpComp.split('.') for tmpComp in comp]])+'_opt.png')
    #     getTandemDuplications.tandemGap_LinePlot(geneDupGapsCount[comp], xmin=0, normalized=True,
    #                                              plotCDF=True,
    #                                              xTitle='Synteny Block size in genes', yTitle='Proportion of blocks',
    #                                              plotTitle='SB size distribution '+str(comp)+'\n(for simulations, best repetition is plotted with mean KS in legend)',
    #                                              saveImageToPath=outPath + 'synDisRel_' + '-'.join(
    #                                                  [compL[0][0] + compL[1][0] for compL in
    #                                                   [tmpComp.split('.') for tmpComp in comp]]) + '_opt.png')

    # Next
    # sList = listOfLoadedScores('Opt_S15R1000_SR/0/', 'score')
    # scList = listOfLoadedScores('Opt_S15R1000_SR/0/', 'comp')
    # sr = myScore.myScore.load('Opt_S15R1000_SR/realData/real.score')

    # errorMean = {}
    # for s in sList:
    #     for (key, specDict) in s.error.iteritems():
    #         if key not in errorMean:
    #             errorMean[key] = {}
    #         for (spec, statDict) in specDict.iteritems():
    #             if spec not in errorMean[key]:
    #                 errorMean[key][spec] = OrderedDict()
    #             for statName, stat in statDict.iteritems():
    #                 if statName not in errorMean[key][spec]:
    #                     errorMean[key][spec][statName] = 0
    #                 errorMean[key][spec][statName] += float(stat)/len(sList)
    #
    # out = ''
    # for (key, specDict) in errorMean.iteritems():
    #     out += key+'\n'
    #     # Head of table:
    #     out += '\t' + '\t'.join(specDict.itervalues().next().keys()) + '\n'
    #     for (spec, statDict) in specDict.iteritems():
    #         out += str(spec) + '\t' + '\t'.join([str(x) for x in statDict.itervalues()]) + '\n'
    #     out += '\n'
    # with open('v81_SR/error.csv', "w") as text_file:
    #     text_file.write(out)

    loadedScores = {}
    source = 'Grid_2/GridOpt/'
    sims = {'\cistat. estimate': source+'0/',
            '\cisim. estimate': source+'1/'}

    propToXs = [None, 'Genes', 'Blocks']
    plotOptReps = [False]
    plotCDFs = [False, True]
    for combo in myScore.myScore.load(source+'realData/real.score').combos:
        for normalized in propToXs:
            for plotOptRep in plotOptReps:
                for plotCDF in plotCDFs:
                    if not plotCDF:
                        plotDF(sims, source+'realData/real.score', combo, savePath='Simulations/Gamma_v83/', cumulated=plotCDF, yLogScale=False, normalizedIn=normalized, plotOptRep=plotOptRep, xmax=100, showPlot=False)
                    else:
                        plotDF(sims, source+'realData/real.score', combo, savePath='Simulations/Gamma_v83/', cumulated=plotCDF, yLogScale=False, normalizedIn=normalized, plotOptRep=plotOptRep, showPlot=False)
    print('Finished DF plots!')

    # combine images
    from PIL import Image
    outPath = 'Simulations/Gamma_v83/'
    oldCombo = ''
    size = (1025, 600)
    for file in os.listdir(outPath):
        pictureNameParts = file.split('_')
        pathPicture = outPath+file
        yPos = 0
        xPos = 0
        if pictureNameParts[0] != oldCombo:
            if oldCombo != '':
                newPicture.save(outName)
            oldCombo = pictureNameParts[0]
            outName = outPath + 'overview_' + pictureNameParts[0] + ".png"
            newPicture = Image.new('RGB', (size[0]*1, size[1]*3))
        elif 'pdf' in file:
            continue
            yPos += 1
            if 'log' in file:
                xPos += 1
        elif 'opt' in file:
            yPos += 3
        elif 'Blocks' in file:
            yPos += 1
        elif 'Genes' in file:
            yPos += 2

        picture = Image.open(pathPicture)
        picture = picture.resize(size, Image.ANTIALIAS)
        newPicture.paste(picture, (xPos*size[0], yPos*size[1]))
        newPicture.save(outName)

    # Comparison of block optimal and gene optimal simulation
    combo = ('Homo sapiens', 'Mus musculus')
    sims = {'\ciStat. Est.': 'Gamma83_Opt/0/',
            '\ciBlock Opt.': 'Gamma83_Opt/7/',
            '\ciGene Opt.': 'Gamma83_ABCOpt/7/'}
    plotDF(sims, 'MagSimus/res/real.score', combo, cumulated=True, normalizedIn='Blocks', plotOptRep=False, yLogScale=False, showPlot=True, xmax=350)
    plotDF(sims, 'MagSimus/res/real.score', combo, cumulated=True, normalizedIn='Genes', plotOptRep=False, yLogScale=False, showPlot=True, xmax=350)

    # Gamma v83 graphics
    speciesTree = myPhylTree.PhylogeneticTree('MagSimus/data/speciesTree.phylTree')
    srdev=myResults.calcSRDevelopment('Gamma83_Opt', speciesTree)
    myResults.plotSRDevelopment(srdev, 'chrInvert', mainTitle='MagSimus input development during optimization:\nInversion rates per mya')
    myResults.plotSRDevelopment(srdev, 'chrTransloc', mainTitle='MagSimus input development during optimization:\nTranslocation rates per mya')

    sca = myScore.myScoreComparisonArchive.createFromDir('Gamma83_Opt/')
    sca.plotScoreDistribution(sca.scores, 'dist-nInv-Diff', mainTitle='Distance development during optimization:\nInversions').show()

    sca.plotScoreDistribution(sca.scores, 'dist-nInv-Diff', helpLines='minmax', mainTitle='Distance development during optimization:\nInversions (dotted lines indicate min and max)').show()
    sca.plotScoreDistribution(sca.scores, 'dist-nTransl-Diff', mainTitle='Distance development during optimization:\nTranslocations').show()
    sca.plotScoreDistribution(sca.scores, 'dist-nTransl-Diff', helpLines='minmax', mainTitle='Distance development during optimization:\nInversions (dotted lines indicate min and max)').show()

    sca.plotScoreDistribution(sca.scores, 'gMeasure-LogRatio').show()


    # mOG.plotDF(sims, 'Grid_42/GridOpt8/realData/real.score', combo, plotCDF = True, normalizedIn='Blocks', plotOptRep=False, yLogScale = False, showPlot = True, xmax = 350, savePath='./')
    # mOG.plotDF(sims, 'Grid_42/GridOpt5/realData/real.score', combo, plotCDF = True, normalizedIn='Genes', plotOptRep=False, yLogScale = False, showPlot = True, xmax = 350, savePath='./')