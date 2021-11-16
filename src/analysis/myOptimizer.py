#!/usr/bin/python

# MagSimus/src/analysis
# python 2.7.8
# Author: Lucas TITTMANN, except the functions executeCondorMode and executeCondorModeOld, which were written by J. Lucas
# Copyright (c) 2015 IBENS/Dyogen Lucas TITTMANN, Joseph LUCAS and Hugues ROEST CROLLIUS
# mail : hrc(at)ens.fr
# Licence: GPL 3.0
# This is free software, you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 (GPL v3) or later and the CeCiLL v2 license in France

# This script provides an optimization framework for MagSimus. It uses the statistics calculated with
# myScore.py as benchmarks for the quality of the simulated genomes.
# It can either do a (sampled) grid search, optimizing a parameter set, or
# sampling parameters from a distribution and optimizing on them in order to do Approximate Bayesian Computation
# MagSimus simulations can be either calculated locally, also multi-core, or be sent to Condor in order to
# use the computer cluster for faster calculation, using the two functions and the API provided by J. Lucas.

import os
import time
import sys
import shutil
import copy
import itertools
import multiprocessing
import numpy as np
import random
import cPickle as pickle
import traceback
from ast import literal_eval

import analysis.myScore as myScore
import utils.myDiags as myDiags
import utils.myPhylTree as myPhylTree
import utils.myTools as myTools
import utils.myFile as myFile
import utils.myCondor as myCondor
import magSimus1

chrNumbers = {
    'Boreoeutheria': 24,
    'Canis lupus familiaris': 39,
    'Euarchontoglires': 25,
    'Gallus gallus': 30,
    'Homo sapiens': 23,
    'Monodelphis domestica': 9,
    'Mus musculus': 20,
    'Theria': 18,
    'Amniota': 21
}

# Not used
# tandemDupProps = {
#     'Gallus gallus': 0.6508,
#     'Theria': 0.6517,
#     'Monodelphis domestica': 0.4426,
#     'Boreoeutheria': 0.6136,
#     'Canis lupus familiaris': 0.3475,
#     'Euarchontoglires': 0.5410,
#     'Homo sapiens': 0.7984,
#     'Mus musculus': 0.8220
# }

# use estimations of translocations and inversions computed from Mazowita instead of the real number of performed events
doNotUseDiffRatio = True

# FIXME I would rather deactivate the breakpoint analyzer to save time. But this makes the optimizer bug...
# With lazyBPAnalyzer = 'True', the error is
# Traceback (most recent call last):
#   File "Grid_42/GridOpt/execMSandScoreOnCondor.py", line 16, in <module>
#     recScore = optimizer.calculateScore(arguments['iSim'], arguments['iRepet'], saveScore = False)
#   File "/kingdoms/dyogen/workspace4/workspace4/jlucas/Libs/MagSimus/src/analysis/myOptimizer.py", line 115, in calculateScore
#     verbose = self.scoreCalculationVerbose, ID = str(iSim))
#   File "/kingdoms/dyogen/workspace4/workspace4/jlucas/Libs/MagSimus/src/analysis/myScore.py", line 1566, in __init__
#     self.calcErrorSplitting()
#   File "/kingdoms/dyogen/workspace4/workspace4/jlucas/Libs/MagSimus/src/analysis/myScore.py", line 1224, in calcErrorSplitting
#     self.error[branch][key]['LSE'] = self.branchesEst[branch]['nLSE'+toI]
# KeyError: 'nLSETransl'
withBA = True  # Needed to compute scores... :(
saveScore = True


def doRepetitionParallel(args):
    # This function is needed to use the multiprocessing capability of myOptimizer
    # The fundamental problem is that native multiprocessing is using pickle, and
    # methods are not directly pickable, whereas functions are.
    # Hence, this function was created to call the real method myOptimizer.doRepetition()
    # The args-input is a list with the first element being an instance of myOptimizer (self)
    # and the other arguments are the arguments for myOptimizer.doRepetition()
    return args[0].doRepetition(*args[1:])


class myOptimizer:
    def __getstate__(self):
        state = self.__dict__.copy()
        state['speciesTree'] = self.speciesTreePath
        state['scoreReal'] = self.realDataPath + 'real.score'
        return state

    def __setstate__(self, newstate):
        newstate['speciesTree'] = myPhylTree.PhylogeneticTree(newstate['speciesTree'])
        newstate['scoreReal'] = myScore.myScore.load(newstate['scoreReal'])
        self.__dict__.update(newstate)

    # Not used
    # def calcSpecRatesFromScore(self, score, speciesTree):
    #     specRates = copy.deepcopy(score.geneEventRatesOnBranches)
    #     for (branch, chrEvents) in score.branchesEst.iteritems():
    #
    #         # Calculate number of Fusions/Fissions per Branch
    #         chrDiff = chrNumbers[branch] - chrNumbers[speciesTree.parent[branch].name]
    #         fissionRate = fusionRate = 0
    #         if chrDiff > 0:
    #             fissionRate = float(chrDiff) / speciesTree.parent[branch].distance
    #         else:
    #             fusionRate = float(abs(chrDiff)) / speciesTree.parent[branch].distance
    #
    #         specRates[branch]['chrInvert'] = score.branchesEst[branch]['nInv']
    #         specRates[branch]['chrTransloc'] = score.branchesEst[branch]['nTransl']
    #         specRates[branch]['chrFusion'] = fusionRate
    #         specRates[branch]['chrFission'] = fissionRate
    #         specRates[branch]['tandemDupRate'] = tandemDupProps[branch]
    #
    #     return (specRates)

    def calculateScore(self, iSim, iRepet, saveScore=False, usePhylDiagDfltParams=False):
        repetPath = self.optimizationPath + str(iSim) + '/' + str(iRepet) + '/'
        ID = str(iSim) + '-' + str(iRepet)
        # use default parameters of phylDiag
        if usePhylDiagDfltParams:
            tandemDupGap = 5
            phylDiagGapMax = 5
            phylDiagPValue = None
            phylDiagTruncationMax = 5
            phylDiagFilter = myDiags.FilterType.InBothGenomes
        else:
            # Same parameters than for the calculation of realScore
            tandemDupGap = self.phylDiagTandemGapMax
            phylDiagGapMax = self.phylDiagGapMax
            phylDiagPValue = self.phylDiagPValue
            phylDiagTruncationMax = self.phylDiagTruncationMax
            phylDiagFilter = self.phylDiagFilter
        score = myScore.myScore(repetPath + 'genes.%s.list.bz2',
                                repetPath + 'ancGenes.%s.list.bz2',
                                self.speciesTreePath,
                                tandemDupGap=tandemDupGap,
                                phylDiagGapMax=phylDiagGapMax,
                                phylDiagPValue=phylDiagPValue,
                                phylDiagTruncationMax=phylDiagTruncationMax,
                                phylDiagFilter=phylDiagFilter,
                                logErrPath=repetPath + ID + '.stderr',
                                verbose=self.scoreCalculationVerbose,
                                ID=ID)
        if saveScore:
            score.save(repetPath + ID + '.score')
        return score

    @staticmethod
    def changeAbsValuesIntoRates(specValues, speciesTree, decimalPlaces=6):
        # Warning: in specValues, values are all absolutes!
        # Even 'geneTandemDupProp' corresponds the absolute number of tandem duplications!!!
        specRates = {}
        for (branch, values) in specValues.iteritems():
            specRates[branch] = {}
            for (valueName, value) in values.iteritems():
                if valueName != 'geneTandemDupProp':
                    specRates[branch][valueName] = round(float(value) / speciesTree.parent[branch].distance,
                                                         decimalPlaces)
            if specValues[branch]['geneDup'] != 0:
                specRates[branch]['geneTandemDupProp'] = round(
                    float(specValues[branch]['geneTandemDupProp']) / specValues[branch]['geneDup'], decimalPlaces)
            else:
                specRates[branch]['geneTandemDupProp'] = 0.0
        return (specRates)

    def _changeMSInputsGridSearch(self, iSim):
        if self._gridPointNumber == 0:
            iteratorNumber = -1
        else:
            iteratorNumber = self._inputParameterGridPoints[self._gridPointNumber - 1]

        self._lastSpecRates = copy.deepcopy(self._originalSpecRates)
        self._lastParameterFiles = copy.deepcopy(self._originalParameterFiles)

        values = self._inputParameterGrid.next()
        iteratorNumber += 1
        while iteratorNumber != self._inputParameterGridPoints[self._gridPointNumber]:
            values = self._inputParameterGrid.next()
            iteratorNumber += 1
        for (i, value) in enumerate(values):
            if i < len(self.specRatesChangeMode[2]):
                # Case if it is a parameterFile parameter
                parameterName = self.specRatesChangeMode[2][i][0]
                if not isinstance(value, float):
                    self._lastParameterFiles[parameterName] = value
                else:
                    self._lastParameterFiles[parameterName] = self._originalParameterFiles[parameterName] * (1 + value)
            else:
                # Case if value is in specRates, e.g. chrFission
                specRatesComboNumber = i - len(self.specRatesChangeMode[2])
                speciesNameNumber = specRatesComboNumber / len(self.specRatesChangeMode[1])
                speciesName = self._originalSpecRates.keys()[speciesNameNumber]
                parameterNameNumber = specRatesComboNumber % len(self.specRatesChangeMode[1])
                parameterName = self.specRatesChangeMode[1][parameterNameNumber][0]
                self._lastSpecRates[speciesName][parameterName] = self._originalSpecRates[speciesName][
                                                                      parameterName] * (
                                                                      1 + value)
        self._gridPointNumber += 1

    def _changeMSInputsABC(self, iSim):
        self._lastSpecRates = copy.deepcopy(self._originalSpecRates)
        self._lastParameterFiles = copy.deepcopy(self._originalParameterFiles)
        for (i, value) in enumerate(self.specRatesChangeMode[2]):
            parameterName = self.specRatesChangeMode[2][i][0]
            boundLower = self.specRatesChangeMode[2][i][1][0]
            boundUpper = self.specRatesChangeMode[2][i][1][1]
            self._lastParameterFiles[parameterName] = random.uniform(boundLower, boundUpper)
            print >> sys.stderr, 'yupp'
            print >> sys.stderr, (parameterName, boundLower, boundUpper, self._lastParameterFiles[parameterName])

    def changeMSInputs(self, iSim, scoreCompArchive):
        if self.specRatesChangeMode[0] == 'constant':
            pass

        elif self.specRatesChangeMode[0] == 'randomwalk':
            # this function may change self.hasConverged
            self._lastSpecRates = self._changeMSInputStepwise(iSim, myOptimizer._calcSRNextValue, self._lastSpecRates,
                                                              scoreCompArchive=scoreCompArchive)
            # TODO
            self._lastParameterFiles = self._lastParameterFiles

        elif self.specRatesChangeMode[0] == 'grid':
            self._changeMSInputsGridSearch(iSim)

        # elif self.specRatesChangeMode[0].lower() == 'abc':
        #     self._changeMSInputsABC(iSim)

        else:
            raise ValueError('Wrong executionMode ("' + str(
                self.specRatesChangeMode[0]) + '") in myOptimizer.changeMSInput for Sim ' + str(iSim))

    def _changeMSInputStepwise(self, iSim, functionHowToChange, oldSpecRates, scoreCompArchive=None,
                               sigmaRatio=0.25):
        oldInputAbsValues = myOptimizer.changeRatesIntoAbsValues(oldSpecRates, self.speciesTree)

        nothingChanged = True  # Test if some value still develops; if not, member hasConverged = True
        newValues = {}
        for (branch, values) in oldInputAbsValues.iteritems():
            newValues[branch] = {}
            for (valueName, value) in values.iteritems():
                if valueName in self.specRatesChangeMode[1]:
                    newValues[branch][valueName] = oldInputAbsValues[branch][valueName]
                else:
                    newValues[branch][valueName] = functionHowToChange(oldInputAbsValues[branch][valueName],
                                                                       valueName,
                                                                       branch,
                                                                       iSim,
                                                                       self.scoreReal,
                                                                       scoreCompArchive=scoreCompArchive,
                                                                       sigmaRatio=sigmaRatio)
                    if newValues[branch][valueName] != oldInputAbsValues[branch][valueName]:
                        nothingChanged = False

            # There can never be more tandemDup than dup, hence limit it:
            newValues[branch]['geneTandemDupProp'] = min(newValues[branch]['geneTandemDupProp'],
                                                         newValues[branch]['geneDup'])

        if nothingChanged:
            print >> sys.stderr, 'Optimization converged !'
            self.hasConverged = True
        newSpecRates = myOptimizer.changeAbsValuesIntoRates(newValues, self.speciesTree)
        return (newSpecRates)

    @staticmethod
    def _calcSRNextValue(inputAbsValues, valueName, branch, iSim, scoreReal, scoreCompArchive=None, sigmaRatio=0):
        # usually Val : is the absolute number corresponding of the lastSpecRates
        # TODO: compute the diff between the given Val nb of events in specRates and the number of performed events
        # within magSimus
        SimID = str(iSim - 1) + '-real'
        scoreName = myOptimizer._calcMSNameToMCAName(valueName)

        # FIXME: its is an unelegant way to select only "inversions" and "translocations"
        if 'branches' in scoreName:
            # recBest = scoreCompArchive.calcBestScore(recScore, species)
            # if recBest['ID'] == SimID:
            # last Sim was best simulation

            # difference compared to the last input rates
            diff = scoreCompArchive.scores[scoreName][branch][SimID]['mean']
            print >> sys.stderr, "scoreCompArchive.scores[%s][%s][%s]['mean'] = %s" % (scoreName, branch, SimID, diff)
            diffRatio = scoreCompArchive.scores[scoreName + 'Ratio'][branch][SimID]['geoMean']
            if doNotUseDiffRatio:  # by default
                # remove the first 'branches-' at the beginning and remove '-Diff' at the end
                scoreName2 = '-'.join(scoreName.split('-')[1:-1])
                # input : input of magSimus
                # simOutput : after magSimus and M2006
                # realOutput : after M2006
                excessOutput = diff
                # realOutput is just used to compute SimOutput
                # remove the '-Diff' at the end of the scoreName
                realOutput = scoreReal.branchesEst[branch][scoreName2]
                simOutput = realOutput + excessOutput
                simInput = inputAbsValues
                print >> sys.stderr, 'valueName=%s' % valueName
                print >> sys.stderr, 'simOutput=%s' % simOutput
                print >> sys.stderr, 'realOutput=%s' % realOutput
                print >> sys.stderr, 'excessOutput=%s' % excessOutput
                print >> sys.stderr, 'oldSimInput=%s' % simInput
                if simOutput < realOutput:
                    assert excessOutput <= 0
                    keptFactor = float(simOutput) / float(simInput) if simInput > 0 else 1.0
                    if keptFactor != 0:
                        assert keptFactor > 0
                        excessInput = excessOutput / keptFactor
                    else:
                        assert keptFactor == 0
                        excessInput = 0
                    assert excessInput <= 0
                    newValue = max(0, simInput - excessInput)
                else:
                    # this sometimes happens due to the LSE part
                    # what we see is that the branches represent branches where the input is 0, but du to the LSE, some
                    # events are put there. In all the seen cases the simOutput is lowly higher than simInput, not more than 3 events.
                    print >> sys.stderr, "Rare case when there are more output rearrangements than input\n" +\
                                         " %s %s : input = %s > output = %s" % (branch, scoreName2, simInput, simOutput)
                    excessInput = excessOutput
                    assert 0.0 <= excessOutput

                    newValue = max(0, simInput - excessInput)
                print >> sys.stderr, 'newValue=%s' % newValue
            else:
                if np.isfinite(diffRatio):
                    newValue = max(0, inputAbsValues - diffRatio)
                else:
                    print >> sys.stderr, "Since diffRatio=%s is not a finite number, we use the value from the logErr of magSimus instead" % diffRatio
        else:
            newValue = inputAbsValues
        return newValue

    def _addNoiseSR(self, iSim, oldSpecRates, sigmaRatio=0.25):
        def functionHowToChange(oldVal, valueName, species, iSim, scoreCompArchive=None, sigmaRatio=0):
            out = max(0, np.random.normal(oldVal,
                                          max(1, oldVal * sigmaRatio)))
            return (out)

        # this function may change self.hasConverged
        newSpecRates = self._changeMSInputStepwise(iSim, functionHowToChange, oldSpecRates,
                                                   sigmaRatio=sigmaRatio)
        return (newSpecRates)

    @staticmethod
    def changeRatesIntoAbsValues(specRates, speciesTree):
        specValues = {}
        for (branch, rates) in specRates.iteritems():
            specValues[branch] = {}
            for (rateName, rate) in rates.iteritems():
                if rateName != 'geneTandemDupProp':
                    specValues[branch][rateName] = round(rate * speciesTree.parent[branch].distance, 0)
            specValues[branch]['geneTandemDupProp'] = round(
                specValues[branch]['geneDup'] * specRates[branch]['geneTandemDupProp'], 0)
        return (specValues)

    def cleanAfterCondor(self):
        os.remove(self.optimizationPath + 'execMSandScoreOnCondor.py')
        os.remove(self.optimizationPath + 'rec.optimizer')

    def cleanAfterLocalMode(self, iSim):
        for iRepet in xrange(1, self.nRepPerInput + 1):
            repetPath = self.optimizationPath + str(iSim) + '/' + str(iRepet) + '/'
            for fileName in os.listdir(repetPath):
                if len(fileName) >= 7 and fileName[-7:] == '.stdout':
                    os.remove(repetPath + fileName)
                elif len(fileName) >= 9 and fileName[-9:] == '.list.bz2':
                    os.remove(repetPath + fileName)

    def createMS1Arguments(self, iSim, iRepet):

        def getExactFileName(path, stringSnippet):
            fileName = [fileName for fileName in os.listdir(path) if stringSnippet in fileName]
            if len(fileName) == 0:
                raise ValueError('File "' + stringSnippet + '" not found in directory ' + path)
            else:
                return (fileName[0])

        simPath = self.optimizationPath + str(iSim) + '/'
        repetPath = simPath + str(iRepet) + '/'
        sys.argv = ['src/magSimus1.py',
                    self.speciesTreePath,
                    '-parameterFile=' + simPath + getExactFileName(simPath, 'parameter'),
                    '-userRatesFile=' + simPath + getExactFileName(simPath, 'specRate'),
                    '-out:genomeFiles=' + repetPath + 'genes.%s.list.bz2',
                    '-out:ancGenesFiles=' + repetPath + 'ancGenes.%s.list.bz2',
                    ('+breakpointAnalyzer' if withBA else '-breakpointAnalyzer'),
                    '-out:empiricalSbsAsGenomes=' + repetPath + 'sbs.genes.%s.%s.list.bz2'
                    ]
        (argsFormatList, optsFormatList) = magSimus1.getDefaultParameters()
        arguments = myTools.checkArgs(argsFormatList, optsFormatList, __doc__, showArgs=False)
        return (arguments)

    @staticmethod
    def createOptimizationDir(pathToSimulation):

        if pathToSimulation[-1] == '/':
            pathToSimulation = pathToSimulation[:-1]

        optimizationNumber = ''
        while os.path.exists(pathToSimulation + optimizationNumber):
            if optimizationNumber == '':
                optimizationNumber = str(1)
            else:
                optimizationNumber = str(int(optimizationNumber) + 1)
        optimizationPath = pathToSimulation + optimizationNumber + '/'
        os.makedirs(optimizationPath)
        return (optimizationPath)

    def fillRealScoreDataAndStoreThemOnHDD(self, speciesTreePath, realScorePath):
        self.realDataPath = self.optimizationPath + 'realData/'
        os.makedirs(self.realDataPath)

        self.speciesTreePath = self.realDataPath + 'speciesTree.phylTree'

        shutil.copy(speciesTreePath, self.realDataPath)

        if realScorePath == '':
            print >> sys.stderr, 'Create real score.'
            self.scoreReal = myScore.myScore(self.realGenomes, self.realAncGenomes,
                                             self.speciesTreePath,
                                             tandemDupGap=self.phylDiagTandemGapMax,
                                             phylDiagGapMax=self.phylDiagGapMax,
                                             phylDiagPValue=self.phylDiagPValue,
                                             phylDiagTruncationMax=self.phylDiagTruncationMax,
                                             phylDiagFilter=self.phylDiagFilter,
                                             verbose=self.scoreCalculationVerbose,
                                             ID='real')
            print >> sys.stderr, 'Score creation finished.'
        else:
            print >> sys.stderr, 'Load real score.'
            self.scoreReal = myScore.myScore.load(realScorePath)
            print >> sys.stderr, 'Loading finished.'
        self.scoreReal.save(self.realDataPath + 'real.score')
        print >> sys.stderr, 'Real score saved as %s' % self.realDataPath + 'real.score'

    def doRepetition(self, optPath, iSim, iRepet):
        repetPath = optPath + str(iRepet) + '/'
        os.makedirs(repetPath)
        stopTimes = [time.time()]
        self.launchMagSimus1(iSim, iRepet)
        stopTimes.append(time.time())
        score = self.calculateScore(iSim, iRepet)
        name = self.optimizationPath + str(iSim) + '/' + str(iRepet) + '/' + str(iSim) + '-' + str(iRepet)
        if saveScore:
            score.save(name + '.score')
        stopTimes.append(time.time())
        comp = myScore.myScoreComparison(score, self.scoreReal)
        stopTimes.append(time.time())
        comp.save(name + '.comp')
        stopTimes.append(time.time())
        print >> sys.stderr, 'RunTimes in s: MS = {0:.1f}, Score = {1:.1f}, ScoreComp = {2:.1f}, Save = {3:.1f}'.format(
            stopTimes[1] - stopTimes[0], stopTimes[2] - stopTimes[1], stopTimes[3] - stopTimes[2],
            stopTimes[4] - stopTimes[3])
        return (comp)

    # def executeCondorModeOld(self):
    #     (startNumber, endNumber) = (0, self.numberOfSimulations)
    #
    #     # ####################################
    #     # 1) Prepare Optimization
    #     # Write inputs for Condor execution to hard drive: the script to be executed and the pickled myOptimizerInstance
    #     self.writeCondorExecScript(self.optimizationPath)
    #     self.save(self.optimizationPath + 'rec.optimizer')
    #
    #     scoreCompArchive = myScore.myScoreComparisonArchive()
    #     self.prepareOptimization(startNumber, endNumber, scoreCompArchive)
    #     print >> sys.stderr, "Jobs are send to Condor."
    #     # ####################################
    #     # 2) Submit To Condor
    #     listOfJobIDs = []
    #     simRepCombos = itertools.product(xrange(startNumber, endNumber + 1), xrange(1, self.nRepPerInput + 1))
    #     for (iSim, iRepet) in simRepCombos:
    #         jobID = myCondor.submit(self.optimizationPath + 'execMSandScoreOnCondor.py',
    #                                 arguments = self.optimizationPath + ' ' + str(iSim) + ' ' + str(iRepet),
    #                                 niceUser = True, universe = "vanilla", log = "condorpy.log",
    #                                 outfile = "condorpy.stdout.%s.log" % "$(Cluster)",
    #                                 errfile = "condorpy.stderr.%s.log" % "$(Cluster)")
    #         listOfJobIDs.append((iSim, iRepet, jobID))
    #         if (iSim*self.nRepPerInput + iRepet) % 1000 == 1:
    #             print >> sys.stderr, "Simulation {}-{} sent to Condor (jobID = {})".format(iSim, iRepet, jobID)
    #
    #     # ####################################
    #     # 3) Wait Until Condor Has Finished
    #     for (iSim, iRepet, jobID) in listOfJobIDs:
    #         stdout, stderr = myCondor.getoutput(jobID)
    #         recRepetPath = self.optimizationPath + '/' + str(iSim) + '/' + str(iRepet) + '/' + str(iSim) + '-' + str(iRepet)
    #         with open(recRepetPath + '.condorstdout', "w") as text_file:
    #             text_file.write(stdout)
    #         with open(recRepetPath + '.condorstderr', "w") as text_file:
    #             text_file.write(stderr)
    #         if (iSim * self.nRepPerInput + iRepet) % 1000 == 1:
    #             print >> sys.stderr, "Simulation {}-{} finished.".format(iSim, iRepet)
    #
    #     self.cleanAfterCondor()

    def executeCondorMode(self):

        # ####################################
        # 1) Prepare Optimization
        # Write inputs for Condor execution to hard drive: the script to be executed and the pickled myOptimizerInstance
        self.writeCondorExecScript(self.optimizationPath, optimised=True)
        self.save(self.optimizationPath + 'rec.optimizer')

        # DEBUG
        print >> sys.stderr, "self.specRatesChangeMode = %s" % str(self.specRatesChangeMode)

        if (self.specRatesChangeMode[0].lower() == 'grid'
            or self.specRatesChangeMode[0].lower() == 'constant'
            or self.specRatesChangeMode[0].lower() == 'abc'):
            (startNumber, endNumber) = (0, self.numberOfSimulations)
        else:
            (startNumber, endNumber) = (0, 0)

        # DEBUG
        print >> sys.stderr, "self.numberOfSimulations = %s" % self.numberOfSimulations
        print >> sys.stderr, "(startNumber, endNumber) = %s" % str((startNumber, endNumber))

        scoreCompArchive = myScore.myScoreComparisonArchive(withBA=withBA)

        print >> sys.stderr, "Beginning: ", startNumber, endNumber, self.numberOfSimulations
        while endNumber <= self.numberOfSimulations:
            # DEBUG
            print >> sys.stderr, startNumber, endNumber, self.numberOfSimulations

            if self.hasConverged:
                break

            # this function changes self.hasConverged
            self.prepareOptimization(startNumber, endNumber, scoreCompArchive, self.scoreReal)

            print >> sys.stderr, "Jobs are sent to Condor."
            simRepCombos = list(itertools.product(xrange(startNumber, endNumber + 1), xrange(1, self.nRepPerInput + 1)))
            executable = self.optimizationPath + 'execMSandScoreOnCondor.py'
            # listOfArguments = [self.optimizationPath + ' ' + str(iSim) + ' ' + str(iRepet) for (iSim, iRepet) in simRepCombos]

            # Ensure that jobs are sent by packs of at most 3000 jobs
            nJobsFinished = 0
            nJobsSentAtOnce = 2000
            cs = myCondor.CondorSubmitor('execMSandScoreOnCondor')
            while nJobsFinished < len(simRepCombos):
                simRepCombos = simRepCombos[nJobsFinished:(nJobsFinished + nJobsSentAtOnce)]
                (iSim1st, iRepet1st) = simRepCombos[0]
                (iSimLast, iRepetLast) = simRepCombos[-1]
                print >> sys.stderr, "A bunch of simulations " + self.optimizationPath + " from " + '/'.join(
                    [str(iSim1st), str(iRepet1st)]) + ' to ' + '/'.join(
                    [str(iSimLast), str(iRepetLast)]) + " are processed."
                listOfArguments = [self.optimizationPath + ' ' + str(iSim) + ' ' + str(iRepet) for (iSim, iRepet) in
                                   simRepCombos]

                #######################
                # 1) Submit To Condor #
                #######################
                listOfJobIds = cs.submit_ManyJobs(executable, listOfArguments, niceUser=True,
                                                  maxSimultaneousJobsInGroup=None)
                simRepetOfJobId = {}
                for idx, jobid in enumerate(listOfJobIds):
                    simRepetOfJobId[jobid] = simRepCombos[idx]

                ##############################
                # 2) get results from Condor #
                ##############################
                for (jobid, stdoutFileName, stderrFileName) in cs.getoutput_ManyJobs(listOfJobIds):
                    (iSim, iRepet) = simRepetOfJobId[jobid]
                    repetPath = self.optimizationPath + str(iSim) + '/' + str(iRepet) + '/' + str(iSim) + '-' + str(
                        iRepet)
                    try:
                        shutil.move(stdoutFileName, repetPath + '.condorstdout')
                        shutil.move(stderrFileName, repetPath + '.condorstderr')
                    except:
                        print >> sys.stderr, 'Problem with Condor results for JobID: ' + str(jobid)

                nJobsFinished += nJobsSentAtOnce

            for iSim in xrange(startNumber, endNumber + 1):
                for iRepet in xrange(1, self.nRepPerInput + 1):
                    # DEBUG
                    print >> sys.stderr, "iSim=%s and iRepet=%s" % (iSim, iRepet)
                    pathToRes = self.optimizationPath + str(iSim) + '/' + str(iRepet) + '/' + str(iSim) + '-' + str(
                        iRepet) + '.comp'
                    if os.path.isfile(pathToRes):
                        comp = myScore.myScoreComparison.load(pathToRes)
                        assert isinstance(comp, myScore.myScoreComparison)
                        # the score comp.ID = iSim-iRepet-real
                        #scoreCompArchive.addComp(comp, namePrefix=self.optimizationPath.split('/')[-3:])
                        scoreCompArchive.addComp(comp)
                    else:
                        print >> sys.stderr, "Warning, simulation+score-comp: %s/%s may not have finished successfully, the path does not exist" % (
                        iSim, iRepet)
                        print >> sys.stderr, "cannot find file %s" % pathToRes
            scoreCompArchive.calcStatistics()

            startNumber += 1
            endNumber += 1

        print >> sys.stderr, 'Condor jobs finished !'
        self.cleanAfterCondor()
        return scoreCompArchive

    def executeLocalMode(self):

        (startNumber, endNumber) = (0, self.numberOfSimulations)
        scoreCompArchive = myScore.myScoreComparisonArchive()

        for iSim in xrange(startNumber, endNumber + 1):
            print >> sys.stderr, 'Starting simulation ' + str(iSim)
            optPath = self.optimizationPath + str(iSim) + '/'
            os.makedirs(optPath)
            self.writeSpecRatesFile(self._lastSpecRates, self.speciesTree,
                                    optPath + self.specRatesBaseName + '-' + str(iSim))
            self.writeParameterFile(self._lastParameterFiles, optPath + self.parameterFileBaseName + '-' + str(iSim))

            if self.nCPU > 1:
                newPool = multiprocessing.Pool(processes=self.nCPU)
                doRepInput = [(self, optPath, iSim, iRepet) for iRepet in xrange(1, self.nRepPerInput + 1)]
                compList = newPool.map(doRepetitionParallel, doRepInput)

                for comp in compList:
                    scoreCompArchive.addComp(comp)
            else:
                for iRepet in xrange(1, self.nRepPerInput + 1):
                    # Do all action for one repetition and return comparison to real score
                    comp = self.doRepetition(optPath, iSim, iRepet)
                    scoreCompArchive.addComp(comp)

            # DEBUG ###################
            #print >> sys.stderr, 'scoreCompBeforeCalcStatsSavedIn %s' % '/home/jlucas/test'
            #scoreCompArchive.save('/home/jlucas/test')
            #scoreCompArchive = myScore.myScoreComparisonArchive.load('/home/jlucas/test')
            scoreCompArchive.calcStatistics()
            # /DEBUG ###################
            if iSim != endNumber:
                self.changeMSInputs(iSim + 1, scoreCompArchive)

            if self.deleteMagSimusOutput:
                self.cleanAfterLocalMode(iSim)

            if self.hasConverged:
                break

        return (scoreCompArchive)

    def getSpecTotals(self, specValues):
        totals = {}
        for (branch, rates) in specValues.iteritems():
            for (rateName, rate) in rates.iteritems():
                if rateName not in totals:
                    totals[rateName] = rate
                else:
                    totals[rateName] += rate
        return (totals)

    def getSpecStandardDeviation(self, specValues):
        totalValues = self.getSpecTotals(specValues)
        nObservations = len(specValues)
        meanValues = {valueName: float(totalValue) / nObservations for (valueName, totalValue) in
                      totalValues.iteritems()}

        sds = {}
        for (branch, values) in specValues.iteritems():
            for (valueName, value) in values.iteritems():
                if valueName not in sds:
                    sds[valueName] = float(value - meanValues[valueName]) ** 2
                else:
                    sds[valueName] += float(value - meanValues[valueName]) ** 2
        sds = {valueName: np.sqrt(sd / (nObservations - 1)) for (valueName, sd) in sds.iteritems()}
        return (sds)

    def launchMagSimus1(self, iSim, iRepet):

        ID = str(iSim) + '-' + str(iRepet)
        repetPath = self.optimizationPath + str(iSim) + '/' + str(iRepet) + '/'

        # Direct StdOut and StdErr to files
        stderrOld = sys.stderr
        stderrFile = open(repetPath + ID + '.stderr', 'wb')
        sys.stderr = stderrFile

        stdoutOld = sys.stdout
        # stdoutFile = open(repetPath + ID + '.stdout', 'wb')
        stdoutFile = open(os.devnull, 'wb')
        sys.stdout = stdoutFile
        msParameters = self.createMS1Arguments(iSim, iRepet)
        try:
            magSimus1.launch(msParameters)
        except:
            exception = traceback.format_exc()
            print >> sys.stderr, 'MagSimus throwed an Exception in ' + ID
            print >> sys.stderr, exception

        stderrFile.close()
        sys.stderr = stderrOld

        stdoutFile.close()
        sys.stdout = stdoutOld

    @staticmethod
    def load(fileName):
        with open(fileName, 'rb') as inFile:
            optimizer = pickle.load(inFile)
        return (optimizer)

    def prepareGridSearch(self, specRatesChangeMode):
        specRatesKeys = self._lastSpecRates.itervalues().next().keys()
        parameterFileKeys = self._lastParameterFiles.keys()

        # DEBUG
        print >> sys.stderr, "len(specRatesChangeMode) = %s" % len(specRatesChangeMode)
        print >> sys.stderr, "specRatesChangeMode = %s" % str(specRatesChangeMode)

        # Check if the input was correct
        if not (len(specRatesChangeMode) == 4
                and isinstance(specRatesChangeMode[1], list)
                and isinstance(specRatesChangeMode[2], list)
                and all([vtuple[0] in specRatesKeys for vtuple in specRatesChangeMode[1]])
                and all([vtuple[0] in parameterFileKeys for vtuple in specRatesChangeMode[2]])
                and all([isinstance(vtuple[1], list) for vtuple in specRatesChangeMode[1]])
                and all([isinstance(vtuple[1], list) for vtuple in specRatesChangeMode[2]])
                and specRatesChangeMode[3] > 0):
            # Case if specRatesChangeMode is somehow miss-specified, in this case there is an exception thrown
            # due to self.specRatesChangeMode = None
            return

        self.specRatesChangeMode = specRatesChangeMode

        nBranches = len(self._originalSpecRates)
        gridPointsPerBranch = [vtuple[1] for vtuple in self.specRatesChangeMode[2]] + nBranches * [vtuple[1] for vtuple
                                                                                                   in
                                                                                                   self.specRatesChangeMode[
                                                                                                       1]]
        gridSizeTotal = np.prod([len(gridPoints) for gridPoints in gridPointsPerBranch])
        self._inputParameterGrid = itertools.product(*gridPointsPerBranch)

        # We should have at least 1 point sampled and at max all the points be sampled
        sampleSize = min(int(np.ceil(gridSizeTotal * self.specRatesChangeMode[3])), gridSizeTotal)
        self.numberOfSimulations = sampleSize
        print >> sys.stderr, 'Grid is created, and number of simulations is automatically set to ' + str(
            self.numberOfSimulations) + ' ( = ' + str(self.numberOfSimulations * self.nRepPerInput) + ' repetitions)'
        self._inputParameterGridPoints = sorted(random.sample(range(gridSizeTotal), sampleSize))
        # _gridPointNumbe indicates which number of the sample was taken last
        self._gridPointNumber = 0

    def parseSpecRatesChangeMode(self, specRatesChangeMode):
        # For convenience, one can enter the specRatesChangeMode as a string, an a default tuple is created here
        # FIXME: please choose one way to write constant
        if (specRatesChangeMode.lower() == 'constant'
            or specRatesChangeMode.lower() == "('constant')"
            or specRatesChangeMode.lower() == '("constant")'):
            self.specRatesChangeMode = ('constant',)
        elif specRatesChangeMode.lower() == 'randomwalk':
            self.specRatesChangeMode = (specRatesChangeMode, [])
        else:
            try:
                specRatesChangeMode = literal_eval(specRatesChangeMode)
            except:
                pass
            if len(specRatesChangeMode) > 0:
                if ((specRatesChangeMode[0].lower() == 'constant' and len(specRatesChangeMode) == 1)
                    or (specRatesChangeMode[0].lower() == 'randomwalk' and len(specRatesChangeMode) == 2)
                    or (specRatesChangeMode[0].lower() == 'abc' and len(specRatesChangeMode) == 3)):
                    self.specRatesChangeMode = specRatesChangeMode
                elif (specRatesChangeMode[0].lower() == 'grid'):
                    self.prepareGridSearch(specRatesChangeMode)

        if self.specRatesChangeMode is None:
            # Case if specRatesChangeMode were not recognized by decision tree above and hence is still at default
            raise ValueError('specRatesChangeMode not well specified !')

    def prepareOptimization(self, startNumber, endNumber, scoreCompArchive, scoreReal):

        for iSim in xrange(startNumber, endNumber + 1):
            if startNumber > 0:
                # this function changes self.hasConverged
                self.changeMSInputs(iSim, scoreCompArchive)

            optPath = self.optimizationPath + str(iSim) + '/'
            os.makedirs(optPath)
            self.writeSpecRatesFile(self._lastSpecRates, self.speciesTree,
                                    optPath + self.specRatesBaseName + '-' + str(iSim))
            self.writeParameterFile(self._lastParameterFiles,
                                    optPath + self.parameterFileBaseName + '-' + str(iSim))

            for iRepet in xrange(1, self.nRepPerInput + 1):
                repPath = optPath + str(iRepet) + '/'
                os.makedirs(repPath)

    @staticmethod
    def _calcMSNameToMCAName(MSName):
        """
        change name from names used in magSimus to names used in myScoreCompsArchive to score this name
        :param MSName: the name in magSimus
        :return: the name in myScoreCompsArchive
        """
        if MSName == 'geneTandemDupProp':
            # important tu return an empty ()
            return ()
        elif MSName == 'chrFusion':
            return 'nChr-Diff'
        elif MSName == 'geneDup':
            return 'geneEventRatesOnBranches-geneDup-LogRatio'
        elif MSName == 'geneLoss':
            return 'geneEventRatesOnBranches-geneLoss-LogRatio'
        elif MSName == 'chrInvert':
            return ('branches-nInv-Diff')
        elif MSName == 'chrFission':
            return ('nChr-Diff')
        elif MSName == 'chrTransloc':
            return ('branches-nTransl-Diff')
        elif MSName == 'geneBirth':
            # FIXME (from Joseph): why do geneBirths are associated with geneDups ?
            return ('geneEventRatesOnBranches-geneDup-LogRatio')
        else:
            raise ValueError('Not defined MagSimus variable name !')

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
                elif ((x[0][:3] == 'chr' and x[0][:22] != 'chr:invDist' and x[0][
                                                                            :29] != 'chr:absOrRelChoiceOfInvLength')
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

        return (inDict)

    @staticmethod
    def readSpecRates(inPath):
        userRatesReader = myFile.myTSV.readTabular(inPath, [str, str, str, float])

        inDict = {}
        for x in userRatesReader:
            (rateName, parentSp, childSp, value) = x
            if childSp not in inDict:
                inDict[childSp] = {rateName: value}
            else:
                inDict[childSp][rateName] = value
        return (inDict)

    def save(self, fileName=None):
        if fileName is None:
            fileName = 'rec.optimizer'
        with open(fileName, 'wb') as outFile:
            pickle.dump(self, outFile, -1)

    # # This function is not used at the moment and may be erased at some point
    # # TODO erase this function and use myCondor API for sending jobs to condor
    # def writeCondorInstructionsForMagSimus(self, fileName, startNumber, endNumber):
    #     if self.magSimusPath.split('.')[-1] == 2:
    #         print >> sys.stderr, 'writeCondorInstructions works only for magSimus1 right now!'
    #         return (None)
    #
    #     out = '#\nExecutable = ' + self.magSimusPath
    #     out += '\nUniverse = vanilla\nOutput = ' + str(startNumber) + '/' + '1' + '/logStd'
    #     out += '\nInput =\nError = ' + str(startNumber) + '/' + '1' + '/logErr'
    #     out += '\nLog = logCondor'
    #     out += '\nArguments = ' + self.speciesTreePath
    #     out += ' -userRatesFile = ' + str(startNumber) + '/' + self.specRatesBaseName + '-' + '1'
    #     out += ' -parameterFile = ' + str(startNumber) + '/' + self.parameterFileBaseName + '-' + '1'
    #     out += ' -out:genomeFiles = ' + str(startNumber) + '/' + '1' + '/' + 'genes.%s.list.bz2'
    #     out += ' -out:ancGenesFiles = ' + str(startNumber) + '/' + '1' + '/' + 'ancGenes.%s.list.bz2'
    #     out += '\nGetEnv = True\nInitialdir = ' + self.optimizationPath
    #     out += '\nshould_transfer_files = NO'
    #     out += '\nrun_as_owner = True'
    #     out += '\nRequirements = ((machine != "bioclust08.bioclust.biologie.ens.fr") & & (machine != "bioclust07.bioclust.biologie.ens.fr"))'
    #     out += '\nNotify_user = tittmann@biologie.ens.fr'
    #     out += '\nNotification = never'
    #     out += '\nNiceUser = False'
    #     out += '\nPriority ='
    #     out += '\nRank = kflops + 1000 * Memory'
    #     out += '\nQueue'
    #
    #     simRepCombos = itertools.product(xrange(startNumber, endNumber + 1), xrange(1, self.nRepPerInput + 1))
    #     for (iSim, iRepet) in simRepCombos:
    #         if iSim == 1 and iRepet == 1:
    #             continue
    #             # The first job was already written above
    #         out += '\n\nOutput = ' + str(iSim) + '/' + str(iRepet) + '/logStd'
    #         out += '\nError = ' + str(iSim) + '/' + str(iRepet) + '/ logErr'
    #         out += '\nArguments = ' + self.speciesTreePath
    #         out += ' -userRatesFile = ' + str(iSim) + '/' + self.specRatesBaseName + '-' + str(iSim)
    #         out += ' -parameterFile = ' + str(iSim) + '/' + self.parameterFileBaseName + '-' + str(iSim)
    #         out += ' -out:genomeFiles = ' + str(iSim) + '/' + str(iRepet) + '/' + 'genes.%s.list.bz2'
    #         out += ' -out:ancGenesFiles = ' + str(iSim) + '/' + str(iRepet) + '/' + 'ancGenes.%s.list.bz2'
    #         out += '\nQueue'
    #
    #     with open(self.optimizationPath + fileName, 'a') as outFile:
    #         outFile.write(out)

    @staticmethod
    # http://stackoverflow.com/questions/3306518/cannot-pass-an-argument-to-python-with-usr-bin-env-python
    # http://stackoverflow.com/questions/17458528/why-does-this-snippet-with-a-shebang-bin-sh-and-exec-python-inside-4-single-q
    def writeCondorExecScript(pathToDir, optimised=False):
        if optimised is True:
            condorScript = "#!/bin/sh\n"
            condorScript += "\'\'\'\'exec python -O -- \"$0\" ${1+\"$@\"} # \'\'\'\n"  # the option -O is for optimised!
        else:
            condorScript = "#!/usr/bin/python\n"
        condorScript += "\n"
        condorScript += "import analysis.myScore as myScore\n"
        condorScript += "import os\n"
        condorScript += "# Cannot use import myOptimizer as ... because pickle starts to bug\n"
        condorScript += "from analysis.myOptimizer import myOptimizer\n"
        condorScript += "import utils.myTools as myTools\n"
        condorScript += "\n"
        condorScript += "__doc__ = ''\n"
        condorScript += "\n"
        condorScript += "arguments = myTools.checkArgs([('optimizationPath', str), ('iSim', int), ('iRepet', int)], [], __doc__, showArgs=False)\n"
        condorScript += "\n"
        condorScript += "optimizer = myOptimizer.load(arguments['optimizationPath'] + 'rec.optimizer')\n"
        condorScript += "optimizer.launchMagSimus1(arguments['iSim'], arguments['iRepet'])\n"
        # score.ID = arguments['iSim'] + '-' + arguments['iRepet']
        condorScript += "score = optimizer.calculateScore(arguments['iSim'], arguments['iRepet'], saveScore=False)\n"
        condorScript += "name = arguments['optimizationPath'] + str(arguments['iSim']) + '/' + str(arguments['iRepet']) + '/' + str(arguments['iSim']) + '-' + str(arguments['iRepet'])\n"
        if saveScore:
            condorScript += "score.save(name + '.score')\n"
        # optimizer.scoreReal.ID = 'real'
        # comp.ID = score.ID + '-' + optimizer.scoreReal.ID
        condorScript += "comp = myScore.myScoreComparison(score, optimizer.scoreReal)\n"
        condorScript += "comp.save(name + '.comp')\n"
        condorScript += "if optimizer.deleteMagSimusOutput:\n"
        condorScript += "   repetPath = optimizer.optimizationPath + str(arguments['iSim']) + '/' + str(arguments['iRepet']) + '/'\n"
        condorScript += "   for fileName in os.listdir(repetPath):\n"
        condorScript += "       if len(fileName) >= 7 and fileName[-7:] == '.stdout':\n"
        condorScript += "           os.remove(repetPath + fileName)\n"
        condorScript += "       elif len(fileName) >= 9 and fileName[-9:] == '.list.bz2':\n"
        condorScript += "           os.remove(repetPath + fileName)\n"
        execScriptName = pathToDir + 'execMSandScoreOnCondor.py'
        with open(execScriptName, "w") as text_file:
            text_file.write(condorScript)

        if os.name == 'posix':
            os.system('chmod a+x ' + execScriptName)

    @staticmethod
    def writeParameterFile(outputDict, outPath):
        outputList = []
        for (parameter, value) in outputDict.iteritems():
            outputList.append(parameter + '\t' + str(value))
        outputList.sort()
        output = '\n'.join(outputList) + '\n'
        with open(outPath, "w") as text_file:
            text_file.write(output)

    @staticmethod
    def writeSpecRatesFile(outputDict, speciesTree, outPath):
        output = ''
        for (branch, rates) in outputDict.iteritems():
            for (rateName, rate) in rates.iteritems():
                output += rateName + '\t' + speciesTree.parent[branch].name + '\t' + branch + '\t' + str(rate) + '\n'

        with open(outPath, "w") as text_file:
            text_file.write(output)

    def launch(self, magSimusPath, chrNumbers, startSpecRatesPath='', startParameterFilePath='',
               numberOfSimulations=100, nRepPerInput=5, executionMode='',
               specRatesChangeMode="('constant',)", deleteMagSimusOutput=False):
        # executionMode = ['', 'parallel', 'condor', int()]
        # '' launches simulation on 1 core
        # 'parallel' launches simulation on recent PC on multiple cores = min(1, numberOfCores - 1)
        # if executionMode is an integer, the optimization is launched parallel on that many cores
        # ATTENTION: parallel doesn't work in the interactive console of PyCharm
        # 'condor' uses condor for simulation calculation

        startTime = time.time()
        print >> sys.stderr, 'Optimization started !'

        self.magSimusPath = magSimusPath
        self.chrNumber = chrNumbers

        self.numberOfSimulations = numberOfSimulations
        self.nRepPerInput = nRepPerInput

        self.specRatesBaseName = startSpecRatesPath.split('/')[-1]
        self.parameterFileBaseName = startParameterFilePath.split('/')[-1]

        self.deleteMagSimusOutput = deleteMagSimusOutput

        if executionMode.lower() == 'parallel':
            self.nCPU = max(1, multiprocessing.cpu_count() - 1)
        elif executionMode.isdigit():
            self.nCPU = int(executionMode)

        if startSpecRatesPath != '' and startParameterFilePath != '':
            self._originalSpecRates = self.readSpecRates(startSpecRatesPath)
            self._originalParameterFiles = self.readParameterFile(startParameterFilePath)
        else:
            self._originalSpecRates = self.calcStartRatesFromScore(self.scoreReal, self.speciesTree)
            self._originalParameterFiles = self.readParameterFile(startParameterFilePath)

        self._lastSpecRates = copy.deepcopy(self._originalSpecRates)
        self._lastParameterFiles = copy.deepcopy(self._originalParameterFiles)

        self.parseSpecRatesChangeMode(specRatesChangeMode)
        print >> sys.stderr, self.specRatesChangeMode

        if executionMode.lower() == 'condor':
            scoreCompArchive = self.executeCondorMode()
            # self.executeCondorModeOld()
        else:
            scoreCompArchive = self.executeLocalMode()

        scoreCompArchive.calcStatistics()
        scoreCompArchive.save(self.optimizationPath + 'scoreComp.archive')
        print >> sys.stderr, 'scoreCompArchive recorded in %s' % self.optimizationPath + 'scoreComp.archive'
        del scoreCompArchive

        print >> sys.stderr, 'Optimization finished after {0}s!'.format(round(time.time() - startTime, 1))

    def __init__(self, pathToSimulation, realGenomesPath, realAncGenomesPath, speciesTreePath, realScorePath='',
                 phylDiagTandemGapMax=5,
                 phylDiagGapMax=5,
                 phylDiagPValue=None,
                 phylDiagTruncationMax=5,
                 phylDiagFilter=myDiags.FilterType.InBothGenomes,
                 scoreCalculationVerbose=False):
        self.speciesTree = myPhylTree.PhylogeneticTree(speciesTreePath)
        self.speciesTreePath = ''
        self.optimizationPath = ''
        self.realDataPath = ''
        self.magSimusPath = ''
        self.scoreReal = ''
        self.realAncGenomes = realAncGenomesPath
        self.realGenomes = realGenomesPath

        self.chrNumber = {}

        self.numberOfSimulations = 0
        self.nRepPerInput = 0

        self.phylDiagTandemGapMax = phylDiagTandemGapMax
        self.phylDiagGapMax = phylDiagGapMax
        self.phylDiagTruncationMax = phylDiagTruncationMax
        self.phylDiagPValue = phylDiagPValue
        self.phylDiagFilter = phylDiagFilter

        self.deleteMagSimusOutput = False
        self.scoreCalculationVerbose = scoreCalculationVerbose

        self.hasConverged = False

        # Internal variables
        # The next two are used as storage for version number of start specRates file,
        # e.g. may store "specRates.v80"
        self.specRatesBaseName = ''
        self.parameterFileBaseName = ''

        # How shall the old specRate / parameterFile be changed for the new wave of simulations
        # Until now, there are 3 different methods implemented:
        # the input is a tuple, with the first entry being the mode of specRate change, and the next ones are additional parameters for the mode
        # ('constant',) -> specRates stay the same over all simulations
        # ('randomwalk', ['chrFusion', 'chrFission'])
        #       The next specRates are based on a normally distributed random variation on the last specRates, while 'chrFusion'/'chrFission' are held constant
        # ('grid', list1, list2, float) -> This creates a mesh of specRateFiles
        #              list1 contains the information for the specRates in form of tuple,
        #                   the first entry in the tuple is the name of the parameter which shall be explored (e.g. ['chrFission', 'geneDup'],
        #                   the second entry are the grid points which shall be explored
        #                       two things are possible: 1) for an input which is a number, e.g. chrFission == 16.0
        #                                                   than the list should consist of floats which indicate percental changes,
        #                                                   e.g. ('chrFission', [-0.5, -0.25, 0, 0.25, 0.5]) -> this creates a grid with 8, 12, 16, 20, 24
        #                                               2) for an input which is a str or bool, e.g. forceClustering = True
        #                                                   than the list should consist of discrete points,
        #                                                   e.g. ('forceClustering', [True, False])
        #              list2 works like dict1, but is used for the parameterFile
        #              float should be a number 0 < x <= 1 which indicates how much of the grid should be sampled (1 for complete)
        self.specRatesChangeMode = None

        # Number of CPUs used; gets set with myOptimizer.launch()
        self.nCPU = 1

        self.optimizationPath = self.createOptimizationDir(pathToSimulation)
        print >> sys.stderr, 'Created optimization directory at ' + self.optimizationPath
        self.fillRealScoreDataAndStoreThemOnHDD(speciesTreePath, realScorePath)

        print >> sys.stderr, 'myOptimizer successfully built.'


class myABC:
    def __init__(self):
        pass

    @staticmethod
    def startABC(arguments):
        print >> sys.stderr, 'ABC mode launched!'
        # create super directory
        abcDir = arguments["pathToSimulation"] + '_' + str(arguments['randomSeed']) + '/'
        abcDir = myOptimizer.createOptimizationDir(abcDir)
        arguments["pathToSimulation"] = abcDir + 'ABCOpt'

        specRatesChangeMode = literal_eval(arguments["specRatesChangeMode"])

        arguments["specRatesChangeMode"] = "('randomwalk', ['chrFusion','chrFission'])"

        # change parameter file
        baseParaFileName = arguments['startParameterFilePath']
        paraFile = myOptimizer.readParameterFile(baseParaFileName)

        for iSeed in xrange(specRatesChangeMode[3]):
            for parameters in specRatesChangeMode[2]:
                if isinstance(parameters[1][0], int) or isinstance(parameters[1][0], float):
                    lowerBound = parameters[1][0]
                else:
                    boundList = parameters[1][0].split('*')
                    lowerBound = paraFile[boundList[0]]
                    if (len(boundList) > 1):
                        lowerBound *= double(boundList[1])

                if isinstance(parameters[1][1], int) or isinstance(parameters[1][1], float):
                    upperBound = parameters[1][1]
                else:
                    boundList = parameters[1][1].split('*')
                    upperBound = paraFile[boundList[0]]
                    if (len(boundList) > 1):
                        upperBound *= double(boundList[1])

                paraFile[parameters[0]] = random.uniform(lowerBound, upperBound)
                print >> sys.stderr, 'Bounds: ' + str(lowerBound) + ', ' + str(upperBound) + ' -> Gamma(1, ' + str(
                    paraFile[parameters[0]]) + ')'

            newParaFilePath = abcDir + '_ABC.'.join(baseParaFileName.split('/')[-1].split('.'))
            myOptimizer.writeParameterFile(paraFile, newParaFilePath)
            arguments['startParameterFilePath'] = newParaFilePath

            newOpt = myOptimizer(arguments["pathToSimulation"], arguments["genesFile"],
                                 arguments["ancGenesFile"], arguments["speciesTreePath"],
                                 realScorePath=arguments["realScorePath"],
                                 scoreCalculationVerbose=arguments["scoreCalculationVerbose"],
                                 )

            newOpt.launch(arguments["magsimus1Path"],
                          chrNumbers,
                          startSpecRatesPath=arguments["startSpecRatesPath"],
                          startParameterFilePath=arguments["startParameterFilePath"],
                          numberOfSimulations=arguments["numberOfSimulations"],
                          nRepPerInput=arguments["nRepPerInput"],
                          executionMode=arguments["executionMode"],
                          specRatesChangeMode=arguments["specRatesChangeMode"],
                          deleteMagSimusOutput=arguments["deleteMagSimusOutput"]
                          )
            del newOpt

        os.remove(newParaFilePath)


class myGrid:
    def __init__(self):
        pass

    @staticmethod
    def startGrid(arguments):
        print >> sys.stderr, 'Grid mode with optimisation launched!'
        # create super directory
        gridDir = arguments["pathToSimulation"] + '_' + str(arguments['randomSeed']) + '/'
        gridDir = myOptimizer.createOptimizationDir(gridDir)
        arguments["pathToSimulation"] = gridDir + 'GridOpt'

        specRatesChangeMode = literal_eval(arguments["specRatesChangeMode"])

        arguments["specRatesChangeMode"] = "('randomwalk', ['chrFusion','chrFission'])"

        # FIXME: "read" parameter file ?
        # change parameter file
        baseParaFileName = arguments['startParameterFilePath']
        paraFile = myOptimizer.readParameterFile(baseParaFileName)

        parameterNames = []
        parameterPoints = []
        for parameters in specRatesChangeMode[2]:
            parameterNames.append(parameters[0])
            parameterPoints.append(parameters[1])
        grid = itertools.product(*parameterPoints)

        for gridPoint in grid:
            print >> sys.stderr, 'Gridpoint: ' + str(gridPoint)
            for i, parameterName in enumerate(parameterNames):
                paraFile[parameterName] = gridPoint[i]

            newParaFilePath = gridDir + '_Grid.'.join(baseParaFileName.split('/')[-1].split('.'))
            myOptimizer.writeParameterFile(paraFile, newParaFilePath)
            arguments['startParameterFilePath'] = newParaFilePath

            newOpt = myOptimizer(arguments["pathToSimulation"], arguments["genesFile"],
                                 arguments["ancGenesFile"], arguments["speciesTreePath"],
                                 realScorePath=arguments["realScorePath"],
                                 scoreCalculationVerbose=arguments["scoreCalculationVerbose"]
                                 )

            newOpt.launch(arguments["magsimus1Path"],
                          chrNumbers,
                          startSpecRatesPath=arguments["startSpecRatesPath"],
                          startParameterFilePath=arguments["startParameterFilePath"],
                          numberOfSimulations=arguments["numberOfSimulations"],
                          nRepPerInput=arguments["nRepPerInput"],
                          executionMode=arguments["executionMode"],
                          specRatesChangeMode=arguments["specRatesChangeMode"],
                          deleteMagSimusOutput=arguments["deleteMagSimusOutput"]
                          )
            del newOpt

        os.remove(newParaFilePath)


if __name__ == '__main__':
    # sys.argv = ["src/analysis/myOptimizer.py",
    #             os.getcwd().replace('\\','/')+'/'+"Opt",
    #             "MagSimus/data/genesST.%s.list.bz2",
    #             "MagSimus/data/ancGenes.%s.list.bz2",
    #             "MagSimus/data/speciesTree.phylTree",
    #             "MagSimus/src/magSimus1.py", "-realScorePath=MagSimus/res/real.score",
    #             "-startSpecRatesPath=MagSimus/data/specRates_MS1.v83",
    #             "-numberOfSimulations=1",
    #             "-nRepPerInput=1",
    #             "-startParameterFilePath=MagSimus/data/parametersG.v84",
    #             "-executionMode=''",
    #             "-specRatesChangeMode=('grid', [], [('chr:invDist', ['vonMises']), ('chr:maxLengthOfInvertedSegments', [1330.0]), ('chr:absOrRelChoiceOfInvLength', ['relative']), ('chr:invDistShape', [19.0, 25.0])])",#-specRatesChangeMode=('abc', [], [('chr:invDistShape', [1, 250])], 4)",#('grid', [], [('chr:invDistMean', [-0.5, 0]), ('chr:invDistShape', [0, 1.])], 1.0)",#('grid', [('chrInvert', [-0.25, 0, 0.25])], [], 0.0001)",#'abc', [], [('chr:invDistMean', [1, 1]), ('chr:invDistShape', [1, 250])], 2)",#('randomwalk', ['chrFusion','chrFission'])",
    #             "-deleteMagSimusOutput=False",
    #             "-scoreCalculationVerbose=True"]

    # # Execution on LibsDyogen in dir Libs:
    # python MagSimus/src/analysis/myOptimizer.py Opt MagSimus/data/genesST.%s.list.bz2 MagSimus/data/ancGenes.%s.list.bz2 MagSimus/data/speciesTree.phylTree MagSimus/src/magSimus1.py -startSpecRatesPath=MagSimus/data/specRates_MS1.v80 -startParameterFilePath=MagSimus/data/parameters.v80 -scoreCalculationVerbose=False -realScorePath=MagSimus/res/real.score -executionMode=condor -nRepPerInput=1 -specRatesChangeMode="('grid', [('chrInvert', [-0.25, 0, 0.25, 0.5])], [], 1.)" -deleteMagSimusOutput=True > optimizer.stdout 2> optimizer.stderr & sleep 2 && tail -f optimizer.stderr

    # python MagSimus/src/analysis/myOptimizer.py Opt MagSimus/data/genesST.%s.list.bz2 MagSimus/data/ancGenes.%s.list.bz2 MagSimus/data/speciesTree.phylTree MagSimus/src/magSimus1.py -startSpecRatesPath=MagSimus/data/specRates_MS1.v83 -startParameterFilePath=MagSimus/data/parametersG.v83 -scoreCalculationVerbose=False -executionMode=condor -numberOfSimulations=10 -nRepPerInput=100 -specRatesChangeMode="('randomwalk', ['chrFusion','chrFission'])" -deleteMagSimusOutput=False >& optimizer.stderr

    # python MagSimus/src/analysis/myOptimizer.py ABC MagSimus/data/genesST.%s.list.bz2 MagSimus/data/ancGenes.%s.list.bz2 MagSimus/data/speciesTree.phylTree MagSimus/src/magSimus1.py -startSpecRatesPath=MagSimus/data/specRates_MS1.v80 -startParameterFilePath=MagSimus/data/parametersG.v80 -scoreCalculationVerbose=False -realScorePath=MagSimus/res/real.score -executionMode=condor -numberOfSimulations=10 -nRepPerInput=50 -specRatesChangeMode="('abc', [], [('chr:invDistMean', [0.0001, 10]), ('chr:invDistShape', [0.0001, 10])], 10)" -deleteMagSimusOutput=False >& optimizer.stderr
    # python MagSimus/src/analysis/myOptimizer.py ABC MagSimus/data/genesST.%s.list.bz2 MagSimus/data/ancGenes.%s.list.bz2 MagSimus/data/speciesTree.phylTree MagSimus/src/magSimus1.py -startSpecRatesPath=MagSimus/data/specRates_MS1.v82 -startParameterFilePath=MagSimus/data/parametersG.v80 -scoreCalculationVerbose=False -executionMode=condor -numberOfSimulations=8 -nRepPerInput=70 -specRatesChangeMode="('abc', [], [('chr:invDistMean', [0.1, 1]), ('chr:invDistShape', ['chr:invDistMean'/100, 'chr:invDistMean'])], 1000)" -randomSeed=29 -deleteMagSimusOutput=False >& optimizer.stderr
    # python MagSimus/src/analysis/myOptimizer.py ABC MagSimus/data/genesST.%s.list.bz2 MagSimus/data/ancGenes.%s.list.bz2 MagSimus/data/speciesTree.phylTree MagSimus/src/magSimus1.py -startSpecRatesPath=MagSimus/data/specRates_MS1.v83 -startParameterFilePath=MagSimus/data/parametersG.v83 -scoreCalculationVerbose=False -executionMode=condor -numberOfSimulations=7 -nRepPerInput=100 -specRatesChangeMode="('abc', [], [('chr:invDistShape', [1, 250])], 1000)" -randomSeed=7 -deleteMagSimusOutput=False >& gamma_v83_abc_7.stderr
    # (python MagSimus/src/analysis/myOptimizer.py ABC MagSimus/data/genesST.%s.list.bz2 MagSimus/data/ancGenes.%s.list.bz2 MagSimus/data/speciesTree.phylTree MagSimus/src/magSimus1.py -startSpecRatesPath=MagSimus/data/specRates_MS1.v83 -startParameterFilePath=MagSimus/data/parametersG.v83 -scoreCalculationVerbose=False -executionMode=condor -numberOfSimulations=7 -nRepPerInput=100 -specRatesChangeMode="('abc', [], [('chr:invDistShape', [1, 250])], 100)" -randomSeed=48 -deleteMagSimusOutput=False > gamma_v83_abc_48.stdout) >& gamma_v83_abc_48.stderr
    # (python MagSimus/src/analysis/myOptimizer.py Grid MagSimus/data/genesST.%s.list.bz2 MagSimus/data/ancGenes.%s.list.bz2 MagSimus/data/speciesTree.phylTree MagSimus/src/magSimus1.py -realScorePath=real.score -startSpecRatesPath=MagSimus/data/specRates_MS1.v84 -startParameterFilePath=MagSimus/data/parametersG.v84 -scoreCalculationVerbose=False -executionMode=condor -numberOfSimulations=7 -nRepPerInput=100 -specRatesChangeMode="('grid', [], [('chr:invDist', ['vonMises']), ('chr:maxLengthOfInvertedSegments', [1330.0]), ('chr:absOrRelChoiceOfInvLength', ['relative']), ('chr:invDistMean', [0.005, 0.01, 0.015]), ('chr:invDistShape', [2.5, 3.0, 3.5])])" -randomSeed=42 -deleteMagSimusOutput=False > mises.stdout) >& mises.stderr

    arguments = myTools.checkArgs([("pathToSimulation", str),
                                   ("genesFile", str),
                                   ("ancGenesFile", str),
                                   ("speciesTreePath", str),
                                   ("magsimus1Path", str)
                                   ],
                                  [("realScorePath", str, ''),
                                   ("startSpecRatesPath", str, ''),
                                   ("startParameterFilePath", str, ''),
                                   ("specRatesChangeMode", str, "('constant',)"),
                                   ("executionMode", str, ''),
                                   ("numberOfSimulations", int, 1),
                                   ("nRepPerInput", int, 2),
                                   ("deleteMagSimusOutput", bool, False),
                                   ("scoreCalculationVerbose", bool, False),
                                   ("randomSeed", int, 42)
                                   ],
                                  __doc__, showArgs=False)

    np.random.seed(arguments['randomSeed'])
    random.seed(arguments['randomSeed'])

    if 'abc' in arguments['specRatesChangeMode'].lower():
        # specRatesChangeMode should be ['abc',
        #                               [[SRCriteria1, [lowerBound, upperBound]],
        #                               [[ParameterCriteria1, [lowerBound, upperBound]],
        #                               nABCSeeds]
        myABC.startABC(arguments)
    elif 'grid' in arguments['specRatesChangeMode'].lower():
        # specRatesChangeMode should be ['grid',
        #                               [[SRCriteria1, [GridPoints]],
        #                               [[ParameterCriteria1, [GridPoints]]]
        myGrid.startGrid(arguments)
    else:
        # Case if normal optimization should be started

        # Initialize myOptimizer instance
        newOpt = myOptimizer(arguments["pathToSimulation"], arguments["genesFile"],
                             arguments["ancGenesFile"], arguments["speciesTreePath"],
                             realScorePath=arguments["realScorePath"],
                             scoreCalculationVerbose=arguments["scoreCalculationVerbose"]
                             )

        newOpt.launch(arguments["magsimus1Path"],
                      chrNumbers,
                      startSpecRatesPath=arguments["startSpecRatesPath"],
                      startParameterFilePath=arguments["startParameterFilePath"],
                      numberOfSimulations=arguments["numberOfSimulations"],
                      nRepPerInput=arguments["nRepPerInput"],
                      executionMode=arguments["executionMode"],
                      specRatesChangeMode=arguments["specRatesChangeMode"],
                      deleteMagSimusOutput=arguments["deleteMagSimusOutput"]
                      )

        # myOptimizer is now pickable, thus parallel executionMode works
        # import cPickle as pickle
        # with open('test.pickle', 'wb') as outFile:
        #     pickle.dump(newOpt, outFile, -1)
        #
        # with open('test.pickle', 'r') as inFile:
        #     newOpt2 = pickle.load(inFile)
