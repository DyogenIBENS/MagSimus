# -*- coding: utf-8 -*-

import collections
import sys

import utils.myFile as myFile
import myEvolutionProbas
from utils import myPhylTree

# Class that manages the user parameters specific to each branches of the species
# tree and each event.
# "chrInvert", "chrTransloc", "chrFusion", "chrFission", "geneBirth",
# "geneDup", "geneLoss".
class BranchParameters:

    def _computeValuesForChildrenOf(self, node):
        if node not in self.speciesTree.items:
            return
        for (child, bLength) in self.speciesTree.items[node]:
            for (paramName, paramSpeVal) in self.paramSpeParamValue.iteritems():
                paramBranchSpeVal = self.paramBranchSpeParamValue[child][paramName]
                if paramName != 'geneTandemDupProp':
                    if self.b_randomAccel:
                        randomFactor =\
                            myEvolutionProbas.randomAccel(self.eventMaxRandomFactor,
                                                          self.eventRandomVonMisesKappa)
                    else:
                        randomFactor = 1.0
                    rate = self.globalFactor * paramSpeVal * paramBranchSpeVal * randomFactor
                    self.rates[child][paramName] = rate
                    self.paramVal[child][paramName] = int(round(bLength * self.rates[child][paramName]))
                    self.randomFactors[child][paramName] = randomFactor
                else:
                    assert paramName == 'geneTandemDupProp'
                    assert - 0.1 < paramBranchSpeVal and paramBranchSpeVal < 1.1
                    assert - 0.1 < paramSpeVal and paramSpeVal < 1.1
                    self.paramVal[child][paramName] = paramBranchSpeVal * paramSpeVal
            self._computeValuesForChildrenOf(child)

    def __init__(self,
                 speciesTree,
                 # event specific parameter value
                 paramSpeParamValue={"geneBirth":1, "geneDup":1, "geneLoss":1, "chrInvert":1, "chrTransloc":1, "chrFusion":1, "chrFission":1},
                 globalFactor=1,
                 paramBranchSpeParamFile=None,
                 b_randomAccel=False,
                 eventRandomVonMisesKappa=2.0,
                 eventMaxRandomFactor=1):
        # Load user defined parameters of genic and chromosomal events as well
        # as the tandemDupProp
        # paramSpeParamValue = {'chrInvert': chrInvertParameter, ...}
        self.paramSpeParamValue = paramSpeParamValue
        assert isinstance(speciesTree, myPhylTree.PhylogeneticTree)
        self.speciesTree = speciesTree
        self.b_randomAccel = b_randomAccel
        self.eventMaxRandomFactor = eventMaxRandomFactor
        self.eventRandomVonMisesKappa = eventRandomVonMisesKappa
        self.globalFactor = globalFactor

        # Value of the parameter specific to each branch
        # Set a default value of 1.0 for all the parameters specific to each branch
        float1 = lambda: 1.0
        self.paramBranchSpeParamValue = collections.defaultdict(lambda: collections.defaultdict(float1))
        if bool(paramBranchSpeParamFile):
            print >> sys.stderr, "Loading branch specific parameters from %s" % paramBranchSpeParamFile
            paramBranchSpeParamReader = myFile.myTSV.readTabular(paramBranchSpeParamFile,
                                                                 [str, str, str, float])
            for line in paramBranchSpeParamReader:
                paramNames = set(line[0].split(","))
                parent = line[1]
                children = set(line[2].split(","))
                assert all(self.speciesTree.isChildOf(child, parent) for child in children)
                value = line[3]
                for child in children:
                    assert set(speciesTree.dicLinks[parent][child]) == set(speciesTree.dicLinks[child][parent])
                    for node in set(speciesTree.dicLinks[parent][child]) - {parent}:
                            for paramName in paramNames:
                                # the ID of the branch is the name of the child
                                if child in self.paramBranchSpeParamValue and paramName in self.paramBranchSpeParamValue[child]:
                                    print >> sys.stderr,\
                                        "The file %s is ambiguous, parameter %s is defined many times for branch %s -> %s" % (paramBranchSpeParamFile, paramName, node, child)
                                self.paramBranchSpeParamValue[child][paramName] = value

        self.rates = collections.defaultdict(lambda: collections.defaultdict(float))
        self.paramVal = collections.defaultdict(lambda: collections.defaultdict(int))
        self.randomFactors = collections.defaultdict(lambda: collections.defaultdict(float))
        self._computeValuesForChildrenOf(speciesTree.root)

    def totalNbEventsBetweenSp1AndSp2(self, sp1, sp2, paramName):
        branchesChilds = set(self.speciesTree.dicLinks[sp1][sp2]) - {self.speciesTree.lastCommonAncestor([sp1, sp2])}
        total_number_of_events = 0
        for child in branchesChilds:
            total_number_of_events += self.paramVal[child][paramName]
        return total_number_of_events

def readParameterFile(arguments):
    parameterFile = arguments['parameterFile']
    if parameterFile != '':
        print >> sys.stderr, 'WARNING - Usage of parameterFile: options in command line are overwritten by values in parameterFile'
        print >> sys.stderr, "Loading parameters from %s " % parameterFile
        userRatesReader = myFile.myTSV.readTabular(parameterFile, [str, str])
        for x in userRatesReader:
            assert len(x) == 2
            # sometimes I don't know why but x[1] = 'something\r'
            # to remove the chomp, the rstrip() method is used
            x = (x[0], x[1].rstrip())
            try:
                if x[0] in {"startChrLengths",  "forceClustering", "b_randomAccel", "printIntermGenomes"}:
                    arguments[x[0]] = eval(x[1])
                elif ((x[0][:3] == 'chr' and x[0][:22] != 'chr:invDist' and x[0][:29] != 'chr:absOrRelChoiceOfInvLength')
                      or x[0][:5] == 'rate:' or x[0][:5] == 'gene:' or x[0] == 'globalFactor'):
                    arguments[x[0]] = float(eval(x[1]))
                else:
                    arguments[x[0]] = x[1]
            except ValueError:
                print >> sys.stderr, 'Wrong parameter for "' + x[0] + '" (e.g. character where should be a number). Default used.'
            except KeyError:
                print >> sys.stderr, 'Parameter "' + x[0] + '" unknown and hence ignored.'

    return(arguments)