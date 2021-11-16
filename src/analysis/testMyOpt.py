#!/usr/bin/python
# -*- coding: utf-8 -*-
# Execute this script in Libs (in the super-directory of MagSimus) !
import analysis.myOptimizer as myOptimizer
import numpy as np
import random
import os
import sys

np.random.seed(42)
random.seed(42)

scriptArgs = myOptimizer.myTools.checkArgs([("testOnCondor", bool)], [('realScorePath', str, '')], '', showArgs=False)
if scriptArgs["testOnCondor"]:
    executionModeText = "-executionMode=condor"
else:
    executionModeText = "-executionMode=''"

sys.argv = ["src/analysis/myOptimizer.py",
            os.getcwd().replace('\\','/')+'/'+"Opt",
            "MagSimus/data/genesST.%s.list.bz2",
            "MagSimus/data/ancGenes.%s.list.bz2",
            "MagSimus/data/speciesTree.phylTree",
            "MagSimus/src/magSimus1.py",
            "-startSpecRatesPath=MagSimus/data/specRates_MS1.v84",
            "-startParameterFilePath=MagSimus/data/parametersVM.v84",
            "-numberOfSimulations=1",
            "-nRepPerInput=1",
            executionModeText,
            "-specRatesChangeMode=('grid', [('chrInvert', [0, 0.25, 0.5])], [], 0.0001)",
            "-deleteMagSimusOutput=False",
            "-scoreCalculationVerbose=True",
            '-realScorePath=%s' % scriptArgs['realScorePath']]

arguments = myOptimizer.myTools.checkArgs([("pathToSimulation", str),
                               ("genesFile", str),
                               ("ancGenesFile", str),
                               ("speciesTreePath", str),
                               ("magsimus1Path", str)
                               ],
                              [("realScorePath", str, ''),
                               ("startSpecRatesPath", str, ''),
                               ("startParameterFilePath", str, ''),
                               ("specRatesChangeMode", str, "('Constant',)"),
                               ("executionMode", str, ''),
                               ("numberOfSimulations", int, 1),
                               ("nRepPerInput", int, 2),
                               ("deleteMagSimusOutput", bool, False),
                               ("scoreCalculationVerbose", bool, False)
                               ],
                              __doc__, showArgs = False)

# Initialize myOptimizer instance
newOpt = myOptimizer.myOptimizer(arguments["pathToSimulation"], arguments["genesFile"],
                     arguments["ancGenesFile"], arguments["speciesTreePath"],
                     realScorePath = arguments["realScorePath"],
                     scoreCalculationVerbose = arguments["scoreCalculationVerbose"]
                     )

newOpt.launch(arguments["magsimus1Path"],
              myOptimizer.chrNumbers,
              startSpecRatesPath = arguments["startSpecRatesPath"],
              startParameterFilePath = arguments["startParameterFilePath"],
              numberOfSimulations = arguments["numberOfSimulations"],
              nRepPerInput = arguments["nRepPerInput"],
              executionMode = arguments["executionMode"],
              specRatesChangeMode = arguments["specRatesChangeMode"],
              deleteMagSimusOutput = arguments["deleteMagSimusOutput"]
              )

mca = myOptimizer.myScore.myScoreComparisonArchive.createFromDir(newOpt.optimizationPath)
# if any([val['mean'] != 0 for val in mca.scores['geneEventRatesOnBranches-geneDup-LogRatio'].itervalues()]):
#     raise ValueError('Wrong duplication number in Score archive!')
# if np.round(mca.rankList['chrSizes_KS-LogStatistic']['Homo sapiens'][0]['Stat'], 6) != -2.021294: # -2.014758:
#     print >> sys.stderr, np.round(mca.rankList['chrSizes_KS-LogStatistic']['Homo sapiens'][0]['Stat'], 6)
#    raise ValueError('Something is wrong with chromosome sizes calculations !')

import os
for root, dirs, files in os.walk(newOpt.optimizationPath, topdown=False):
    for name in files:
        os.remove(os.path.join(root, name))
    for name in dirs:
        os.rmdir(os.path.join(root, name))
os.rmdir(newOpt.optimizationPath)

print >> sys.stderr, 'myOptimizer is probably working.'
