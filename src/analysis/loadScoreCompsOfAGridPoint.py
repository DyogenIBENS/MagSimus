#!/usr/bin/python

import analysis.myScore as myScore
import sys
import os

from utils import myMaths, myTools



#simu = '4'

# uniform (absolute) -> rho(dim2) = 3.11(1.14)
# TODO
# distribName = 'uniform_F4'
# simu = '4'
# gridOpt = 3000
# opt = ''

#vonmises (mu=0.001, kappa=2) -> 1.8(1.45)
# distribName = 'vonMises_mu0001_kappa_2'
# gridOpt = 300
# opt = ''
# gamma (shape=0.2, scale=300) -> 1.77(1.44)
# distribName = 'gamma_shape02_scale300'
# gridOpt = 1001
# opt = ''

# vonmises (mu=0.001, kappa=10) -> 1.79(1.45)
# distribName = 'vonMises_mu0001_kappa_10'
# gridOpt = 41
# opt = 7
# repet = 3
# #opt = 14
# gamma (shape=0.7, scale=10) -> 1.86(1.43)
# distribName = 'gamma_shape07_scale10'
# gridOpt = 31
# opt = 1

# distribName = 'gamma_shape01_scale300'
# gridOpt = 4000
# opt = ''

# # TODO
distribName = 'gamma_shape01_scale800_F3'
gridOpt = 9000
simu = '10'
opt = ''



# distribName = 'gamma_shape01_scale1000'
# gridOpt = 6000
# opt = ''
# distribName = 'gamma_shape005_scale2000'
# gridOpt = 7000
# opt = ''

repet = '1'

os.chdir('/home/jlucas/Workspace/Libs')






outDir = '/home/jlucas/Desktop/%s' % distribName
os.mkdir(outDir)

dirSim = 'Grid_%s/GridOpt%s/%s/' % (gridOpt, opt, simu)

#draw the WHM
bashCommand = "/home/jlucas/Libs/PhylDiag/src/wholeGenomeHomologyMatrix.py %s/genes.Homo.sapiens.list.bz2 %s/genes.Mus.musculus.list.bz2 %s/ancGenes.Euarchontoglires.list.bz2 -out:imageFileName=%s" \
              % (dirSim + repet, dirSim+repet, dirSim+repet, outDir + '/WHM.svg')
import subprocess
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output = process.communicate()[0]
print >> sys.stderr, 'WHM drawed'

# draw the distribution of cs lengths :
import analysis.myOptimizerGraphics as mOG
mOG.loadedScores = {}
combo = ('Homo.sapiens', 'Mus.musculus')
#combo = ('Mus musculus', 'Gallus gallus')
#combo = ('Gallus gallus', 'Mus musculus')
sims = {'sim' + simu: dirSim}
mOG.plotDF(sims, 'Grid_%s/GridOpt/realData/real.score' % gridOpt, combo, cumulated=False, plotOptRep=True,
           yLogScale=False, showPlot=False, xmax=50, savePath=outDir + '/', asHists=True)
print >> sys.stderr, 'distrib of cs length drawed'





sys.argv = ['', 'Grid_%s/GridOpt%s/realData/real.score' % (gridOpt, opt), 'Grid_%s/GridOpt%s/%s' % (gridOpt, opt, simu), 'Grid_%s/GridOpt%s/%s/scorArchive.score' % (gridOpt, opt, simu)]
print >> sys.stderr, sys.argv
assert len(sys.argv) == 4
pathToRealScore = sys.argv[1]
pathToThe100Simus = sys.argv[2]
pathToSavedScoreCompArchive = sys.argv[3]
withBA = True
sr = myScore.myScore.load(pathToRealScore)
scoreCompArchive = myScore.myScoreComparisonArchive(withBA=withBA)
listSimScores = []
listSimScoreComps = []
for i in range(1,101):
    pathToComp = pathToThe100Simus + '/' + str(i) + '/%s-%s.comp' % (simu, i)
    pathToScore = pathToThe100Simus + '/' + str(i) + '/%s-%s.score' % (simu, i)
    try:
        ss = myScore.myScore.load(pathToScore)
    except:
        print >> sys.stderr, 'Warning: no score for repet %s' % i
        continue
    ss.ID = str(simu) + '-' + str(i)
    listSimScores.append(ss)
    scoreComp = myScore.myScoreComparison(ss, sr)
    listSimScoreComps.append(scoreComp)
    #scoreComp = myScore.myScoreComparison.load(pathToComp)
    scoreCompArchive.addComp(scoreComp)
scoreCompArchive.calcStatistics()
scoreCompArchive.save(pathToSavedScoreCompArchive)

nbChrRatioDiff_geoMean = scoreCompArchive.scores['nChr-DiffRatio']['geoMean']['geoMean']
nbChrRatioDiff_geoSd = scoreCompArchive.scores['nChr-DiffRatio']['geoMean']['geoSd']
print >> sys.stderr, 'nb of chromosomes RatioDiff = %s (%s)' %  (nbChrRatioDiff_geoMean, nbChrRatioDiff_geoSd)

nbChrDiff_Hs_mean = scoreCompArchive.scores['nChr-Diff']['Homo.sapiens']['mean']
nbChrDiff_Hs_sd = scoreCompArchive.scores['nChr-Diff']['Homo.sapiens']['sd']
print >> sys.stderr, 'nb of chromosomes Diff Hs = %s (%s)' % (nbChrDiff_Hs_mean, nbChrDiff_Hs_sd)

distribChrLengthsRatioDiff_geoMean = scoreCompArchive.scores['ks_chrLengths-DiffRatio']['geoMean']['geoMean']
distribChrLengthsRatioDiff_geoSd = scoreCompArchive.scores['ks_chrLengths-DiffRatio']['geoMean']['geoSd']
print >> sys.stderr, 'distribution of chromosomes RatioDiff = %s (%s)' % (distribChrLengthsRatioDiff_geoMean, distribChrLengthsRatioDiff_geoSd)
distribChrLengthsRatioDiff_Hs_geoMean = scoreCompArchive.scores['ks_chrLengths-DiffRatio']['Homo.sapiens']['geoMean']
distribChrLengthsRatioDiff_Hs_geoSd = scoreCompArchive.scores['ks_chrLengths-DiffRatio']['Homo.sapiens']['geoSd']
print >> sys.stderr, 'distribution of chromosomes RatioDiff Hs = %s (%s)' % (distribChrLengthsRatioDiff_Hs_geoMean, distribChrLengthsRatioDiff_Hs_geoSd)


nbSbsRatioDiff_geoMean = scoreCompArchive.scores['dist-nSB-DiffRatio']['geoMean']['geoMean']
nbSbsRatioDiff_geoSd = scoreCompArchive.scores['dist-nSB-DiffRatio']['geoMean']['geoSd']
print >> sys.stderr, 'nb of cs RatioDiff = %s (%s)' % (nbSbsRatioDiff_geoMean, nbSbsRatioDiff_geoSd)
nbSbsDiff_Hs_Mm_mean = scoreCompArchive.scores['dist-nSB-Diff'][('Homo.sapiens', 'Mus.musculus')]['mean']
nbSbsDiff_Hs_Mm_sd = scoreCompArchive.scores['dist-nSB-Diff'][('Homo.sapiens', 'Mus.musculus')]['sd']
print >> sys.stderr, 'nb of cs Diff Hs-Mm= %s (%s)' % (nbSbsDiff_Hs_Mm_mean, nbSbsDiff_Hs_Mm_sd)

distribSbsLengthsRatioDiff_geoMean = scoreCompArchive.scores['ks_sbsLengths-DiffRatio']['geoMean']['geoMean']
distribSbsLengthsRatioDiff_geoSd = scoreCompArchive.scores['ks_sbsLengths-DiffRatio']['geoMean']['geoSd']
print >> sys.stderr, 'distribution of cs RatioDiff = %s (%s)'  % (distribSbsLengthsRatioDiff_geoMean, distribSbsLengthsRatioDiff_geoSd)
distribSbsLengthsRatioDiff_Hs_Mm_geoMean = scoreCompArchive.scores['ks_sbsLengths-DiffRatio'][('Homo.sapiens', 'Mus.musculus')]['geoMean']
distribSbsLengthsRatioDiff_Hs_Mm_geoSd = scoreCompArchive.scores['ks_sbsLengths-DiffRatio'][('Homo.sapiens', 'Mus.musculus')]['geoSd']
print >> sys.stderr, 'distribution of cs RatioDiff Hs-Mm = %s (%s)'  % (distribSbsLengthsRatioDiff_Hs_Mm_geoMean, distribSbsLengthsRatioDiff_Hs_Mm_geoSd)

l =  [float(nbChrRatioDiff_geoMean), float(distribChrLengthsRatioDiff_geoMean), float(nbSbsRatioDiff_geoMean), float(distribSbsLengthsRatioDiff_geoMean)]
#l =  [float(distribChrLengthsRatioDiff_geoMean), float(distribSbsLengthsRatioDiff_geoMean)]
l = [myMaths.ratioAbs(v) for v in l]
generalscore_geoMean = myMaths.geoMean(l)
generalscore_geoSd = myMaths.geoSd(l)
print >> sys.stderr, 'general RatioDiff = %s (%s)' % (generalscore_geoMean, generalscore_geoSd)

listOfTableElements =[('$\\rho$', generalscore_geoMean,generalscore_geoSd),
                      ('$\\rho(c)$', nbChrRatioDiff_geoMean, nbChrRatioDiff_geoSd),
                      ('$\\rho(\\gamma)$', distribChrLengthsRatioDiff_geoMean, distribChrLengthsRatioDiff_geoSd),
                      ('$\\rho(b)$', nbSbsRatioDiff_geoMean, nbSbsRatioDiff_geoSd),
                      ('$\\rho(\\beta)$', distribSbsLengthsRatioDiff_geoMean, distribSbsLengthsRatioDiff_geoSd),
                      ('$\\Delta(c^{h})$', nbChrDiff_Hs_mean, nbChrDiff_Hs_sd),
                      ('$\\rho(\\gamma^{h})$', distribChrLengthsRatioDiff_Hs_geoMean, distribChrLengthsRatioDiff_Hs_geoSd),
                      ('$\\Delta(b^{h,s})$', nbSbsDiff_Hs_Mm_mean, nbSbsDiff_Hs_Mm_sd),
                      ('$\\rho(\\beta^{h,s})$', distribSbsLengthsRatioDiff_Hs_Mm_geoMean, distribSbsLengthsRatioDiff_Hs_Mm_geoSd)]
import tabulate as Tabulate
from tabulate import tabulate as tab
if u'$' in Tabulate.LATEX_ESCAPE_RULES:
    del(Tabulate.LATEX_ESCAPE_RULES[u'$'])
if u'\\' in Tabulate.LATEX_ESCAPE_RULES:
    del(Tabulate.LATEX_ESCAPE_RULES[u'\\'])
if u'^' in Tabulate.LATEX_ESCAPE_RULES:
    del(Tabulate.LATEX_ESCAPE_RULES[u'^'])
if u'{' in Tabulate.LATEX_ESCAPE_RULES:
    del(Tabulate.LATEX_ESCAPE_RULES[u'{'])
if u'}' in Tabulate.LATEX_ESCAPE_RULES:
    del(Tabulate.LATEX_ESCAPE_RULES[u'}'])
headers = []
table = [['mean'],['s.d']]
for (h, m, sd) in listOfTableElements:
    headers.append(h)
    table[0].append('%.2f' % m)
    table[1].append('%.2f' % sd)
with open (outDir + '/table.tex', 'w') as f:
    print >> f, tab(table, headers, tablefmt="latex")
print tab(table, headers, tablefmt="latex")


