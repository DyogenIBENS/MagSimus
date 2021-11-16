# To extract comp and specRates/paramFile from all output into a new
# target dir with same folder structure

import os
import shutil
import errno
import myScore
# import sys
# sys.path.append('C:/cygwin64/home/Luc/Libs/MagSimus/src/analysis')

source = 'ABC/'
target = 'ABC_light/'
# source = 'Gamma_v83/'
# target = 'Gamma_v83_light/'
os.mkdir(target)

if source == 'ABC' or source == 'ABC/':
    if source == 'ABC':
        source += '/'
    sources = [recSource for recSource in os.listdir('.') if (recSource[:3]=='ABC' and recSource[-5:] != 'light' and os.path.isdir(recSource))]
    sourceList = sorted([source+'/'+recSource+'/' for source in sources for recSource in os.listdir(source) if os.path.isdir(source+'/'+recSource)])
    print(sourceList)
else:
    sourceList = [source]

for (iSource, source) in enumerate(sourceList):
    recTarget = target+source[:-1].replace('/', '-') + '/'
    try:
        os.mkdir(recTarget)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
    for recFile in os.listdir(source):
        if 'archive' in recFile:
            shutil.copy2(source+recFile, recTarget+recFile)
        elif 'real' in recFile and os.path.isdir(source+recFile):
            try:
                os.mkdir(recTarget+recFile)
            except OSError as exception:
                if exception.errno != errno.EEXIST:
                    raise
            shutil.copy2(source+recFile+'/real.score', recTarget+recFile+'/real.score')
        elif os.path.isdir(source+recFile) and (recFile == '0' or recFile == '7'):
            try:
                os.mkdir(recTarget+recFile)
            except OSError as exception:
                if exception.errno != errno.EEXIST:
                    raise
            for recFile2 in os.listdir(source+recFile):
                if 'specRates' in recFile2 or 'param' in recFile2:
                    targetPath = recTarget + recFile+ '/' + recFile2
                    if not os.path.exists(targetPath):
                        shutil.copy2(source + recFile + '/' + recFile2,
                                     targetPath)
                elif os.path.isdir(source+recFile+ '/' + recFile2):
                    try:
                        os.mkdir(recTarget + recFile+'/' + recFile2)
                    except OSError as exception:
                        if exception.errno != errno.EEXIST:
                            raise
                    for recFile3 in os.listdir(source + recFile+'/' + recFile2):
                        if 'comp' in recFile3:
                            targetPath = recTarget + recFile + '/' + recFile2 + '/' +recFile+'-'+recFile2+'.comp'
                            if not os.path.exists(targetPath):
                                shutil.copy2(source + recFile + '/' + recFile2+ '/' + recFile3,
                                             targetPath)
                        elif 'score' in recFile3:
                            targetPath = recTarget + recFile + '/' + recFile2 + '/' +recFile+'-'+recFile2+'.score'
                            if not os.path.exists(targetPath):
                                shutil.copy2(source + recFile + '/' + recFile2+ '/' + recFile3,
                                             targetPath)


# # To rename compID and specRates after fusing two directories
# source = 'InvTransl_Sim10_Rep1000_SR/'
# target = 'InvTransl_Sim10_Rep1000_SR/'
# for file in os.listdir(source):
#     if 'real' not in file and os.path.isdir(source+file):
#         # os.rename(source+file, target+str(int(file)+10))
#         print(file)
#         for recFile2 in os.listdir(source + file):
#             if 'param' in recFile2 or 'specRates' in recFile2:
#                 os.rename(source + file + '/' + recFile2, target + file + '/' + recFile2.split('-')[0] + '-' + file)
#             else:
#                 for recFile3 in os.listdir(source + file + '/' + recFile2):
#                     if 'comp' in recFile3:
#                         fileNameOld = source + file + '/' + recFile2+ '/' + recFile3
#                         fileNameNew = target + file + '/' + recFile2 + '/' +file+'-'+recFile2+'.comp'
#                         sc = myScore.myScoreComparison.load(fileNameOld)
#                         os.remove(fileNameOld)
#                         sc.ID = file + '-real'
#                         print(sc.ID)
#                         sc.save(fileNameNew)

# Get all ABC directories and add them to myScoreComparisonArchive
# import os
# import analysis.myScore as myScore
# allOpt = [recDir+'/' for recDir in os.listdir('.') if recDir[:4] == 'ABC_' and recDir != 'ABC_light']
# allSim = sorted([Opt+recSim for Opt in allOpt for recSim in os.listdir(Opt) if os.path.isdir(Opt+recSim)])
# mca = myScore.myScoreComparisonArchive.createFromDir(allSim, '/7-')
# mca.save('ABC.archive')
#
# mca.scores['error-M2006_Transl-LSE+M2006+PhylDiag']['Homo sapiens']
# myScore.myScoreComparisonArchive.plotScoreDistribution(mca.scores, 'dist-nInv-DiffEst', markComps = [])

res = {}
for dir in xrange(1,101):
    recScore = myScore.myScore.load('Gamma83_Opt/7/'+str(dir)+'/7-'+str(dir)+'.score')
    if len(res) == 0:
        res = recScore.chrSizes
    else:
        for spec in res:
            res[spec] = [res[spec][i] + val for (i, val) in enumerate(recScore.chrSizes[spec])]

for spec in res:
        res[spec] = [float(val)/sum(res[spec]) for val in res[spec]]

out = ''
for spec in res:
    line = ''
    recSum = 0
    for i, val in enumerate(sorted(res[spec], reverse=True)):
        recSum += val
        line += spec.lower().replace(' ','_') + ','+str(recSum)+','+str(float(i+1)/len(res[spec]))+'\n'
    out += line

out

sra=myScore.myScore('Simulations/data78/genesST/cleanGenesST/genesST.%s.list.bz2',
                   'Simulations/data78/ancGenes025/ancGenes.%s.list.bz2',
                   'Simulations/amniota-good.nwk', verbose=True)
print(1)

import os
source = "Grid_421/"
text = ''
for dir1 in os.listdir(source):
    if os.path.isdir(source+dir1):
        for dir2 in xrange(0,8):
            dir2 = str(dir2)
            if os.path.exists(source+dir1+'/'+dir2+'/'):
                for dir3 in xrange(1,101):
                    dir3 = str(dir3)
                    if not os.path.exists(source+dir1+'/'+dir2+'/'+dir3+'/'+dir2+'-'+dir3+'.comp'):
                        text += source+dir1+'/'+dir2+'/'+dir3+'/'+dir2+'-'+dir3+'.comp\n'

