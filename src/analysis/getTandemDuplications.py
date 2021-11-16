#!/usr/bin/python

# MagSimus/src/analysis
# python 2.7.8
# Author: Lucas TITTMANN, except for the function "nans". There, the source is written within the function.
# Copyright (c) 2015 IBENS/Dyogen Lucas TITTMANN and Hugues ROEST CROLLIUS
# mail : hrc(at)ens.fr
# Licence: GPL 3.0, except for function "nans"
# This is free software, you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 (GPL v3) or later and the CeCiLL v2 license in France

# This script analyses the duplication distances in modern genomes,
# and estimates the duplication distances in ancestral genomes.
# The results can be returned or written as csv-tables
# Plot functions for graphical representations are included

__doc__ = ""

import sys
import itertools
from collections import Counter
import numpy as np
import time
import copy
import warnings
import matplotlib.pyplot as plt
from itertools import cycle

import utils.myMapping as myMapping
import utils.myTools as myTools
import utils.myPhylTree as myPhylTree
import utils.myFile as myFile
import utils.myLightGenomes as myLightGenomes

def nans(shape, dtype=float):
    # Jorge Israel Pena, https://stackoverflow.com/questions/1704823/initializing-numpy-matrix-to-something-other-than-zero-or-one
    a = np.empty(shape, dtype)
    a.fill(np.nan)
    return a

def calcGeneGapNan(pos1, pos2):
    # Returns the distance between two Genes
    # pos1 = (chromosome, positionOnChromosome)
    # Returns np.finfo(np.float32).max if Genes are on different chromosome
    # Returns np.nan if either of the positions is np.nan
    if pos1 is np.nan or pos2 is np.nan:
        return(np.nan)
    if pos1[0] != pos2[0]:
        return(np.finfo(np.float32).max)
    else:
        return(abs(pos1[1] - pos2[1])-1)

def calcGenePosition(gene, genome):
    # 'gene' is the name of the gene (a string)
    # genome can be either a dict of lightGenomes or a dict of str
    # If gene not in genome, returns None

    if isinstance(genome.itervalues().next()[0], basestring):
        # Genome contains gene names as strings
        for (chr, genes) in genome.iteritems():
            pos = genes.index(gene) if gene in genes else None
            if pos is not None:
                return ((chr, pos))
    else:
        # Genome contains genes as LightGenome objects
        for (chr, genes) in genome.iteritems():
            for (pos, recGene) in enumerate(genes):
                if gene == recGene.n:
                    return((chr, pos))
    print(gene + ' not found in associated genome !')
    return(np.nan)

def calcDistMatrix(genesDict, genome, labels = ''):
    # 'genesDict' is a dict of lists {oldGeneName: modernGeneName} to be searched in the modern 'genome'
    # Returns None if less than 2 geneFamilies in genome

    # genesDict = modernGenesDict
    # genome = genomes[desc]
    # labels = ancestorGeneNames

    if len([len(geneList) for geneList in genesDict.itervalues() if len(geneList) > 0]) < 2:
        return((None,None))

    # Returns also non if for less than 2 genes the position on genome can be calculated

    if labels == '':
        outLabels = sorted(genesDict.keys())
    else:
        outLabels = labels

    genePositions = {}
    geneFamiliesWithAtLeastOneModernGene = []
    for (iOldGeneName, oldGeneName) in enumerate(outLabels):
        atLeastOneGeneFromFamilyInGenome = False
        genes = genesDict[oldGeneName]
        for gene in genes:
            pos = calcGenePosition(gene, genome)
            if pos is not np.nan:
                genePositions[gene] = pos
                atLeastOneGeneFromFamilyInGenome = True
        if atLeastOneGeneFromFamilyInGenome:
            geneFamiliesWithAtLeastOneModernGene.append(iOldGeneName)

    if len(geneFamiliesWithAtLeastOneModernGene) < 2:
        # print >> sys.stderr, 'Not enough genes to create distance matrix!'
        return ((None, None))

    distMat = nans([len(outLabels), len(outLabels)], np.float32)
    for i in xrange(0, len(geneFamiliesWithAtLeastOneModernGene)-1):
        for j in xrange(i+1, len(geneFamiliesWithAtLeastOneModernGene)):
            iGene = geneFamiliesWithAtLeastOneModernGene[i]
            jGene = geneFamiliesWithAtLeastOneModernGene[j]
            iGeneName = outLabels[iGene]
            jGeneName = outLabels[jGene]
            modernGenePairs = itertools.product(genesDict[iGeneName], genesDict[jGeneName])
            shortestDist = np.finfo(np.float32).max
            for genePair in modernGenePairs:
                if genePair[0] not in genePositions or genePair[1] not in genePositions:
                    # This is the case if the gene is in the ancFamilies indicated as being
                    # in the modern genome, yet it is not in our modern genome
                    # This can have multiple causes, e.g. filtering for non-coding genes
                    # or exclusion of genes on the Y-chromosome
                    continue

                posi = genePositions[genePair[0]]
                posj = genePositions[genePair[1]]
                if posi is np.nan or posj is np.nan: print('Problem with ' + str(genePair))
                newGap = calcGeneGapNan(posi, posj)
                shortestDist = min(shortestDist, newGap)
                if shortestDist == 0:
                    break
            distMat[iGene, jGene] = shortestDist

    return(distMat, outLabels)

def calcDupHistory(distMat, matLabels):
    n = len(matLabels)
    duplications = {}

    chains = {i: [i,i] for i in xrange(n)}

    for i in xrange(n-1):

        if np.isnan(distMat).all():
            # in distMatrix is no information left (all entries are np.nan)
            break

        aminAbs = np.nanargmin(distMat)
        aminRow = aminAbs / distMat.shape[0]
        aminCol = aminAbs % distMat.shape[1]

        duplications[(matLabels[aminRow], matLabels[aminCol])] = distMat[aminRow, aminCol]

        if chains[aminRow][0] == aminRow and chains[aminCol][0] == aminCol:
            chains[aminRow][0] = aminCol
            chains[aminCol][0] = aminRow
            distMat[aminRow, aminCol] = np.nan
        elif chains[aminRow][0] != aminRow and chains[aminCol][0] == aminCol:
            chains[aminRow][1] = aminCol
            aminRowEnd = chains[aminRow][0]
            chains[aminRowEnd][0] = aminCol
            chains[aminCol][0] = aminRowEnd
            distMat[aminRow, aminCol] = np.nan
            distMat[aminRowEnd, aminCol] = np.nan
            distMat[aminCol, aminRowEnd] = np.nan
            distMat[aminRow, :] = np.nan
            distMat[:, aminRow] = np.nan
        elif chains[aminRow][0] == aminRow and chains[aminCol][0] != aminCol:
            chains[aminCol][1] = aminRow
            aminColEnd = chains[aminCol][0]
            chains[aminColEnd][0] = aminRow
            chains[aminRow][0] = aminColEnd
            distMat[aminCol, aminRow] = np.nan
            distMat[aminColEnd, aminRow] = np.nan
            distMat[aminRow, aminColEnd] = np.nan
            distMat[aminCol, :] = np.nan
            distMat[:, aminCol] = np.nan
        elif chains[aminRow][0] != aminRow and chains[aminCol][0] != aminCol:
            aminRowEnd = chains[aminRow][0]
            aminColEnd = chains[aminCol][0]
            chains[aminColEnd][0] = aminRowEnd
            chains[aminRowEnd][0] = aminColEnd
            distMat[aminCol, aminRow] = np.nan
            distMat[aminColEnd, aminRowEnd] = np.nan
            distMat[aminRowEnd, aminColEnd] = np.nan
            distMat[aminRow, :] = np.nan
            distMat[:, aminRow] = np.nan
            distMat[aminCol, :] = np.nan
            distMat[:, aminCol] = np.nan
        else:
            raise ValueError("You shouldn't land here!")

    return(duplications)

def getModernGenes(prefix, family):
    # Returns all modern genes in gene family of specified species
    out = [modernGene for modernGene in family.dns if prefix == modernGene[:len(prefix)] and modernGene[len(prefix)+1].isdigit()]
    return(out)

def filterDeNovoGenes(genomeDesc, geneFamilies, speciesStr, parentStr):
    # Filter out all deNovo genes (e.g. if in Theria, filter out all deNovo genes between Amniota and HS)
    # and, if speciesStr is an internal node, all strict tandem duplications which arrived after the species
    # genomeDesc is the genome of the modern species (e.g. Homo sapiens) which is used as reference

    # Indicate all deNovo genes as None
    recGenomeOldNames = myMapping.labelWithFamNames(genomeDesc, geneFamilies[parentStr])
    recGenomeNoDeNovo = {}
    # Delete all Nones
    for (chr, oldNames) in recGenomeOldNames.iteritems():
        recGenomeNoDeNovo[chr] = [genomeDesc[chr][id] for (id, oldName) in enumerate(oldNames) if oldName.n is not None]

    out = recGenomeNoDeNovo

    # Delete direct tandemDups which arrived after the species, in order to minimize over estimation of gap
    if speciesStr in geneFamilies:
        # Case if species is an internal node (ancestor)

        recGenomeNoDeNovoLG = myLightGenomes.LightGenome(recGenomeNoDeNovo)
        recGenomeNoDeNovoOldNames = myMapping.labelWithFamNames(recGenomeNoDeNovoLG, geneFamilies[speciesStr])
        recGenomeNoDeNovoNoTD = {}
        for (chr, oldNames) in recGenomeNoDeNovoOldNames.iteritems():
            recGenes = []
            lastOldName = ''
            for (iOldName, oldName) in enumerate(oldNames):
                if oldName.n != lastOldName:
                    recGenes.append(oldName)
                lastOldName = oldName.n
            recGenomeNoDeNovoNoTD[chr] = recGenes

        out = recGenomeNoDeNovoNoTD

    return(myLightGenomes.LightGenome(out))

def calcModernSpeciesGeneDistances(child, genomes, speciesTree):
    recGenome = genomes[child[0]]
    recGenomeFiltered = filterDeNovoGenes(recGenome, geneFamilies, child[0], speciesTree.parent[child[0]].name)
    recGeneDupDistances = {}
    for (geneFamily, modernGenes) in child[1].iteritems():
        modernGenesDict = {modernGene:[modernGene] for modernGene in modernGenes}
        (distMat, matLabels) = calcDistMatrix(modernGenesDict, recGenomeFiltered)
        recGeneDupDistances[geneFamily] = calcDupHistory(distMat, matLabels)

    return(recGeneDupDistances)

def calcAncestorGeneDistances(child, genomes, geneFamilies, speciesTree):
    # Transform every ancestral gene which were duplicated into its family
    geneDupFamilies = {}
    for (geneFamily, modernGenes) in child[1].iteritems():
        if len(modernGenes) > 1:
            # Consider only families which were duplicated

            # For every (ancestor)GeneFamily, add the geneFamily where the name of oldGene was
            geneDupFamilies[geneFamily] = []
            for modernGene in modernGenes:
                for family in geneFamilies[child[0]]:
                    if modernGene == family.fn:
                        geneDupFamilies[geneFamily].append(family)
                        break

    # Calculate gene distance matrices for every modern descendant
    modernDescendants = list(speciesTree.getTargetsSpec(child[0]))
    recGeneDupDistances = {}
    possibleDuplicationCombinations = 0
    progress = 0
    for (oldGeneFamily, modernGeneFamilies) in geneDupFamilies.iteritems():
        # Print progress
        progress = progress + 1
        progressPourcent = round(100.0 * progress / len(geneDupFamilies), 2)
        printNow = progress % 50 == 0
        if printNow: print('GeneDupFamilies Progress: ' + str(progressPourcent) + '%')

        possibleDuplicationCombinations += len(modernGeneFamilies) - 1

        # distArchive is a dict with {species: distanceMatrixForSpeciesAsNumpyArray}
        distArchive = {}
        # ancestorGeneNames are the column labels (and row labels, resp.) of the matrices in distArchive
        ancestorGeneNames = [recFam.fn for recFam in modernGeneFamilies]

        for desc in modernDescendants:
            # desc stands for descendend, i.e. for Theria this can be Monodelphis, Human, Mouse, ...

            # Filter out all deNovo genes and strict tandemDuplications after recent species
            recGenomeFiltered = filterDeNovoGenes(genomes[desc], geneFamilies, child[0], speciesTree.parent[child[0]].name)
            # Create dictionary were the geneFamilies are replaced by all modernGeneNames of desc (e.g. ENSG0... for Human)
            prefixDesc = genomes[desc].values()[0][0].n.split('0')[0]
            modernGenesDict = {}
            for modernGeneFamily in modernGeneFamilies:
                modernGenesDict[modernGeneFamily.fn] = getModernGenes(prefixDesc, modernGeneFamily)

            (recDistMat, recLabel) = calcDistMatrix(modernGenesDict, genomes[desc], ancestorGeneNames)
            if recDistMat is not None:
                # recDistMat is non if there are less than two geneFamilies represented in the modern genome
                # this can be due to geneDeletion in the line
                distArchive[desc] = recDistMat

        if len(distArchive) > 0:
            # Case if at least for one species a distance matrix has been calculated
            with warnings.catch_warnings():
                # We have to catch warnings here because np.nanmin throws an annoying
                # warnings if one axis is only NaNs; yet we perfectly fine deal with this case
                # The only thing which must be looked at is if all distMat are full of Nans
                # This case is also handled in calcDupHistory
                warnings.simplefilter("ignore")
                minDistMatOverSpecies = np.nanmin(distArchive.values(), axis=0)
            recGeneDupDistances[oldGeneFamily] = calcDupHistory(minDistMatOverSpecies, ancestorGeneNames)
        else:
            recGeneDupDistances[oldGeneFamily] = {}

    print(str(len([combo for combos in recGeneDupDistances.values() for combo in combos])) + ' of ' + str(possibleDuplicationCombinations) + ' possible duplication distances estimated.')

    return(recGeneDupDistances)

def calcDupOnOtherChrRates(geneDupGapsCount):
    dupRates = {}
    for (recBranch, dupCounter) in geneDupGapsCount.iteritems():
        totalSum = sum(dupCounter.values())
        if np.finfo(np.float32).max in dupCounter:
            distantDups = dupCounter[np.finfo(np.float32).max]
        else:
            distantDups = 0
        dupRates[recBranch] = {'DistantDups': distantDups, 'DistantDupRatePerDup': float(distantDups)/totalSum}
    return(dupRates)

def calcTandemDupsFromDistDict(geneDupGapsCount, gapMax = 0, cumulated = True):
    tandemDupRates = {}
    for (recBranch, dupCounter) in geneDupGapsCount.iteritems():
        totalSum = sum(dupCounter.values())
        tandemDups = 0
        for (gap, dups) in dupCounter.iteritems():
            if cumulated:
                withinGapMax = gap <= gapMax
            else:
                withinGapMax = gap == gapMax
            if gap != np.finfo(np.float32).max and withinGapMax:
                tandemDups += dups
        tandemDupRates[recBranch] = {'TandemDups': tandemDups, 'TandemDupRatePerDup': float(tandemDups)/totalSum}
    return(tandemDupRates)

def tandemGap_LinePlot(tandemGapDict, CI = None, saveImageToPath='', normalized = False, plotCDF = True, xmax = None, ymax = None,
                       xmin = 0, ymin = 0, yLogScale = False,
                       plotTitle = 'Distance of duplication by branch (o: duplications on other chromosome)',
                       yTitle = 'Number of duplications', xTitle = 'Gap between duplications', showPlot = True):
    # Use \ci in KeyName in tandemGapDict to plot CI
    plt.ioff()

    tandemGapDict = copy.deepcopy(tandemGapDict)
    #lines = ["-", "--", "-.", ":"]
    lines = ["-"]
    colors = ["r", "g", "b", 'c', 'm', 'y', 'k']
    linecycler = cycle(lines)
    colorcycler = cycle(colors)
    labels = []
    legendEntries = []
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

    branchColor = {}
    for recBranch in tandemGapDict:
        if normalized:
            branchSum = sum(tandemGapDict[recBranch].values())
            tandemGapDict[recBranch] = {key: float(value)/branchSum for (key, value) in tandemGapDict[recBranch].iteritems()}
        recColor = next(colorcycler)
        recLineType = next(linecycler)
        if 'Gene' in recBranch:
            recColor = 'blue'
            recLineType = "--"
        elif 'Stat' in recBranch:
            recColor = 'red'
            recLineType = "-."
        elif 'Block' in recBranch:
            recColor = 'green'
            recLineType = ".."
        if recBranch == 'Reality':
            recColor = 'black'
            recLineType = "-"
        elif '\ci' in recBranch:
            if recLineType == '-':
                recLineType = next(linecycler)
            ax.fill_between(CI[recBranch][0], CI[recBranch][1], CI[recBranch][2], facecolor=recColor, interpolate=True, linewidth=0.0, alpha=0.1)

        branchColor[recBranch] = recColor
        recTupleList = sorted([(key, value) for (key, value) in tandemGapDict[recBranch].iteritems() if key != np.finfo(np.float32).max])
        if plotCDF:
            CDF = []
            oldCDF = 0
            for (key, pdf) in recTupleList:
                oldCDF += pdf
                CDF.append((key, oldCDF))
            recTupleList = CDF
        thisPlot, = ax.plot([key for (key, value) in recTupleList], [value for (key, value) in recTupleList],
                            recLineType, color = recColor, linewidth=4.0)
        legendEntries.append(thisPlot)
        labels.append(recBranch.replace('\ci',''))
        if ymaxAuto:
            ymax = max(ymax, max([value for (key, value) in recTupleList]))
        if xmaxAuto:
            xmax = max(xmax, max([key for (key, value) in recTupleList]))

    for (recBranch, distCounter) in tandemGapDict.iteritems():
        if np.finfo(np.float32).max in distCounter.keys():
            ax.plot(xmax, distCounter[np.finfo(np.float32).max], color=branchColor[recBranch], marker='o', linewidth=4.0)

    myPlotTitle = plotTitle
    ax.set_title(myPlotTitle, fontsize=20)
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
        fig.savefig(saveImageToPath, dpi = 150)
    if showPlot:
        plt.show()
    plt.close(fig)

def writeTandemDupCountsToCSV(tandemDupCount, saveTo, sep='\t'):
    allInterestingGaps = set()
    latestCountInBranch = {}
    out = ''
    for (branch, dupDistDict) in tandemDupCount.iteritems():
        allInterestingGaps = allInterestingGaps | set(dupDistDict.keys())
        latestCountInBranch[branch] = 0
        out += branch + sep
    out = out[:-len(sep)]
    out += '\n'
    allInterestingGaps = sorted(list(allInterestingGaps))

    for recGap in allInterestingGaps:
        if recGap == np.finfo(np.float32).max:
            gapStr = 'otherChr'
        else:
            gapStr = str(int(recGap))
        newLine = gapStr + sep
        for branch in msTDs['DuplicationGapsCount']:
            if recGap in tandemDupCount[branch]:
                latestCountInBranch[branch] += msTDs['DuplicationGapsCount'][branch][recGap]
            newLine += str(latestCountInBranch[branch]) + sep
        newLine = newLine[:-len(sep)]
        newLine += '\n'

        out += newLine

    with open(saveTo, 'w') as outFile:
        outFile.write(out)

if __name__ == '__main__':
    # sys.argv = ['src/analysis/getTandemDuplications.py', 'Simulations/data78/amniota-good.nwk', 'Simulations/data78/amniotaGood_geneEvents.txt', 'Simulations/data78/genesST/cleanGenesST/genesST.%s.list.bz2', 'Simulations/data78/ancGenes025/ancGenes.%s.list.bz2', 'Simulations/data78/amniotaGood_tandemDups']
    # sys.argv = ['src/analysis/getTandemDuplications.py', 'Simulations/data78/magsimus.nwk', 'Simulations/data78/magsimus_geneEvents.txt', 'Simulations/data78/genesST/cleanGenesST/genesST.%s.list.bz2', 'Simulations/data78/ancGenes025/ancGenes.%s.list.bz2', 'Simulations/data78/magsimus_tandemDups']

    arguments = myTools.checkArgs(
        [("speciesTree.phylTree", file),
         ("geneEvents", file),
         ("genesFile", str),
         ("ancGenesFile", str),
         ("saveTo", str)],
        [],
        __doc__
    )

    speciesTree = myPhylTree.PhylogeneticTree(arguments["speciesTree.phylTree"])

    genomes = {}
    for species in speciesTree.listSpecies:
        genomes[species] = myLightGenomes.LightGenome(arguments["genesFile"] % speciesTree.fileName[species])

    geneFamilies = {}
    for ancestor in speciesTree.listAncestr:
        geneFamilies[ancestor] = myLightGenomes.Families(arguments["ancGenesFile"] % speciesTree.fileName[ancestor])

    # Recursive function to load informations from geneEvents
    def loadRecursivlyGenomes(father):
        global alldata
        alldata[father] = []
        for (child, _) in speciesTree.items.get(father, []):
            t = f.readline().split("\t")
            # trans: {..., oldGeneName1: [newGeneName1], ..., oldGeneName2:
            # [newGeneName2: [newGeneName2a, newGeneName2b]]}
            # oldGeneName1 has been conserved in one copy from 'father' to 'child'
            # oldGeneName2 has been duplicated and has two gene duplicates in 'child'
            # print >> sys.stderr, child + str(len(t))
            trans = eval(t[0])
            # deleted: set([..., geneName, ...])
            deleted = eval(t[1])
            # gained: set([..., geneName, ...])
            gained = eval(t[2])
            alldata[father].append([child, trans, deleted, gained])
            loadRecursivlyGenomes(child)

    alldata = {}
    f = myFile.openFile(arguments["geneEvents"], "r")
    # The first row are Amniota gene names, which are not important for the calculation
    # Hence, we can just pop them and then forget them
    eval(f.readline())
    loadRecursivlyGenomes(speciesTree.root)
    f.close()

    # # Get Species name for each line
    # alldata = {}
    # lineNumber = 0
    # def createRecursivlyGenomesDict(father):
    #     global alldata
    #     global lineNumber
    #     alldata[father] = []
    #     for (child, _) in speciesTree.items.get(father, []):
    #         alldata[father].append([lineNumber, child])
    #         lineNumber += 1
    #         createRecursivlyGenomesDict(child)
    #
    # createRecursivlyGenomesDict(speciesTree.root)
    #
    # with myFile.openFile(arguments["geneEvents"], "r") as f:
    # content = f.readlines()
    # for line in content:
    #     print(len(line.split('\t')))
    # Load data normally
    # Doesn't work as in geneEvents.txt there is no name of species...
    # alldata = {}
    # for (child, parent) in speciesTree.parent.iteritems():
    # if parent.name not in alldata:
    #         alldata[parent.name] = [[child]]
    #     else:
    #         alldata[parent.name].append([child])

    geneDupGaps = {}
    for (spec, childs) in alldata.iteritems():
        if len(childs) > 0:
            # It is an internal node
            for child in childs:
                # child is one of the branches under the node
                print('##############################################')
                print('##   Now processing ' + child[0])
                print('##############################################')

                # Select only duplicated genes
                duplicatedGenesFamilies = {}
                for (geneFamily, modernGenes) in child[1].iteritems():
                    if len(modernGenes)>1:
                        duplicatedGenesFamilies[geneFamily] = modernGenes
                child[1] = duplicatedGenesFamilies

                if child[0] in genomes:
                    # Child is modern species
                    recGeneDupDistances = calcModernSpeciesGeneDistances(child, genomes, speciesTree)
                else:
                    # Child is internal node (ancestor)
                    recGeneDupDistances = calcAncestorGeneDistances(child, genomes, geneFamilies, speciesTree)
                geneDupGaps[child[0]] = recGeneDupDistances

    print('Finished!')

    geneDupGapsCount = {}
    for recBranch in geneDupGaps:
        recDists = [dist for subdic in geneDupGaps[recBranch].itervalues() for dist in subdic.itervalues()]
        geneDupGapsCount[recBranch] = dict(Counter(recDists))

    objectToSave = {'DuplicationGapsCount': geneDupGapsCount, 'DuplicationGaps': geneDupGaps}

    import cPickle as pickle
    with open(arguments['saveTo']+'.pickle', 'w') as outFile:
        pickle.dump(objectToSave, outFile, -1)

    # import cPickle as pickle
    # msTDs = pickle.load(open("Simulations/data78/magsimus_tandemDups.pickle", "rb"))
    # geneDupGapsCount = msTDs['DuplicationGapsCount']
    # geneDupGaps = msTDs['DuplicationGaps']

    writeTandemDupCountsToCSV(geneDupGapsCount, arguments['saveTo'] + '.csv')
    # writeTandemDupCountsToCSV(msTDs['DuplicationGapsCount'], arguments['saveTo'] + '.csv')

    # plotCurves(geneDupGapsCount, xmax = 10)
    # plotCurves(geneDupGapsCount, xmax = 10, plotCDF = False)
    tandemGap_LinePlot(geneDupGapsCount, xmax = 30, normalized = True, plotCDF=False)
    tandemGap_LinePlot(geneDupGapsCount, xmax = 30, normalized = True, plotCDF=True)

    nTandemDs = calcTandemDupsFromDistDict(geneDupGapsCount, gapMax = 10)
    nDistantDs = calcDupOnOtherChrRates(geneDupGapsCount)
    #
    # # Compare in Human with observed distances in modern genome
    # (distDict, genomeFilt) = myMapping.calcDupDist(genomes['Homo sapiens'], geneFamilies['Amniota'])
    # myMapping.calcTandemDupFromDist(distDict, 10, True)
    #
    # branchesInHumanLine = ['Homo sapiens', 'Euarchontoglires', 'Boreoeutheria', 'Theria']
    # sum([nTandemDs[branch]['TandemDups'] for branch in branchesInHumanLine])

    for branch in nTandemDs:
        print("'" + branch + "': " + str(round(nTandemDs[branch]['TandemDupRate'],4)) + "," )
