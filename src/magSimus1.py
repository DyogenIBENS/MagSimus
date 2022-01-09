#! /usr/bin/python
# -*- coding: utf-8 -*-

__doc__ = """
Simulate the evolution of an initial and artificial ancestral genome through
a species tree.
Events of evolution are:
    - genic events: de novo births, duplications (tandem and dispersed), losses
    - chromosomal rearrangements: fusions, fissions, translocations and inversions
Also simulates:
    - chromosomal rearrangements constraints
    - incomplete sequencing
    - 2X,6X coverage
"""

import sys
import copy
import random
import collections

# https://pypi.python.org/pypi/bidict/0.1.4
import bidict
import operator

import utils.myFile as myFile
import utils.myLightGenomes as myLightGenomes
from utils.myLightGenomes import OGene as OGene
import utils.myTools as myTools
import utils.myPhylTree as myPhylTree
import utils.myProteinTree as myProteinTree

import libs.myEvents as myEvents
import libs.myMagSimusTools as myMagSimusTools
import libs.myEvolutionProbas as myEvolutionProbas
import libs.myManagerOfEvents as myManagerOfEvents
import libs.myBreakpointsAnalyser as myBreakpointsAnalyser
import utils.myMaths as myMaths

import re
patternAlphas = re.compile('[^a-zA-Z]')

# Simulate incomplete sequencing (2X, 6X)
def simulateIncompleteSequencing(genome, node, speciesTree):
    assert isinstance(genome, myLightGenomes.LightGenome)

    if node in speciesTree.listSpecies:
        # Simulate incomplete sequencing
        if node in speciesTree.lstEsp2X:
            print >> sys.stderr, "Applying 2x coverage on %s ..." % node,
            # Breaks between 9000 and 11000 times a random chromosome
            # The longer the chromosome, the higher the chance it is broken
            for _ in xrange(random.randint(9000, 11000)):
                # Chromosome break
                (c, x) = myEvolutionProbas.randomIntergene(genome, removedExtremityLength=1)
                newC = myLightGenomes.newChromName(genome)
                genome[newC] = genome[c][x:]
                del genome[c][x:]
                assert len(genome[c]) > 0
            # In a 2X only 2/3 of the genome is sequenced
            selectedChroms = genome.keys()
            random.shuffle(selectedChroms)
            selectedChroms = set(selectedChroms[:(len(genome)*2)/3])
            genome = dict([(_c, _chrom) for (_c, _chrom) in genome.iteritems() if _c in selectedChroms])
            print >> sys.stderr, "OK (%d chromosomes)" % len(genome.keys())

        elif node in speciesTree.lstEsp6X:
            print >> sys.stderr, "Applying 6x coverage on %s ..." % node,
            for _ in xrange(random.randint(2000, 4000)):
                (c, x) = myEvolutionProbas.randomIntergene(genome, removedExtremityLength=1)
                newC = myLightGenomes.newChromName(genome)
                genome[newC] = genome[c][x:]
                del genome[c][x:]
                assert len(genome[c]) > 0
            print >> sys.stderr, "OK (%d chromosomes)" % len(genome)

        else:
            assert node in speciesTree.lstEspFull
    else:
        assert node in speciesTree.listAncestr
    return genome

def performEvents(iniGenome,
                  childName,
                  totalNbDups,
                  paramValues,
                  idx2eventType,
                  probaTandemDupSameOrientation=0.75,
                  invDist="vonMises",
                  absOrRelChoiceOfInvLength='absolute',
                  # maxL is the weighted average of chromsome lengths in the extant real genomes corresponding to the
                  # extant simulated genomes
                  # cf src/analysis/weightedAverageOfExtantChromLengths.py for the default value
                  maxL=1330,
                  invDistMean=0.01,
                  invDistShape=3,
                  lowerCap=0,
                  upperCap=sys.maxint,
                  fissionChrSamplingAlgo="simpleRandomSampling",
                  fusionChrSamplingAlgo="simpleRandomSampling",
                  inversionChrSamplingAlgo="simpleRandomSampling",
                  translocChrSamplingAlgo="simpleRandomSampling",
                  breakpointAnalyzer=True):

    lengthsOfInversions = []

    # This option discard to delete a gene that has been duplicated to ensure the good number of visible gene
    # duplications between two nodes.
    # Warning: be aware that this option will discard a tandem-duplication followed by a deletion of the initial gene
    # that yield a monogenic-inversion-like gene order if the duplicates has an inverted orientation.
    unalllowDeleteDupGenes = True

    lastGeneId = 0
    nbChrsIni = len(iniGenome.keys())
    nbGenesIni = sum([len(chrom) for chrom in iniGenome.values()])

    # New genome that evolves owing to genic and chromosomal events
    # Deep copy
    genome = copy.deepcopy(iniGenome)

    setEvts = {"chrInvert", "chrTransloc", "chrFusion", "chrFission", "geneDup", "geneLoss", "geneBirth"}
    assert setEvts == set(idx2eventType.inv.keys())
    setParams = setEvts | {'geneTandemDupProp'}
    assert setParams == set(paramValues.keys())


    setAncGenesLost = set()
    # paramValues: nb theoretical events + tandem duplications proportion
    # paramValues contains 4 int values, each value corresponds to a number of
    # events for a specific chromosomal rearrangement, and one float value in
    # [0, 1] that corresponds to the ratio of nbTandemDup/nbDup
    moe = myManagerOfEvents.ManagerOfEvents(paramValues,
                                            idx2eventType)
    bra = myBreakpointsAnalyser.BreakpointsAnalyser(iniGenome, isLazy=(not breakpointAnalyzer))
    geneTandemDupProp = paramValues['geneTandemDupProp']
    namesNewlyBornedGenes = set()
    newGeneNames = set()
    nbTandemDups = 0
    nbDispersedDups = 0
    nbEmptyChromosomesRemoved = 0
    duplicatesOf = {}
    lengthPropSampling='lengthPropSampling'

    nbEffEvents = dict((idEvt, 0) for idEvt in idx2eventType.keys())

    progressBar = myTools.ProgressBar(len(moe), step=1)
    evtDone = 0
    # if there are no event
    newGenome = genome
    for idEvt in moe:
        # Check that there is at least one chromosome and that all chromosomes contains at least one gene
        # Using the implicit booleanness of the empty list
        # bool([]) = False
        # bool([1,2]) = True
        if bool(genome.keys()) is False:
            raise ValueError('Empty genome, too many gene deletions')
        else:
            assert all(len(genome[c]) > 0 for c in genome.keys())

        # De novo gene birth (also called an origination)
        if idEvt == moe.geneBirthId:
            (infos, flag) = moe.handleEvent(idEvt, genome,
                                            upperCap=upperCap,
                                            geneBirthChrSamplingAlgo=lengthPropSampling)

            if infos is not None:
                (c, x, s) = infos
                newGeneName = childName[:2] + '_' + str(lastGeneId + 1)
                lastGeneId += 1
                assert flag == 'OK'
                print >> sys.stdout, "GeneBirth\t%s" % newGeneName
                # c -> c[:x] + [(newGeneId, s)] + c[x:]
                newGenome = myEvents.performInsertNewGene(genome, (newGeneName, c, x, s), keepOriginalGenome=False)
                namesNewlyBornedGenes.add(newGeneName)
                newGeneNames.add(newGeneName)

        # Gene duplication
        elif idEvt == moe.geneDupId:
            # The duplications of newly born families will be counted as de novo gene births
            # This is due to the way families are defined (from one ancestor to the other)
            (infos, flag) = moe.handleEvent(idEvt, genome,
                                            upperCap=upperCap,
                                            unallowedParentGeneNames=namesNewlyBornedGenes,
                                            probaTandemDupSameOrientation=probaTandemDupSameOrientation,
                                            geneTandemDupProp=geneTandemDupProp,
                                            parentGeneChrSamplingAlgo=lengthPropSampling,
                                            childDispDupChrSamplingAlgo=lengthPropSampling)
            if infos is not None:
                assert flag in ['TandemDuplication', 'DispersedDuplication']
                ((c1, x1), (c2, x2, s2)) = infos
                fatherGeneName = genome[c1][x1].n
                if flag == 'TandemDuplication':
                    nbTandemDups += 1
                elif flag == 'DispersedDuplication':
                    nbDispersedDups += 1
                flag = 'OK'
                totalNbDups[fatherGeneName] += 1
                newGeneName = fatherGeneName + myProteinTree.getDupSuffix(totalNbDups[fatherGeneName], upper=False)
                print >> sys.stdout, "GeneDup\t%s" % newGeneName
                # c -> c2[:x] + [(newGeneId, s2)] + c2[x:]
                newGenome = myEvents.performInsertNewGene(genome, (newGeneName, c2, x2, s2), keepOriginalGenome=False)
                if fatherGeneName not in duplicatesOf:
                    duplicatesOf[fatherGeneName] = []
                duplicatesOf[fatherGeneName].append(newGeneName)
                newGeneNames.add(newGeneName)

        # Gene Loss
        elif idEvt == moe.geneLossId:
            # unable to delete a gene that just was born or a newly paralog gene
            (infos, flag) = moe.handleEvent(idEvt, genome,
                                            lowerCap=lowerCap,
                                            unallowedGeneNames=newGeneNames | (set(duplicatesOf) if unalllowDeleteDupGenes else set()),
                                            geneLossChrSamplingAlgo=lengthPropSampling)
            if infos is not None:
                assert flag in ['OK', 'RemoveEmptyChromosome']
                if flag == 'RemoveEmptyChromosome':
                    # FIXME prevent deleting a specific chr because of a gene deletion
                    print >> sys.stderr, "Following a gene deletion, an empty chromosome has been removed"
                    nbEmptyChromosomesRemoved += 1
                    flag = 'OK'
                (c, idxG) = infos
                lostGeneName = genome[c][idxG].n
                print >> sys.stdout, "GeneLoss\t%s" % lostGeneName
                setAncGenesLost.add(lostGeneName)
                # c -> c[:x] + c[x+1:]
                newGenome = myEvents.performGeneLoss(genome, (c, idxG), keepOriginalGenome=False)
                # Warning!!! be very careful if the genome had one empty chromosome remove to not recreate it by calling
                # newGenome[c]. Since newGenome is a defaultdict(list), this call automatically creates an empty chromosome
                if c in newGenome:
                    # intergene of the new adjacency
                    x = idxG
                    bra.recordNewAdjacencyFromNewIntergene(newGenome, (c, x), ignoredGeneNames=newGeneNames)
                bra.updateGeneContentAfterGeneLoss(lostGeneName)

        # It is a fission
        elif idEvt == moe.chrFissionId:
            (infos, flag) = moe.handleEvent(idEvt, genome,
                                            lowerCap=lowerCap,
                                            fissionChrSamplingAlgo=fissionChrSamplingAlgo,
                                            ignoredGeneNames=newGeneNames)
            if infos is not None:
                assert flag == 'OK'
                (c, x) = infos
                # assert that the breakpoint is not at an extremity of the chromosome
                assert x != 0 and x != len(genome[c])
                print >> sys.stdout, "Fission\t%s:%d=(%s)" % (c, x, genome.getIntervStr(c, x))
                bra.recordBreakpoint((c, x), genome, ignoredGeneNames=newGeneNames)
                # c -> c=c[:x] & newC=c[x:]
                newGenome = myEvents.performFission(genome, (c, x), keepOriginalGenome=False)

        # It is a fusion
        elif idEvt == moe.chrFusionId:
            (infos, flag) = moe.handleEvent(idEvt, genome,
                                            upperCap=upperCap,
                                            fusionChrSamplingAlgo=fusionChrSamplingAlgo,
                                            ignoredGeneNames=newGeneNames)
            if infos is not None:
                assert flag == 'OK'
                ((c1, s1), (c2, s2)) = infos
                x1 = len(genome[c1]) - 1 if s1 == +1 else 0
                x2 = 0 if s2 == +1 else len(genome[c2]) - 1
                print >> sys.stdout, "Fusion\t%s:%d=(%s) with %s:%d=(%s)" % (c1, x1, genome.getIntervStr(c1, x1),
                                                                             c2, x2, genome.getIntervStr(c2, x2))
                oldNbChrom = len(genome.keys())
                # c1 & c2 -> c1 = [s1*c1 + s2*c2]
                lenOldChromC1 = len(genome[c1])
                lenOldChromC2 = len(genome[c2])
                newGenome = myEvents.performFusion(genome, ((c1, s1), (c2, s2)), keepOriginalGenome=False)
                bra.recordNewAdjacencyFromNewIntergene(newGenome, (c1, lenOldChromC1), ignoredGeneNames=newGeneNames)
                assert len(newGenome.keys()) == oldNbChrom - 1, infos
                assert c2 not in newGenome
                assert len(newGenome[c1]) == lenOldChromC1 + lenOldChromC2

        # It is an inversion of a chromosomal segment
        elif idEvt == moe.chrInvertId:
            (infos, flag) = moe.handleEvent(idEvt, genome,
                                            invDist=invDist,
                                            invDistMean=invDistMean,
                                            invDistShape=invDistShape,
                                            inversionChrSamplingAlgo=inversionChrSamplingAlgo,
                                            ignoredGeneNames=newGeneNames,
                                            absoluteOrRelativeChoiceOfLength=absOrRelChoiceOfInvLength,
                                            maxL=maxL)

            if infos is not None:
                assert flag == 'OK'
                (c, x1, x2) = infos
                print >> sys.stdout,  "Inversion\t%s:%d=(%s)-%d=(%s)" % (c, x1, genome.getIntervStr(c, x1),
                                                                         x2, genome.getIntervStr(c, x2))
                assert x2 > x1
                lenInv = x2 - x1
                lenUp = x1
                lenDown = len(genome[c]) - x2
                assert all(l >= 0 for l in (lenUp, lenInv, lenDown)), (lenUp, lenInv, lenDown)
                assert lenInv > 0, lenInv
                assert lenInv < len(genome[c]), (lenInv, len(genome[c]))
                assert not (lenUp == 0 and lenDown == 0)
                # statistics on the length of inversions
                lengthsOfInversions.append(lenInv)
                # TODO an option recordInvOfLen1
                #recordInvOfLen1 = True
                #if lenInv > 1 or (recordInvOfLen1 and lenInv == 1):
                bra.recordBreakpoint((c, x1), genome, ignoredGeneNames=newGeneNames)
                bra.recordBreakpoint((c, x2), genome, ignoredGeneNames=newGeneNames)
                # c -> c = [c[:x1] -1*c[x1:x2] + c[x2:]]
                newGenome = myEvents.performInversion(genome, (c, x1, x2), keepOriginalGenome=False)
                #if lenInv > 1 or (recordInvOfLen1 and lenInv == 1):
                bra.recordNewAdjacencyFromNewIntergene(newGenome, (c, x1), ignoredGeneNames=newGeneNames)
                bra.recordNewAdjacencyFromNewIntergene(newGenome, (c, x2), ignoredGeneNames=newGeneNames)
                #else:
                # TODO
                # update
                # bra.conservedAncestralBlocks
                # bra.setAncChromExtremities
                # bra.remainingAncestralGenome
                # bra.setGenesAtAncChromExtremities
                # bra.setRecoveredAncAdjacencies
                # bra.nbRecoveredAncAdjs

        # It is a reciprocal translocation
        elif idEvt == moe.chrTranslocId:
            (infos, flag) = moe.handleEvent(idEvt, genome,
                                            lowerCap=lowerCap,
                                            upperCap=upperCap,
                                            translocChrSamplingAlgo=translocChrSamplingAlgo,
                                            ignoredGeneNames=newGeneNames)
            if infos is not None:
                assert flag == 'OK'
                ((c1, s1, lenTrans1), (c2, s2, lenTrans2)) = infos
                lenChromArm1 = len(genome[c1]) - lenTrans1
                lenChromArm2 = len(genome[c2]) - lenTrans2
                breakpointIdxOnC1 = lenChromArm1 if s1 == 1 else lenTrans1
                breakpointIdxOnC2 = lenChromArm2 if s2 == 1 else lenTrans2
                print >> sys.stdout, \
                    "Translocation\t%s, %dth intergene: %s and %s, %dth intergene: %s" % \
                    (c1, breakpointIdxOnC1, genome.getIntervStr(c1, breakpointIdxOnC1),
                     c2, breakpointIdxOnC2, genome.getIntervStr(c2, breakpointIdxOnC2))
                bra.recordBreakpoint((c1, breakpointIdxOnC1), genome, ignoredGeneNames=newGeneNames)
                bra.recordBreakpoint((c2, breakpointIdxOnC2), genome, ignoredGeneNames=newGeneNames)
                # c1 & c2 -> c1 = (s1*c1)[:-lenTrans1] + (s2*c2)[-lenTrans2:]  &  c2 = (s2*c2)[:-lenTrans2] + (s1*c1)[-lenTrans1:]
                newGenome = myEvents.performReciprocalTranslocation(genome, ((c1, s1, lenTrans1), (c2, s2, lenTrans2)),
                                                                    keepOriginalGenome=False)
                bra.recordNewAdjacencyFromNewIntergene(newGenome, (c1, lenChromArm1), ignoredGeneNames=newGeneNames)
                bra.recordNewAdjacencyFromNewIntergene(newGenome, (c2, lenChromArm2), ignoredGeneNames=newGeneNames)
        else:
            raise ValueError("The event is not in %s." % idx2eventType.keys())


        if flag == 'OK':
            #print >> sys.stderr, "performed %s" % idx2eventType[idEvt]
            nbEffEvents[idEvt] += 1
            assert all([lowerCap <= len(chrom) <= upperCap for chrom in newGenome.values()])
            genome = newGenome
            evtDone += 1
            progressBar.printProgressIn(sys.stderr, evtDone)
        else:
            # If the event has not been performed, print the reason and how it
            # has been managed.
            print >> sys.stderr, flag

    nbEffEvents = dict([(idx2eventType[idEvtt], nbEvt) for (idEvtt, nbEvt) in nbEffEvents.iteritems()])

    # DEBUG assertions
    # assert that theoretically requested events have all been performed
    assert {n: v for n, v in paramValues.iteritems() if n != 'geneTandemDupProp'} == nbEffEvents, \
        "paramValues=%s\nnbEffEvents=%s" % (paramValues, nbEffEvents)
    # Conservation law of chromosomes
    nbFusion = nbEffEvents['chrFusion']
    nbFission = nbEffEvents['chrFission']
    nbChrs = len(newGenome.keys())
    assert nbChrs - nbChrsIni == nbFission - nbFusion - nbEmptyChromosomesRemoved, "%s - %s == %s - %s - %s" % (nbChrs, nbChrsIni, nbFission, nbFusion, nbEmptyChromosomesRemoved)
    # Conservation law of genes
    nbBirths = nbEffEvents['geneBirth']
    nbDups = nbEffEvents['geneDup']
    nbLosses = nbEffEvents['geneLoss']
    assert nbTandemDups + nbDispersedDups == nbDups, "%s + %s == %s" % (nbTandemDups, nbDispersedDups, nbDups)
    nbGenes = sum(len(chrom) for chrom in newGenome.values())
    assert nbGenes - nbGenesIni == nbBirths + nbDups - nbLosses, "%s - %s == %s + %s -%s" % (nbGenes, nbGenesIni, nbBirths, nbDups, nbLosses)

    # draw the a posteriori distribution of synteny block lengths (not including gene deletions nevertheless)
    # import matplotlib.pyplot as plt
    # plt.figure()
    # binwidth = 1
    # (n, bins, patches) = plt.hist(lengthsOfInversions, bins=range(min(lengthsOfInversions), max(lengthsOfInversions) + binwidth, binwidth))
    # plt.xlim((0, 300))
    # plt.show()
    return (newGenome, nbEffEvents, nbTandemDups, duplicatesOf, setAncGenesLost, totalNbDups, bra)

def duplicatesDescendingFrom(geneName, duplicatesOf):
    setOfDuplicates = set()
    if geneName in duplicatesOf:
        for dgn in duplicatesOf[geneName]:
            assert dgn not in setOfDuplicates
            setOfDuplicates.update({dgn} | duplicatesDescendingFrom(dgn, duplicatesOf))
    return setOfDuplicates

def buildFamilyOf(geneName, duplicatesOf, subFamilies):
    descendants = set()
    # family of the positional ortholog
    subFamily = subFamilies.getFamilyByName(geneName, default=None)
    if subFamily is not None:
        descendants |= set(subFamily.dns)
    # families of the duplicates, dgn: duplicate gene name
    nbDupsOnThisBranch = 0
    for dgn in duplicatesDescendingFrom(geneName, duplicatesOf):
        nbDupsOnThisBranch += 1
        descendants |= {dgn}
        subFamily = subFamilies.getFamilyByName(dgn, default=None)
        assert subFamily is not None  # since a duplicated gene is not allowed to be deletedfamilies
        if subFamily is not None:
            descendants |= set(subFamily.dns)
    family = myLightGenomes.Family(geneName, {geneName} | descendants)
    return family

# Recursive evolution/build of genomes
def evolve(speciesTree,
           node,
           genome,
           idx2eventType,
           arguments,
           branchSpecificParameters,
           totalNbDups):
    assert isinstance(genome, myLightGenomes.LightGenome)
    # Only does something for some extent genomes that are marked as badly
    # assembled
    genome = simulateIncompleteSequencing(genome, node, speciesTree)
    if node in speciesTree.items:  # If node is an ancestral species
        setAncGenesLostInAllChildren = set()
        # For each child
        families = myLightGenomes.Families()
        # ini families
        assert all(len(chrom) > 0 for chrom in genome.values())
        for g in [g for chrom in genome.values() for g in chrom]:
            families.addFamily(myLightGenomes.Family(g.n, set()))
        for (child, bLength) in speciesTree.items[node]:
            speciesTree.printBranchName(child, stream=sys.stderr)

            paramValues = branchSpecificParameters.paramVal[child]
            (newGenome, nbEffEvts, nbDupTandem, duplicatesOf, setAncGenesLost, totalNbDups, bra) = \
                performEvents(genome,
                              child,
                              totalNbDups,
                              paramValues,
                              idx2eventType,
                              probaTandemDupSameOrientation=arguments["probaTandemDupSameOrientation"],
                              invDist=arguments["chr:invDist"],
                              absOrRelChoiceOfInvLength=arguments["chr:absOrRelChoiceOfInvLength"],
                              # maxL is the weighted average of chromsome lengths in the extant real genomes corresponding to the
                              # extant simulated genomes
                              # cf src/analysis/weightedAverageOfExtantChromLengths.py for the default value
                              maxL=arguments["chr:maxLengthOfInvertedSegments"],
                              invDistMean=arguments["chr:invDistMean"],
                              invDistShape=arguments["chr:invDistShape"],
                              lowerCap=arguments["chr:lengthLowerCap"],
                              upperCap=arguments["chr:lengthUpperCap"],
                              fissionChrSamplingAlgo=arguments["chr:fissionSamplingAlgo"],
                              fusionChrSamplingAlgo=arguments["chr:fusionSamplingAlgo"],
                              inversionChrSamplingAlgo=arguments["chr:inversionSamplingAlgo"],
                              translocChrSamplingAlgo=arguments["chr:translocSamplingAlgo"],
                              breakpointAnalyzer=arguments["breakpointAnalyzer"]
                              )
            setAncGenesLostInAllChildren.update(setAncGenesLost)

            table = []
            for (idEvt, evt) in sorted(idx2eventType.iteritems(), key=lambda x: x[0]):
                nbEffEvt = nbEffEvts[evt]
                table.append(["#%s" % evt,
                              "= %s" % nbEffEvt,
                              "rate = %.2f" % (float(nbEffEvt)/bLength),
                              "randomFactor = %.2f" % branchSpecificParameters.randomFactors[child][evt]])
            propTandemDup = float(nbDupTandem)/nbEffEvts["geneDup"] if nbEffEvts["geneDup"] != 0 else None
            propTandemDupStr = "%s" % ("%.2f" % (propTandemDup*100) if propTandemDup is not None else 'NTR')
            table.append(["#propTandemDup",
                          "= " + propTandemDupStr + "%",
                          "",
                          ""])
            table.append(["#chromosomes",
                          "= %s" % len(newGenome),
                          "",
                          ""])
            table.append(["#genes",
                          "= %s" % sum([len(chrom) for chrom in newGenome.values()]),
                          "",
                          ""])
            myTools.printTable(table, sys.stderr)

            print >> sys.stderr, "chrom lengths = %s" % str(sorted([len(chrom) for chrom in newGenome.values()]))

            if arguments["breakpointAnalyzer"]:
                bra.printStats(stream=sys.stderr)

                if arguments['out:empiricalSbsAsGenomes'] != 'None':
                    sortedConservedAncestralBlocks = myLightGenomes.LightGenome()
                    for c, chrom in sorted(bra.conservedAncestralBlocks.iteritems(), key=lambda x: len(x[1]), reverse=True):
                        sortedConservedAncestralBlocks[c] = chrom
                    # sort sbs by decreasing lengths
                    sortedConservedAncestralBlocks.printIn(myFile.openFile(arguments['out:empiricalSbsAsGenomes'] % (str(node), str(child)), 'w'))

            # Recursive call, continue evolution
            (subFamilies, totalNbDups) = evolve(speciesTree,
                                                child,
                                                newGenome,
                                                idx2eventType,
                                                arguments,
                                                branchSpecificParameters,
                                                totalNbDups)
            for (id, family) in enumerate(families):
                families[id].dns.update(buildFamilyOf(family.fn, duplicatesOf, subFamilies).dns)

        # # Descendant genes (without duplication) in child species are added in families
        # # warning set(duplicatesOf.keys()) is not a subset of set([fam.fn for fam in families]), since a new
        # # duplicate may have been duplicated!
        # for family in families:
        #     cptDupsOnThisBranch = 0
        #     dns = set()
        #     for (subFamilies, duplicatesOf) in listOfSubFamilies:
        #         (tmpfamily, nbGDupsOnThisBranch) = buildFamilyOf(family.fn, duplicatesOf, subFamilies)
        #         if len(tmpfamily.dns) > 0:
        #         positionalOrtholog = tmpfamily.dns[0]
        #         assert family.fn == positionalOrtholog
        #         dns |= set(tmpfamily.dns[1:])
        #         cptDupsOnThisBranch += nbGDupsOnThisBranch
        # assert cptDupsOnThisBranch == paramValues['geneDup'], '%s = %s' % (cptDupsOnThisBranch, paramValues['geneDup'])

        # for family in families:
        #     if len(family.dns) == 0:
        #         # if the gene has been removed in its two children
        #         assert family.fn in setAncGenesLostInAllChildren
        #         # X aft means deleted after ...
        #         family.dns.append(family.fn + '_Xaft_' + speciesTree.speciesAcronym(node))
        #         assert len(family.dns) == 1

        # assert all([len(family.dns) > 0 for family in families])
        # print >> sys.stderr, "Warning, %s genes of %s have no extant genes" % (len([family for family in families if len(family.dns) == 0]), node)
        # for family in families:
        #     if len(family.dns) == 0:
        #         print >> sys.stderr, "Warning, family %s in %s has no extant genes" % (family.fn, node)
        print >> sys.stderr, "Writing %s families ..." % node,
        families.printIn(myFile.openFile(arguments["out:ancGenesFiles"] % node, 'w'))
        print >> sys.stderr, "OK"
    else:
         # the node is a leaf, an extant species
         # assert node in speciesTree.listSpecies
         families = myLightGenomes.Families()
         for (c, chrom) in genome.iteritems():
              assert len(chrom) > 0, "At least one non-empty chromosome at the end of simulation"
              # newChrom = []
              for (idxG, g) in enumerate(chrom):
                   # modern names
                   # assert isinstance(c, str)
                   # modernGeneName = speciesTree.speciesAcronym(node) + '_' + c + '_' + str(idxG)
                   # modernOGene = OGene(modernGeneName, g.s)
                   # newChrom.append(modernOGene)
                   # families.addFamily(myLightGenomes.Family(g.n, [modernOGene.n]))
                   #families.addFamily(myLightGenomes.Family(g.n, { speciesTree.fileName[node] + g.n}))
                   families.addFamily(myLightGenomes.Family(g.n, {g.n}))
              # genome[c] = newChrom

    species = speciesTree.fileName[node]
    # Write output genome
    print >> sys.stderr, "Writing %s genome ..." % node,
    # Return the name of the node in a format consistent with filenames, e.g.
    # genes.<s>.list.bz2
    genome.printIn(myFile.openFile(arguments["out:genomeFiles"] % species, 'w'))
    print >> sys.stderr, "OK"
    return (families, totalNbDups)

# Ensure that all chromosome lengths are well in [lowerCap, upperCap],
# if not, make necessary edition of chromLens
@myTools.deprecated
def ensureChromLensAbocveLowerCapAndBelowUpperCap(chromLens, lowerCap, upperCap):
    assert lowerCap != 1 and upperCap != sys.maxint
    idxChromsBelowLowerCap = [idxC for (idxC, chromLen) in enumerate(chromLens) if chromLen < lowerCap]
    idxChromsAboveUpperCap = [idxC for (idxC, chromLen) in enumerate(chromLens) if chromLen > upperCap]
    nbGenesInTooLongChromosomes = sum(idxChromsAboveUpperCap[idxC] for idxC in idxChromsAboveUpperCap)
    nbGenesNeeededInTooShortChromosomes = sum(idxChromsBelowLowerCap[idxC] for idxC in idxChromsBelowLowerCap)
    nbOfFreeGenes = 0
    if nbGenesInTooLongChromosomes >= nbGenesNeeededInTooShortChromosomes:
        for idxC in idxChromsAboveUpperCap:
            assert upperCap < chromLens[idxC]
            chromLens[idxC] = upperCap
            nbOfFreeGenes += chromLens[idxC] - upperCap
        for idxC in idxChromsBelowLowerCap:
            assert chromLens[idxC] < lowerCap
            chromLens[idxC] = lowerCap
            nbOfFreeGenes -= lowerCap - chromLens[idxC]
        if nbOfFreeGenes > 0:
            nbOfFreeGenesToSharePerChrom = 0
            while True and len(idxChromsBelowUpperCap) > 0:
                idxChromsBelowUpperCap = [idxC for (idxC, chromLen) in enumerate(chromLens)
                                          if chromLen + nbOfFreeGenesToSharePerChrom < upperCap]
                nbOfFreeGenesToSharePerChrom = nbOfFreeGenes / len(idxChromsBelowUpperCap)
            if len(idxChromsBelowUpperCap) == 0:
                raise ValueError("lowerCap=(%s genes) and upperCap(=%s genes) constraints for chromosome lengths" % (lowerCap, upperCap),
                                 "are not reachable with the initial number of chromosomes = %s" % len(chromLens))
            else:
                for idxC in idxChromsBelowUpperCap:
                    chromLens[idxC] += nbOfFreeGenesToSharePerChrom
                    nbOfFreeGenes -= nbOfFreeGenesToSharePerChrom
                indexChrMin, minChrLen = min(enumerate(chromLens), key=operator.itemgetter(1))
                chromLens[indexChrMin] = minChrLen + nbOfFreeGenes
    return chromLens

@myTools.deprecated
def createStartGenomeOld(arguments):
    # The initial and artificial genome
    # A float in [0,1] for each chromosome
    tmp = [random.random() for _ in xrange(arguments["iniChrNumber"])]
    # Number of genes by units of tmp
    facteur = arguments["iniGeneNumber"] / sum(tmp)
    # List, where elements' values are equal to the nb of genes in each chr
    # +1 is here to ensure that all chromosomes contain at least one gene
    chromLens = [int(x * facteur) + 1 for x in tmp]
    # Some gene may be added at the end of the last chromosome to be sure to have
    # the exact number of genes specified by the user
    indexChrMin, minChrLen = min(enumerate(chromLens), key=operator.itemgetter(1))
    chromLens[indexChrMin] = minChrLen + (arguments["iniGeneNumber"] - sum(chromLens))
    # Warning, so far the initial genome is not sorted

    lowerCap = arguments['chr:lengthLowerCap'],
    upperCap = arguments['chr:lengthUpperCap']
    if lowerCap != 1 and upperCap != sys.maxint:
        chromLens = ensureChromLensAbocveLowerCapAndBelowUpperCap(chromLens, lowerCap, upperCap)

    chromLens.sort()
    # chrlen.reverse()
    maxGeneId = 0
    genome = myLightGenomes.LightGenome()
    for (i, nb) in enumerate(chromLens):
        # The initial genome is filled [Chr_0,..., Chr_{N-1}] with
        # Chr_0 = [(0,+-1), ..., ((nbGenesChr_0-1),+-1)],
        # the Chr_1 = [(nbGenesChr_0, +-1), ..., (nbGenesChr_1, +-1)]
        genome[str(i)] = []
        for x in range(nb):
            genome[str(i)].append(OGene(maxGeneId, myEvolutionProbas.randomStrand()))
            maxGeneId += 1
    if not all([len(chr) >= 2 for chr in genome.values()]):
        print >> sys.stderr, """WARNING ! There is one random chromosome of only one gene at the beginning of the simulation. This is not a standard case and it may lead to abnormal number of fissions."""
    return(genome, maxGeneId)

def calcNewStartChromosomeLengths(chrLengthsAmniota = [0.14848911, 0.08729771, 0.07206847, 0.06486425, 0.05821389, 0.04958041, 0.04812293, 0.04775569, 0.04684266, 0.04534981, 0.04443612, 0.04057354, 0.03909136, 0.03725483, 0.03410035, 0.03160362, 0.02785170, 0.02502704, 0.02203661, 0.01555877, 0.01388112],
                                nAmniotaGenes = 20292):
    # Distribution of values by Hare-Niemeyer method
    # old (chicken) default: [2015, 1323, 1168, 1071, 921, 764, 521, 510, 491, 439, 413, 405, 365, 352, 347, 342, 341, 323, 311, 287, 254, 253, 250, 248, 222, 182, 172, 123, 46]
    quota = float(sum(chrLengthsAmniota)) / nAmniotaGenes
    amniotaTheoreticGenesPerChr = [chrLength/quota for chrLength in chrLengthsAmniota]
    amniotaGenesPerChr = [int(tGPC) for tGPC in amniotaTheoreticGenesPerChr]
    genesLeft = nAmniotaGenes - sum(amniotaGenesPerChr)
    decimalPlaces = [tGPC - int(tGPC) for tGPC in amniotaTheoreticGenesPerChr]
    orderOfDecimalPlaces = sorted(range(len(decimalPlaces)), key=decimalPlaces.__getitem__, reverse=True)
    for iChr in xrange(genesLeft):
        recChr = orderOfDecimalPlaces[iChr]
        amniotaGenesPerChr[recChr] += 1
    return(amniotaGenesPerChr)

def createStartGenome(arguments, iniAncestorName):
    chromLens = []
    for val in arguments["startChrLengths"]:
        try:
            chromLen = int(val)
            assert chromLen > 0
            chromLens.append(chromLen)
        except:
            raise ValueError('Please ensure that the startChrLengths are integers > 0')
    chromLens.sort()
    maxGeneId = 0
    genome = myLightGenomes.LightGenome()
    for (i, nb) in enumerate(chromLens):
        # The initial genome is filled [Chr_0,..., Chr_{N-1}] with
        # Chr_0 = [(0,+-1), ..., ((nbGenesChr_0-1),+-1)],
        # the Chr_1 = [(nbGenesChr_0, +-1), ..., (nbGenesChr_1, +-1)]
        genome[str(i)] = []
        for x in range(nb):
            genome[str(i)].append(OGene(iniAncestorName[:2] + '_' + str(i) + '_'+ str(x + 1), myEvolutionProbas.randomStrand()))
    lowerCap = arguments['chr:lengthLowerCap']
    upperCap = arguments['chr:lengthUpperCap']
    if not all([lowerCap <= chromLen <= upperCap for chromLen in chromLens]):
        print >> sys.stderr, "lowerCap = %s and upperCap = %s" % (lowerCap, upperCap)
        print >> sys.stderr, "distribution of chrom lengths  = %s" % chromLens
        raise ValueError('At least one chrom length (chromLen) in startChrLengths does not verify (lowerCap <= chromLen <= upperCap)')
    if not all([len(chr) >= 2 for chr in genome.values()]):
        print >> sys.stderr, """WARNING ! There is one random chromosome of only one gene at the beginning of the simulation. This is not a standard case and it may lead to abnormal number of fissions."""
    return genome

def launch(arguments):
    arguments = myMagSimusTools.readParameterFile(arguments)
    if len(arguments["seed"]) > 0:
        random.seed(eval(arguments["seed"]))
        import numpy as np
        np.random.seed(eval(arguments["seed"]))
        print >> sys.stderr, "Random number generator internal state initialised from arguments"

    myTools.printArguments(arguments, sys.stderr)
    # Check if mean for vonMises distribution was well chosen:
    if arguments["chr:invDist"].lower() == 'vonMises' and arguments["chr:invDistMean"] + 0.5 > 1:
        raise ValueError("vonMisesMean shouldn't be 0.5 < mean <= 1")
    # Load phylogenetic tree, the species tree
    speciesTree = myPhylTree.PhylogeneticTree(arguments["speciesTree.phylTree"])
    # This must match with arguments names
    idx2eventType = bidict.bidict(
        {1: "geneBirth", 2: "geneDup", 3: "geneLoss", 4: "chrInvert", 5: "chrTransloc", 6: "chrFusion", 7: "chrFission"})
    # event specific parameter value
    evtSpeParamValue = {}
    for param in idx2eventType.values() + ['geneTandemDupProp']:
        evtSpeParamValue[param] = arguments[param]

    branchSpecificParameters = \
        myMagSimusTools.BranchParameters(speciesTree,
                                         paramSpeParamValue=evtSpeParamValue,
                                         eventRandomVonMisesKappa=arguments['rate:eventRandomVonMisesKappa'],
                                         paramBranchSpeParamFile=arguments['userRatesFile'],
                                         b_randomAccel=arguments['b_randomAccel'],
                                         eventMaxRandomFactor=arguments['rate:eventMaxRandomFactor'],
                                         globalFactor=arguments['globalFactor'])
    genome = createStartGenome(arguments, speciesTree.root)

    # Start the simulation using a recursive function
    # Many parameters are implicitly used in evolve

    # dictionnary with the count of all the duplications of a specific gene, it is used to assign names to duplicates
    totalNbDups = collections.defaultdict(int)
    evolve(speciesTree, speciesTree.root, genome, idx2eventType, arguments, branchSpecificParameters, totalNbDups)

def getDefaultParameters():
    argsFormatList = [  # Species tree
                        ("speciesTree.phylTree", file),
                        ]

    optsFormatList = [  # File with listed parameters and values
                        ("parameterFile", str, ""),  # File with custom branch specific and chrom-rearrangement specific factors
                        ("userRatesFile", str, ""),  # Global factor of chrom-rearrangements
                        # If it is 2 times higher there are 2 times more chrom-rearrangements
                        ("globalFactor", float, 1.),

                        # ################################
                        # initial genome configuration  #
                        #################################
                        # Default is chicken configuration
                        ("startChrLengths", tuple, (3013, 1772, 1462, 1316, 1181, 1006, 977, 969, 951, 920, 902, 823, 793, 756, 692, 641, 565, 508, 447, 316, 282)),

                        ###################################
                        # Chromosomal rearrangement rates #
                        ###################################
                        # Factors specific to each type of chrom-rearrangement
                        # Warning: these names should correspond with the names of events in
                        # idx2evenType and also with the names of events in the userRatesFile
                        ("chrInvert", float, 1.), ("chrTransloc", float, 1.), ("chrFusion", float, 1.), ("chrFission", float, 1.),

                        ###################################
                        # Gene event rates                #
                        ###################################
                        ("geneLoss", float, 1.), ("geneBirth", float, 1.), ("geneDup", float, 1.), ('geneTandemDupProp', float, 1.), ('probaTandemDupSameOrientation', float, 0.75),

                        ###################################################
                        # Distribution of rearranged chromosomal segments #
                        ###################################################
                        # Choose sampling algorithm for chromosomal rearrangements
                        # Alternatives are "lengthPropSampling" or "simpleRandomSampling"
                        ("chr:fissionSamplingAlgo", str, "lengthPropSampling"),
                        ("chr:fusionSamplingAlgo", str, "simpleRandomSampling"),
                        # FIXME the sampling of chromosomes that undergo inversions seems proportional to length
                        # because it corresponds to the formation of chromosomal loops
                        ("chr:inversionSamplingAlgo", str, "lengthPropSampling"),
                        ("chr:translocSamplingAlgo", str, "simpleRandomSampling"),

                        # Non viable events are passed if True
                        # If False, all events specified in specRates are effective events
                        # recorded in the lineage and extant genomes have a constant nb of
                        # chromosomes in all simulations.
                        # ("passNonViableEvents", bool, False),

                        # Lower cap for chromosomes in genes
                        # If chromosome length falls below due to fission or translocation,
                        # the event is either postponed or ignored
                        # (depending on chr:forceChrNumber)
                        ("chr:lengthLowerCap", int, 1),  # Upper cap for chromosomes in genes
                        # If exceeded due to fusion or translocation, the event is either
                        # postponed or ignored (depending on chr:forceChrNumber)
                        ("chr:lengthUpperCap", int, sys.maxint),  # Parameters for the shape of the distribution of the lengths of
                        # rearranged segments (for inversions and translocations)
                        ("chr:invDist", str, "gamma"),  # can either be 'vonMises', 'uniform' or 'gamma'
                        # vonMises is defined in myMaths.randomValue.myVonMises.
                        # Both vonMises and Gamma are using two parameters.
                        # The pseudo-mean of the modified von-Mises should be within 0 < mean <= 0.5
                        # If 'gamma' is chosen, this represents the scale parameter theta
                        ("chr:invDistMean", float, 13.60),
                        # Concentration of the bell-shape of the von Mises
                        # The higher Kappa, the higher and thinner, the bell-shape.
                        # If 'gamma' is chosen, this represents the shape parameter k
                        ("chr:invDistShape", float, 1.0),
                        # If 'gamma' is chosen, this represents the scale-parameter theta
                        ('chr:absOrRelChoiceOfInvLength', str, 'absolute'),
                        # maxLengthOfInvertedSegments is the weighted average of chromsome lengths in the extant real genomes corresponding to the
                        # extant simulated genomes
                        # cf src/analysis/weightedAverageOfExtantChromLengths.py for the default value
                        ('chr:maxLengthOfInvertedSegments', int, 1330),
                        ###################
                        # Gene clustering #
                        ###################
                        # Ignored for MagSimus 1
                        # Ratio of the desired number of clusterised genes
                        # nb of clusterised genes / nb of genes candidate to the clustering
                        # Rq, 'clusteredGenesRatio' is dominant over 'forceClustering'.
                        # I.e. geneClusteredRatio=0.75 => ~75% of the genes are clusterised even
                        # if 'forceClustering' = True
                        ("geneClusteredRatio", float, 0.99),
                        # Force the clustering of genes that have a scoreSameEvent = 0 one after
                        # the other, after genes already clusterised with the gene of referencee.
                        # Clusterise also genes without history at the level of extant species.
                        # Rq : if arguments["forceClustering"] and
                        #       arguments["gene:clusteredGenesRatio"] == 1:
                        #       mkClusterByScore returns blocks of len = 1
                        #       In this case the clustering return one block with all the genes.
                        ("forceClustering", bool, True),
                        #############################
                        # Random variation of rates #
                        #############################
                        # Seed to initialise the random generator
                        ("seed", str, ""),
                        ("b_randomAccel", bool, False),
                        # Concentration of the bell-shape of the von Mises
                        ("rate:eventMaxRandomFactor", float, 1.5),
                        # The higher Kappa, the higher and thinner, the bell-shape.
                        ("rate:eventRandomVonMisesKappa", float, 2.),

                        # Boolean to write intermediary genomes during the evolution along a
                        # branch.
                        ("printIntermGenomes", bool, False),

                        ################
                        # Output Files #
                        ################
                        # Files containing the simulated genomes at each node and leaf of the
                        # species tree.
                        ("out:genomeFiles", str, "res/simu1/genes/genes.%s.list.bz2"),  # File containing gene families at each node and leaf of the
                        # species tree.
                        ("out:ancGenesFiles", str, "res/simu1/ancGenes/ancGenes.%s.list.bz2"),

                        # File to load an initial ancestral genome instead of designing it.
                        # Since the clustering process takes time, his can be used for a faster
                        # execution after the design of the first ancestor.
                        ("in:firstAncestralGenome", str, None),

                        # Boolean to set the breakpoint analyser into an active or lazy state, without the B.A.. The
                        # lazy state allows magSimus to run faster.
                        ("breakpointAnalyzer", bool, True),
                        ('out:empiricalSbsAsGenomes', str, 'res/simu1/sbs.genes.%s.%s.list.bz2')
                        ]
    return(argsFormatList, optsFormatList)

# Arguments
if __name__ == '__main__':

    (argsFormatList, optsFormatList) = getDefaultParameters()
    arguments = myTools.checkArgs(argsFormatList, optsFormatList, __doc__, showArgs=False)

    launch(arguments)
