# -*- coding: utf-8 -*-

import os
import sys
import copy
import random

import myEvolutionProbas
import utils.myLightGenomes as myLightGenomes
from utils.myLightGenomes import OGene as OGene

class MagSimusException(Exception):
    pass

class EmptyGenome(MagSimusException):
    pass

class LessThanTwoChromosomes(MagSimusException):
    # happens when (number of genes == 0)
    pass

class NoChromosomeWithAtLeastTwoGenes(MagSimusException):
    # happens when not any(len(chrom) >= 2 for chrom in genome.values())
    pass

class TooMuchTime(MagSimusException):
    pass

class NoAllowedChromosomes(MagSimusException):
    pass

class NotEnoughGenesInAllowedChromosomes(MagSimusException):
    pass

# Exceptions due to unviable events

class ChromosomeTooLong(MagSimusException):
    # happens when any(len(chrom) > upperCap for chrom in genome.values())
    # def __init__(self, message, chrName):
    #     # name of the chromosome too long
    #     Exception(self, message)
    #     self.chrName = chrName
    pass

class ChromosomeTooShort(MagSimusException):
    # happens when any(len(chrom) < lowerCap for chrom in genome.values())
    # def __init__(self, message, chrName):
    #     Exception(self, message)
    #     # name of the chromosome too short
    #     self.chrName = chrName
    pass

# Turn over the chromosomal segment if needed
def applyStrand(l, strand):
    assert isinstance(l, list)
    if strand == +1:
        newL = l
    else:
        newL = []
        for gene in reversed(l):
            # gene.s = - gene.s, can't set attribute of named tuples
            newS = - gene.s if gene.s is not None else None
            newL.append(OGene(gene.n, newS))
    return newL

# insert new gene 'geneName' on chromosme 'c' within the intergene of index 'x' (in [0, len(c)]) with strand 's'
# x = 0 corresponds to the left extremity of the chromosome 'c'
# x= len(c) corresponds to the right extremity of chromosome 'c'
# c -> c[:x] + [(geneName, s)] + c[x:]
def performInsertNewGene(genome, (geneName, c, x, s), updateGenomeDict=False, keepOriginalGenome=False):
    assert isinstance(genome, myLightGenomes.LightGenome)
    assert 0 <= x <= len(genome[c])
    if keepOriginalGenome:
        newGenome = copy.copy(genome)
        newGenome[c] = copy.copy(genome[c])
    else:
        newGenome = genome
    newGenome[c].insert(x, OGene(geneName, s))
    if updateGenomeDict:
        assert newGenome.withDict
        for (idx, g) in enumerate(newGenome[c][x:]):
            newGenome.g2p[g.n] = myLightGenomes.GeneP((c, x+idx))
    return newGenome

# c -> c[:x] + c[x+1:]
def performGeneLoss(genome, (c, x), updateGenomeDict=False, keepOriginalGenome=False):
    assert isinstance(genome, myLightGenomes.LightGenome)
    if keepOriginalGenome:
        newGenome = copy.copy(genome)
        newGenome[c] = copy.copy(genome[c])
    else:
        newGenome = genome
    geneNameLost = newGenome[c][x].n
    del newGenome[c][x]
    if len(newGenome[c]) == 0:
        del newGenome[c]
    if updateGenomeDict:
        assert newGenome.withDict
        del newGenome.g2p[geneNameLost]
        if c in newGenome:
            for (idx, g) in enumerate(newGenome[c][x:]):
                newGenome.g2p[g.n] = myLightGenomes.GeneP(c, x+idx)
    return newGenome

# fission on chromosome 'c', inner intergene of index 'x'
# x is in [1, len(c)-1]
# c -> c=c[:x] & newC=c[x:]
def performFission(genome, (c, x), updateGenomeDict=False, keepOriginalGenome=False):
    assert isinstance(genome, myLightGenomes.LightGenome)
    if x == 0 or x == len(genome[c]):
            # do not create a new chromosome, do nothing
            newGenome = genome
    else:
        if keepOriginalGenome:
            newGenome = copy.copy(genome)
            newGenome[c] = copy.copy(genome[c])
        else:
            newGenome = genome
        newC = myLightGenomes.newChromName(genome)
        newGenome[newC] = newGenome[c][x:]
        del newGenome[c][x:]
        if updateGenomeDict:
            assert newGenome.withDict
            for (idx, g) in enumerate(newGenome[newC]):
                newGenome.g2p[g.n] = myLightGenomes.GeneP(newC, idx)
    return newGenome

# fuse chromosome 'c1' with strand 's1' with chromosme 'c2' with strand 's2' to yield chromosome 'c1' = s1c1 + s2c2
# FIXME would it be more simple to just have: c1 & c2 -> c1 = [c1 + s2*c2], at least for updateGenomeDict it is simpler...
# c1 & c2 -> c1 = [s1*c1 + s2*c2]
def performFusion(genome, ((c1, s1), (c2, s2)), updateGenomeDict=False, keepOriginalGenome=False):
    assert isinstance(genome, myLightGenomes.LightGenome)
    if keepOriginalGenome:
        newGenome = copy.copy(genome)
        newGenome[c1] = copy.copy(genome[c1])
        newGenome[c2] = copy.copy(genome[c2])
    else:
        newGenome = genome
    newChromLeft = applyStrand(newGenome[c1], s1)
    newChromRight = applyStrand(newGenome[c2], s2)
    newChrom = newChromLeft + newChromRight
    newGenome[c1] = newChrom
    del newGenome[c2]
    if updateGenomeDict:
        assert newGenome.withDict
        for (idx, g) in enumerate(newGenome[c1]):
            newGenome.g2p[g.n] = myLightGenomes.GeneP(c1, idx)
    return newGenome

# invert the segment between the intergene of index 'x1' and the intergene of index 'x2' of chromosome 'c'
# x1 and x2 are both in [0, len(c)]
# the length of the inverted segment in genes is equal to abs(x1 - x2)
# c -> c = [c[:x1] -1*c[x1:x2] + c[x2:]]
def performInversion(genome, (c, x1, x2), updateGenomeDict=False, keepOriginalGenome=False):
    assert isinstance(genome, myLightGenomes.LightGenome)
    if keepOriginalGenome:
        newGenome = copy.copy(genome)
        newGenome[c] = copy.copy(genome[c])
    else:
        newGenome = genome
    newGenome[c][x1:x2] = applyStrand(newGenome[c][x1:x2], -1)
    if updateGenomeDict:
        assert newGenome.withDict
        for (idx, g) in enumerate(newGenome[c][x1:x2]):
            newGenome.g2p[g.n] = myLightGenomes.GeneP(c, x1+idx)
    return newGenome

# c1 & c2 -> c1 = (s1*c1)[:-lenTrans1] + (s2*c2)[-lenTrans2:]  &  c2 = (s2*c2)[:-lenTrans2] + (s1*c1)[-lenTrans1:]
#
# if s = +1:
#    s*c = c
#    (s*c)[:-lenTrans] = c[:-lenTrans]
#    (s*c)[-lenTrans:] = c[-lenTrans:]
# if s = -1:
#    s*c = inverted(c)
#    (s*c)[:-lenTrans] = inverted(c[lenTrans:])
#    (s*c)[-lenTrans:] = inverted(c[:lenTrans])
# FIXME:
# Idea: Because of updateGenomeWithDict I could save time to rather have:
# c1 & c2 -> c1 = [c1[:-lenTrans1] s2*c2[-lenTrans2:]] +   &   c2 = [c2[:-lenTrans2] + s1*c1[-lenTrans1:]]
# and update only the last idexes
# but this would never exchange the left extremities
def performReciprocalTranslocation(genome, ((c1, s1, lenTrans1), (c2, s2, lenTrans2)),
                                   updateGenomeDict=False, keepOriginalGenome=False):
    assert isinstance(genome, myLightGenomes.LightGenome)
    if keepOriginalGenome:
        newGenome = copy.copy(genome)
        newGenome[c1] = copy.copy(genome[c1])
        newGenome[c2] = copy.copy(genome[c2])
    else:
        newGenome = genome
    newC1 = applyStrand(newGenome[c1], s1)
    newC2 = applyStrand(newGenome[c2], s2)
    # translocated segments are from the 'lenTrans'th gene from the end, to the end.
    if lenTrans1 > 0:
        trans1 = newC1[-lenTrans1:]
        del newC1[-lenTrans1:]
    else:
        trans1 = []
    if lenTrans2 > 0:
        trans2 = newC2[-lenTrans2:]
        del newC2[-lenTrans2:]
    else:
        trans2 = []
    # this commented next two lines may require less ram but are less readable
    #genome[c2][:0] = applyStrand(r1, myEvolutionProbas.randomStrand())
    #genome[c1][:0] = applyStrand(r2, myEvolutionProbas.randomStrand())
    newGenome[c1] = newC1 + trans2
    newGenome[c2] = newC2 + trans1
    if updateGenomeDict:
        assert newGenome.withDict
        for (idx, g) in enumerate(newGenome[c2]):
            newGenome.g2p[g.n] = myLightGenomes.GeneP(c2, idx)
        for (idx, g) in enumerate(newGenome[c1]):
            newGenome.g2p[g.n] = myLightGenomes.GeneP(c1, idx)
    return newGenome

# De novo gene birth (also called an origination)
#FIXME check for the viable upperCap limit for chromosome lengths
def prepareGeneBirth(genome,
                     upperCap=sys.maxint,
                     geneBirthChrSamplingAlgo='lengthPropSampling'):
    assert isinstance(genome, myLightGenomes.LightGenome)
    assert geneBirthChrSamplingAlgo in ['simpleRandomSampling', 'lengthPropSampling']
    c = myEvolutionProbas.randomChromosome(genome, samplingAlgo=geneBirthChrSamplingAlgo, minNbGenes=0)
    x = random.randint(0, len(genome[c]))
    # FIXME: before insering the new gene, check for the upperCap limit, the up
    # limit of viable chromosome lengths.
    s = myEvolutionProbas.randomStrand()
    if len(genome[c]) >= upperCap:
        raise ChromosomeTooLong
    else:
        infos = (c, x, s)
        flag = 'OK'
    return (infos, flag)


# Gene duplication
#FIXME check for the viable upperCap limit for chromosome lengths
def prepareGeneDuplication(genome,
                           unallowedParentGeneNames=set([]),
                           geneTandemDupProp=0.6,
                           probaTandemDupSameOrientation=0.75,
                           upperCap=sys.maxint,
                           parentGeneChrSamplingAlgo='lengthPropSampling',
                           childDispDupChrSamplingAlgo='lengthPropSampling'):
    assert isinstance(genome, myLightGenomes.LightGenome)
    assert parentGeneChrSamplingAlgo in ['simpleRandomSampling', 'lengthPropSampling']
    assert childDispDupChrSamplingAlgo in ['simpleRandomSampling', 'lengthPropSampling']
    # FIXME is it too much ?
    iterMax = 1000
    cptIter = 0
    foundParentGene = False
    while not foundParentGene and cptIter <= iterMax:
        c1 = myEvolutionProbas.randomChromosome(genome, samplingAlgo=parentGeneChrSamplingAlgo, minNbGenes=1)
        x1 = random.randint(0, len(genome[c1])-1)
        parentGeneName = genome[c1][x1].n
        if parentGeneName not in unallowedParentGeneNames:
            foundParentGene = True
        cptIter += 1
    if not foundParentGene:
        raise TooMuchTime('Too long to find allowed Gene')
    else:
        # FIXME: before inserting the duplicate, check for the upperCap limit, the
        # up limit of viable chromosome lengths.
        # The new copy is inserted
        if random.random() <= geneTandemDupProp:
            flag = 'TandemDuplication'
            # Tandem duplicate is placed either on the left or on
            # the right of the original gene
            (c2, x2) = random.choice([(c1, x1), (c1, x1+1)])
            if random.random() <= probaTandemDupSameOrientation:
                s2 = genome[c1][x1].s
            else:
                s2 = - genome[c1][x1].s
        else:
            flag = 'DispersedDuplication'
            c2 = myEvolutionProbas.randomChromosome(genome, samplingAlgo=childDispDupChrSamplingAlgo, minNbGenes=0)
            x2 = random.randint(0, len(genome[c2]))
            s2 = myEvolutionProbas.randomStrand()

        if len(genome[c2]) >= upperCap:
            raise ChromosomeTooLong
        else:
            infos = ((c1, x1), (c2, x2, s2))
    return (infos, flag)


# Gene Loss
#FIXME check for the viable lowerCap limit for, lowt limit for the length of viable
def prepareGeneLoss(genome,
                    unallowedGeneNames=set([]),
                    lowerCap=0,
                    geneLossChrSamplingAlgo='lengthPropSampling'):
    assert isinstance(genome, myLightGenomes.LightGenome)
    # FIXME is it too much ?
    iterMax = 1000
    cptIter = 0
    foundGene = False
    while not foundGene and cptIter <= iterMax:
        c = myEvolutionProbas.randomChromosome(genome, samplingAlgo=geneLossChrSamplingAlgo, minNbGenes=1, ignoredGeneNames=unallowedGeneNames)
        x = myEvolutionProbas.random.randint(0, len(genome[c])-1)
        geneName = genome[c][x].n
        # Sometimes the lost gene must be a gene of the father and cannot be a duplicated gene
        if geneName not in unallowedGeneNames:
            foundGene = True
        cptIter += 1
    if not foundGene:
        raise TooMuchTime('Too long to find allowed Gene')
    else:
        if len(genome[c]) <= lowerCap:
            raise ChromosomeTooShort
        if len(genome[c]) <= 1:
            # in our case
            assert len(genome[c]) == 1
            if len(genome.keys()) == 1:
                # this should not happen if the number of deletions is constrained by the number of genes in the parent
                raise EmptyGenome
            else:
                flag = 'RemoveEmptyChromosome'
        else:
            flag = 'OK'
        infos = (c, x)
    return (infos, flag)


# Break a chromosome
def prepareFission(genome,
                   fissionChrSamplingAlgo='simpleRandomSampling',
                   ignoredGeneNames=set(),
                   lowerCap=0):
    assert isinstance(genome, myLightGenomes.LightGenome)
    # Choose a chromosome
    if all([len(chrom) == 1 for chrom in genome.values()]):
        raise NoChromosomeWithAtLeastTwoGenes('Only chromosomes of one gene')
    c = myEvolutionProbas.randomChromosome(genome, samplingAlgo=fissionChrSamplingAlgo,
                                           ignoredGeneNames=ignoredGeneNames, minNbGenes=2)
    # A random intergene location, among 'inner' intergenes, in [1,len(chr)-1]
    x = random.randint(1, len(genome[c])-1)
    chromosomesTooShort = (x < lowerCap or (len(genome[c]) - x) < lowerCap)
    if chromosomesTooShort:
        raise ChromosomeTooShort('Chromosomes too short')
    else:
        infos = (c, x)
        flag = 'OK'
    return (infos, flag)


# Fuse two chromosomes
def prepareFusion(genome,
                  upperCap=sys.maxint,
                  fusionChrSamplingAlgo='simpleRandomSampling',
                  ignoredGeneNames=set()):
    assert isinstance(genome, myLightGenomes.LightGenome)
    isinstance(genome, myLightGenomes.LightGenome)
    assert fusionChrSamplingAlgo in ['simpleRandomSampling', 'lengthPropSampling']
    if len(genome.keys()) == 1:
        raise LessThanTwoChromosomes('Only one chromosome')
    # Two chromosomes that will be fused, both must contain at least one considered gene
    c1 = myEvolutionProbas.randomChromosome(genome, samplingAlgo=fusionChrSamplingAlgo,
                                            ignoredGeneNames=ignoredGeneNames, minNbGenes=1)
    c2 = myEvolutionProbas.randomChromosome(genome, samplingAlgo=fusionChrSamplingAlgo,
                                            unallowedChrs={c1}, ignoredGeneNames=ignoredGeneNames, minNbGenes=1)
    # both chromosomes are put end to end, they may also be inversed
    newChromLen = len(genome[c1]) + len(genome[c2])
    _ChromosomeTooLong = (upperCap < newChromLen)
    if _ChromosomeTooLong:
        raise ChromosomeTooLong
    else:
        s1 = myEvolutionProbas.randomStrand()
        s2 = myEvolutionProbas.randomStrand()
        infos = ((c1, s1), (c2, s2))
        flag = 'OK'
    return (infos, flag)


# Invert a chromosomal segment
def prepareInversion(genome,
                     invDist="vonMises",
                     invDistMean=0.01,
                     invDistShape=3,
                     inversionChrSamplingAlgo='simpleRandomSampling',
                     ignoredGeneNames=set(),
                     absoluteOrRelativeChoiceOfLength='absolute',
                     # maxL is the weighted average of chromsome lengths in the extant real genomes corresponding to the
                     # extant simulated genomes
                     # cf src/analysis/weightedAverageOfExtantChromLengths.py for the default value
                     maxL=1330):
    assert isinstance(genome, myLightGenomes.LightGenome)
    assert inversionChrSamplingAlgo in ['simpleRandomSampling', 'lengthPropSampling']
    # Host chromosome
    if all(len(chrom) == 1 for chrom in genome.values()):
        raise NoChromosomeWithAtLeastTwoGenes('Only chromosomes of one gene')

    infos = myEvolutionProbas.randomSlice(genome, invDist, invDistMean, invDistShape,
                                          inversionChrSamplingAlgo='lengthPropSampling',
                                          removedExtremityLength=0,
                                          ignoredGeneNames=ignoredGeneNames,
                                          atLeastNbGenes=1,
                                          absoluteOrRelativeChoiceOfLength=absoluteOrRelativeChoiceOfLength,
                                          maxL=maxL)
    # infos = (ci, x1, x2)
    # with c a chromosome key
    # x1 and x2 correspond to extremal intergenes of the rearranged
    # chromosomal segment. x1 and x2 are ints in [0,len(chr)].
    # 0 and len(chr) correspond to extremities.
    #   * x1=0 corresponds to the left extremity of chr.
    #   * x1=1 corresponds to the intergene between the 1st gene and the 2nd gene.
    #   * x1=len(chr) corresponds the right extremity of chr.
    flag = 'OK'
    return (infos, flag)


# Perform a reciprocal translocation
def prepareReciprocalTranslocation(genome,
                                   translocChrSamplingAlgo="simpleRandomSampling",
                                   lowerCap=0,
                                   upperCap=sys.maxint,
                                   ignoredGeneNames=set()):
    assert isinstance(genome, myLightGenomes.LightGenome)
    isinstance(genome, myLightGenomes.LightGenome)
    assert translocChrSamplingAlgo in ['simpleRandomSampling', 'lengthPropSampling']
    if len(genome.keys()) < 2:
        raise LessThanTwoChromosomes
    else:
        # Choose two (different chromosomes) ...
        c1 = myEvolutionProbas.randomChromosome(genome, samplingAlgo=translocChrSamplingAlgo,
                                                ignoredGeneNames=ignoredGeneNames, minNbGenes=2)
        c2 = myEvolutionProbas.randomChromosome(genome, samplingAlgo=translocChrSamplingAlgo,
                                                ignoredGeneNames=ignoredGeneNames, minNbGenes=2,
                                                unallowedChrs={c1})
        # random.randint(a, b) returns a random integer N such that a <= N <= b.
        # FIXME is it too much ?
        iterMax = 1000
        cptIter = 0
        foundTransSegs = False
        while not foundTransSegs and cptIter <= iterMax:
            # len(genome[c1]) '- 1' because we do not want the whole genome to be translocated
            lenTrans1 = random.randint(0, len(genome[c1]) - 1)
            lenTrans2 = random.randint(0, len(genome[c2]) - 1)
            # bool = (there is at least one considered gene name in the translocated segment trans1)
            oneConsideredGeneInTrans1 = myEvolutionProbas.chrSegIsAllowed(genome[c1][:lenTrans1],
                                                                          ignoredGeneNames=ignoredGeneNames,
                                                                          minNbGenes=1)
            # bool = (there is at least one considered gene name in the translocated segment trans2)
            oneConsideredGeneInTrans2 = myEvolutionProbas.chrSegIsAllowed(genome[c2][:lenTrans2],
                                                                          ignoredGeneNames=ignoredGeneNames,
                                                                          minNbGenes=1)
            if oneConsideredGeneInTrans1 or oneConsideredGeneInTrans2:
                foundTransSegs = True
            if cptIter > 0:
                print >> sys.stderr, "Reciprocal Transl loop " + str(cptIter)
            cptIter += 1

        if not foundTransSegs:
            raise TooMuchTime('Too long to find allowed Gene')
        lenChromArm1 = len(genome[c1]) - lenTrans1
        lenChromArm2 = len(genome[c2]) - lenTrans2
        lenNewChrom1 = lenChromArm1 + lenTrans2
        lenNewChrom2 = lenChromArm2 + lenTrans1
        _chromosomesTooShort = (lenNewChrom1 < lowerCap or lenNewChrom2 < lowerCap)
        _chromosomesTooLong = (upperCap < lenNewChrom1 or upperCap < lenNewChrom2)
        if _chromosomesTooLong:
            raise ChromosomeTooLong
        elif _chromosomesTooShort:
            raise ChromosomeTooShort
        else:
            s1 = myEvolutionProbas.randomStrand()
            s2 = myEvolutionProbas.randomStrand()
            flag = 'OK'
            infos = ((c1, s1, lenTrans1), (c2, s2, lenTrans2))
    return (infos, flag)