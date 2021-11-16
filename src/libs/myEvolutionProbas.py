# -*- coding: utf-8 -*-

import random
import math
import sys

import utils.myLightGenomes as myLightGenomes
import utils.myMaths as myMaths
import utils.myTools as myTools
import myEvents as myEvents

def chrSegIsAllowed(chrom,
                    ignoredGeneNames=set(),
                    minNbGenes=1):
    nbConsideredGenes = 0
    for (gn, _) in chrom:
        if gn not in ignoredGeneNames:
            nbConsideredGenes += 1
            if minNbGenes <= nbConsideredGenes:
                return True
    return False


def chrsConsidered(genome,
                   ignoredGeneNames=set(),
                   minNbGenes=1,
                   unallowedChrs=set()):
    res = set()
    for c in set(genome.keys()) - unallowedChrs:
        if chrSegIsAllowed(genome[c], ignoredGeneNames=ignoredGeneNames, minNbGenes=minNbGenes):
            res.add(c)
    return res


# A random value between maxAccel^-1 and maxAccel^1
def randomAccel(rateMaxAccel, vmKappa):
    # WARNING!! in Python version 2.7.3, vonmisesvariate(0,K) variates from
    # -pi to +pi. This bug (as in the official documentation it goes well from
    #  0 to 2*pi) was fixed in Pyhton 2.7.4
    # FIXME : We would like the expected value to be equal to  0, but it is
    # not exactly the case
    if sys.version_info <= (2, 7, 3):
        alpha = random.vonmisesvariate(0, vmKappa) / math.pi
    else :
        alpha = random.vonmisesvariate(math.pi, vmKappa) / math.pi - 1

    if alpha > 1:
        raise ValueError("res=%s should be lower than %s" % (alpha, 1))

    res = math.pow(rateMaxAccel, alpha)
    return(res)


# Return +1 or -1 randomly
def randomStrand():
    return random.choice([-1, 1])


# FIXME: to chose gene duplicated and gene loss randomChromosome may be inappropriate...
# DEprecated: User would rather write the two lines rather than using this extremely simple function
@myTools.deprecated
def randomGene(genome,
               samplingAlgo='lengthPropSampling'):
    assert samplingAlgo in ['simpleRandomSampling', 'lengthPropSampling']
    c = randomChromosome(genome, samplingAlgo=samplingAlgo, minNbGenes=1)
    x = random.randint(0, len(genome[c])-1)
    return (c, x)


# A random chromosome, the longer the chromosome the higher is the probability
# to be chosen.
# Underlying hypothesis: a long chromosome has more chances to undergo a fission
# than a small chromosome.
def randomChromosome(genome,
                     samplingAlgo='simpleRandomSampling',
                     ignoredGeneNames=set(),
                     minNbGenes=1,
                     unallowedChrs=set(),
                     iterMax=25):
    assert isinstance(genome, myLightGenomes.LightGenome)
    assert samplingAlgo in ['simpleRandomSampling', 'lengthPropSampling']
    allowedChrs = chrsConsidered(genome,
                                 ignoredGeneNames=ignoredGeneNames,
                                 minNbGenes=minNbGenes,
                                 unallowedChrs=unallowedChrs)
    if len(allowedChrs) == 0:
        print >> sys.stderr, 'No allowed chromosomes left.'
        print >> sys.stderr, 'samplingAlgo: ' + str(samplingAlgo)
        print >> sys.stderr, 'allowedChrs ' + str(allowedChrs)
        print >> sys.stderr, 'unallowedChrs: ' + str(unallowedChrs)
        print >> sys.stderr, 'minNbGenes: ' + str(minNbGenes)
        print >> sys.stderr, 'chrLengths: ' + str([len(chr) for chr in genome.itervalues()])
        print >> sys.stdout, "Problem: NoAllowedChromosomes"
        return(None)
        # raise myEvents.NoAllowedChromosomes
    elif all(chromLen < minNbGenes for chromLen in allowedChrs):
        print >> sys.stdout, "Problem: NotEnoughGenesInAllowedChromosomes"
        return(None)
        # raise myEvents.NotEnoughGenesInAllowedChromosomes
    cptIter = 0
    if samplingAlgo == 'simpleRandomSampling':
        while cptIter <= iterMax:
            cptIter += 1
            c = random.choice(list(allowedChrs))
            if len(genome[c]) >= minNbGenes:
                return c
            else:
                print >> sys.stderr, 'randomChromosome simpleRandomSampling Loop' + str(cptIter)
        print >> sys.stdout, "Problem: Too unprobable to find a convenient random chromosome"
        return(None)
        # raise myEvents.TooMuchTime('Too unprobable to find a convenient random chromosome')
    elif samplingAlgo == 'lengthPropSampling':
        # The iteration over genome is always in the same order since
        # myLightGenomes.LightGenome inherit from myTools.DefaultOrderedDict
        cByIndex = dict([(i, c) for (i, c) in enumerate(allowedChrs)])
        # return an random index in the list, considering the values of the elements
        # that indicate relative probabilities
        relativeProbabilitiesToBeChosen = [len(genome[c]) for c in allowedChrs]
        rCidx = myMaths.randomValue.bisectChooser(relativeProbabilitiesToBeChosen)
        cptIter = 0
        while cptIter <= iterMax:
            cptIter += 1
            i = rCidx()
            c = cByIndex[i]
            if len(genome[c]) >= minNbGenes:
                return c
            else:
                print >> sys.stderr, 'randomChromosome lengthPropSampling Loop' + str(cptIter)
        print >> sys.stdout, "Problem: Too unprobable to find a convenient random chromosome"
        return(None)
        # raise myEvents.TooMuchTime('Too unprobable to find a convenient random chromosome')
    else:
        raise ValueError("Unknown 'samplingAlgo' value. 'samplingAlgo' should either be equal to 'simpleRandomSampling' "
                         "or 'lengthPropSampling'")

def distGamma(scale, shape):
    import numpy.random as nr
    theta = scale  # mean
    kappa = shape
    x = nr.gamma(kappa, theta)
    return(x)

# Before designing the distribution from the von Mises distribution the
# gamma distribution was tried by M. Muffato.
# l = int((len(chrom)-1) * random.gammavariate(arguments["chr:invertGammaAlpha"], arguments["chr:invertGammaBeta"]))

def distVonMises(vmMean, vmKappa):
    # FIXME, median is ok but not the mean, the mean of the modified von Mises.<
    x = myMaths.randomValue.myVonMises(vmMean, vmKappa)
    return x

def randomSlice(*args, **kwargs):
    absoluteOrRelativeChoiceOfLength = kwargs.get('absoluteOrRelativeChoiceOfLength', 'absolute')
    del kwargs['absoluteOrRelativeChoiceOfLength']
    if absoluteOrRelativeChoiceOfLength == 'absolute':
        return randomSliceAbsolute(*args, **kwargs)
    elif absoluteOrRelativeChoiceOfLength == 'relative':
        # maxL only has a meaning when lengths of inverted segments are chosen with absolute lengths
        # if kwargs['maxL'] is not None:
        #     print >> sys.stderr, "Warning: maxL is not None, since absoluteOrRelativeChoiceOfLength='relative', " \
        #                          "maxL is considered as None"
        del kwargs['maxL']
        return randomSliceRelative(*args, **kwargs)
    else:
        raise ValueError('absoluteOrRelativeChoiceOfLengths is either \'absolute\' or \'relative\'')

# A random chromosomal segment in the genome, using relative distributions
def randomSliceRelative(genome, ChromSegDistribution, vmMean, vmKappa,
                        inversionChrSamplingAlgo='lengthPropSampling',
                        removedExtremityLength=0,
                        ignoredGeneNames=set(),
                        atLeastNbGenes=1):
    # FIXME is it too much ?
    iterMax = 1000
    cptIter = 0
    while cptIter <= iterMax:
        if cptIter > 0:
            print >> sys.stderr, "Inversion loop " + str(cptIter)
        # 1) pick a chromosome
        ci = randomChromosome(genome,
                              samplingAlgo=inversionChrSamplingAlgo,
                              minNbGenes=2,
                              ignoredGeneNames=ignoredGeneNames)
        c = genome[ci]
        assert chrSegIsAllowed(c, ignoredGeneNames=ignoredGeneNames, minNbGenes=atLeastNbGenes)
        # nb Intergene = the first intergene on the left of the first gene + the second intergene between the first
        # gene and the second gene, ....
        nbIntergenes = len(c) + 1
        if not (0 < atLeastNbGenes <= nbIntergenes - 2 * removedExtremityLength):
            raise myEvents.ChromosomeTooShort()

        # 2) chose a length of chromosome that is relative to the length of the chosen chromosome
        lenC = len(c)
        # l is at most equal to lenC - 1, otherwise the whole chromosome is "inverted"
        maxL = lenC - 2 * removedExtremityLength if removedExtremityLength > 0 else lenC - 1
        if ChromSegDistribution.lower() == 'uniform':
            l = random.randint(atLeastNbGenes, maxL)
        else:
            if ChromSegDistribution.lower() == 'gamma':
                # FIXME please give appropriate arguments
                # This is the case when the Mazowita2005-gamma distribution should be used
                # Usually it should get its own function to be proper, but at this phase of
                # alpha-testing, we just tweek it in here
                x = distGamma(vmMean, vmKappa)
            else:
                x = distVonMises(vmMean, vmKappa)
            # x is in [0, 1]
            l = int((maxL - atLeastNbGenes) * x) + atLeastNbGenes
            # l is in [atLeastNbGenes, maxL]
        # print >> sys.stderr, str(atLeastNbGenes)
        # print >> sys.stderr, str(l)
        # print >> sys.stderr, str(maxL)
        if not atLeastNbGenes <= l <= maxL:
            cptIter += 1
            continue
        # x1 is the position of the 1st intergene:
        #       - x=0 corresponds the intergene at the left of the 1st gene at the
        #         extremity of the chrom.
        #       - x=1 corresponds the intergene between the 1st and the 2nd gene.
        x1 = random.randint(0, lenC - l)
        # x1 is in [|0, lenC - l|]
        x2 = x1 + l
        if (x1 != x2 and x1 >= 0 and x2 <= lenC and (x1 != 0 or x2 != lenC) and
            chrSegIsAllowed(c[x1:x2], ignoredGeneNames=ignoredGeneNames, minNbGenes=atLeastNbGenes)):
            return (ci, x1, x2)
        else:
            cptIter += 1
    raise myEvents.TooMuchTime('Too long to find allowed Slice')

def randomSliceAbsolute(genome, ChromSegDistribution, vmMean, vmKappa,
                        inversionChrSamplingAlgo='lengthPropSampling',
                        removedExtremityLength=0,
                        ignoredGeneNames=set(),
                        atLeastNbGenes=1,
                        # maxL is the weighted average of chromsome lengths in the extant real genomes corresponding to the
                        # extant simulated genomes
                        # cf src/analysis/weightedAverageOfExtantChromLengths.py for the default value
                        maxL=1330):
    # FIXME is it too much ?
    iterMax1 = 100
    cptIter1 = 0
    iterMax2 = 10
    while cptIter1 <= iterMax1:
        cptIter2 = 0
        # 1) first yield an inverted segment length
        if ChromSegDistribution.lower() == 'uniform':
            # x is in [0,1)
            x = random.random()
            l = int((maxL - atLeastNbGenes) * x + atLeastNbGenes)
            # l is in [atLeastNbGenes, maxL]
        elif ChromSegDistribution.lower() == 'gamma':
            # FIXME please give appropriate arguments
            # This is the case when the Mazowita2005-gamma distribution should be used
            # Usually it should get its own function to be proper, but at this phase of
            # alpha-testing, we just tweek it in here
            l = int(distGamma(vmMean, vmKappa)) + atLeastNbGenes
            # l is in [atLeastNbGenes, Infinity]
        elif ChromSegDistribution.lower() == 'vonmises':
            # x is in [0,1)
            x = distVonMises(vmMean, vmKappa)
            l = int((maxL - atLeastNbGenes) * x + atLeastNbGenes)
            # l is in [atLeastNbGenes, maxL]
        else:
            raise ValueError('invDist should be in [\'uniform\', \'gamma\', \'vonMises\']')

        # 2) then pick a chromosome that is long enough
        while cptIter2 <= iterMax2:
            # lenC is at least equal to l + 1, otherwise the whole chromosome is inverted
            ci = randomChromosome(genome,
                                  samplingAlgo=inversionChrSamplingAlgo,
                                  minNbGenes=(l + 1),
                                  ignoredGeneNames=ignoredGeneNames)
            if ci is None:
                # Case if it was not possible to find chromosome with good size
                # Draw again from gamma
                break
            c = genome[ci]
            lenC = len(c)
            assert chrSegIsAllowed(c, ignoredGeneNames=ignoredGeneNames, minNbGenes=atLeastNbGenes)
            # nb Intergene = the first intergene on the left of the first gene + the second intergene between the first
            # gene and the second gene, ....
            nbIntergenes = lenC + 1
            if not (0 < atLeastNbGenes <= nbIntergenes - 2 * removedExtremityLength):
                raise myEvents.ChromosomeTooShort()
            x1 = random.randint(0, lenC - l)
            x2 = x1 + l
            if (x1 != x2 and x1 >= 0 and x2 <= lenC and (x1 != 0 or x2 != lenC) and
                chrSegIsAllowed(c[x1:x2], ignoredGeneNames=ignoredGeneNames, minNbGenes=atLeastNbGenes)):
                return (ci, x1, x2)

            cptIter2 += 1
        cptIter1 += 1
    raise myEvents.TooMuchTime('Too long to find allowed Slice: Chromosomes are too small, chosen slices are too long '
                               'to find chromomes that can undergo an inversion')


# Return a list of 'n' chromosome lengths.
# The sum of the chrom lengths is 's'.
# Each length is at least equal to 'm'.
@myTools.deprecated
def randomSizes(n, s, m):
    while True:
        tmp = [random.random() for _ in xrange(n)]
        facteur = s / sum(tmp)
        lengths = [max(int(x*facteur), m) for x in tmp]
        # if there remain some genes to put in chromosomes ...
        lengths[-1] += s-sum(lengths)
        if lengths[-1] >= m:
            return lengths
