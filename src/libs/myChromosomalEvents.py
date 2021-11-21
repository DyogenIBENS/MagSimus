# -*- coding: utf-8 -*-

import collections
import random
import sys

import myEvolutionProbas
import myMagSimusTools

# Turn over the chromosomal segment if needed
def applyStrand(l, strand):
    print "with strand", strand,
    return l if strand > 0 else [(gene, -s) for (gene, s) in reversed(l)]


# Chromosomal rearrangements
def performChromEvents(newGenome, nbThEvents, ChromSegVM_Mean, ChromSegVM_Kappa,
                       idx2eventType, translocSamplingAlgo, lowerCap, upperCap,
                       fissionSamplingAlgo, forceChrNumber,
                       printIntermGenomes=False, prefix=''):
    # ChromSegVM_Mean = arguments['chr:invertVonMisesMean'];
    # ChromSegVM_Kappa = arguments['chr:invertVonMisesKappa'];
    # translocSamplingAlgo = arguments['chr:translocSamplingAlgo'];
    # lowerCap = arguments['chr:sizeLowerCap'];
    # upperCap = arguments['chr:sizeUpperCap'];
    # fissionSamplingAlgo = arguments['chr:fissionSamplingAlgo'];
    # forceChrNumber = arguments['chr:forceChrNumber']

    assert set(["chrInvert", "chrTransloc", "chrFusion", "chrFission"]) == \
        set(idx2eventType.inv.keys())
    # idx2NbEffEvents = collections.defaultdict(int)
    nbEffEvents = {i: 0 for i in set(idx2eventType.inv.keys())}
    # niniC = len(newGenome)

    if printIntermGenomes:
        myMagSimusTools.buildFamId(newGenome)
        myMagSimusTools.printCurrentGenome(newGenome, prefix=prefix)

    # nbThEvents: nb theoretical events
    # nbThEvents contains 4 int values, each value corresponds to a number of
    # events for a specific chromosomal rearrangement.

    # List of events to perform
    # 'events' is a list of ints in [0,3], the number of each int 'i' is equal to
    # nbThEvents[i].
    events = []
    for (i, evt) in idx2eventType.iteritems():
        events.extend([i] * nbThEvents[evt])
    random.shuffle(events)
    reinsertEvtCnt = 0
    numberOfTriedReinsertionsBeforeBreak = 1000
    while len(events) > 0:
        idEvt = events.pop(0)

        print prefix, "rearrangement", idx2eventType[idEvt],
        # It is an inversion
        if idEvt == idx2eventType.inv['chrInvert']:
            # Inverted chromosomal segment
            (c, x1, x2) = myEvolutionProbas.randomSlice(newGenome, ChromSegVM_Mean, ChromSegVM_Kappa)
            print "%d:%d=(%s)-%d=(%s)" % (c+1,
                                          x1, myMagSimusTools.getIntervStr(newGenome, c, x1),
                                          x2, myMagSimusTools.getIntervStr(newGenome, c, x2)),
            # x1 and x2 correspond to extremal intergenes of the rearranged
            # chromosomal segment. x1 and x2 are ints in [0,len(chr)].
            # 0 and len(chr) correspond to extremities.
            #   * x1=0 corresponds to the left extremity of chr.
            #   * x1=1 corresponds to the intergene between the 1st gene and the 2nd
            # gene.
            #   * x1=len(chr) corresponds the right extremity of chr.
            # New chromosome with the inverted segment in the middle
            newGenome[c][x1:x2] = applyStrand(newGenome[c][x1:x2], -1)

        # It is a translocation
        # TODO : change the functioning of the translocation, for
        # the moment it is not really a 'translocation' but a 'transposition'.
        elif idEvt == idx2eventType.inv['chrTransloc']:
            if len(newGenome) < 2:
                # If there are not enough chromosomes, the fusion is
                # postponed, i.e. added at the end of the event vector.
                # Underlying Hypothesis: there may be fissions that will yield
                # more chromosomes before applying the last fission.
                if all([e == idx2eventType.inv['chrFusion']
                        or e==idx2eventType.inv['chrTransloc'] for e in events]):
                    print(events)
                    raise ValueError('too many fusions. Not enough chromosomes')
                if forceChrNumber:
                    events.append(idx2eventType.inv["chrFusion"])
                    reinsertEvtCnt += 1
                continue
            else:
                # The translocated chromosomal segment
                # A destination location is chosen
                # FIXME : The choice of the destination chromosome should be
                # independent of the size of the chromosome. This might match better
                # with known mechanisms of translocations.

                # Choose two (different chromosomes) ...
                if (translocSamplingAlgo == "lengthPropSampling"):
                    # ... based on their lengths
                    c1 = myEvolutionProbas.randomChromosome(newGenome, minNbGenes=2)
                    c2 = myEvolutionProbas.randomChromosome(newGenome[:c1] + newGenome[(c1+1):], minNbGenes=2)
                elif (translocSamplingAlgo == "simpleRandomSampling"):
                    # ... independently of their lengths
                    (c1, c2) = random.sample(range(len(newGenome)), 2)
                else:
                    raise ValueError('Wrong sampling algorithm for Translocations. Please choose "lengthPropSampling" or "simpleRandomSampling".')

                # if (c2 >= c1):
                #     c2 += 1
                b1 = random.randint(0, len(newGenome[c1]))
                b2 = random.randint(0, len(newGenome[c2]))
                sizesAreGood = (lowerCap <= b1 + len(newGenome[c2]) - b2 and upperCap > b1 + len(newGenome[c2]) - b2
                                and lowerCap <= b2 + len(newGenome[c1]) - b1 and upperCap > b2 + len(newGenome[c1]) - b1)

                if sizesAreGood:
                    print(
                        "Exchange of %d:(0:%d) with intergenes %s and %d:(0:%d) with intergenes %s" %
                        (
                            c1, b1, myMagSimusTools.getIntervStr(newGenome, c1, b1), c2, b2,
                            myMagSimusTools.getIntervStr(newGenome, c2, b2)))

                    r1 = newGenome[c1][:b1]
                    del newGenome[c1][:b1]
                    r2 = newGenome[c2][:b2]
                    del newGenome[c2][:b2]
                    newGenome[c2][:0] = applyStrand(r1, myEvolutionProbas.randomStrand())
                    newGenome[c1][:0] = applyStrand(r2, myEvolutionProbas.randomStrand())
                elif forceChrNumber:
                    reinsertEvtCnt += 1
                    if reinsertEvtCnt > numberOfTriedReinsertionsBeforeBreak and \
                        all([e == idx2eventType.inv['chrFusion']
                             or e == idx2eventType.inv['chrTransloc'] for e in events]):
                        raise ValueError(
                            'too many Fusions, all chromosomes are over upperCap.')
                    elif reinsertEvtCnt > numberOfTriedReinsertionsBeforeBreak and \
                            all([e == idx2eventType.inv['chrFission']
                                 or e == idx2eventType.inv['chrTransloc'] for e in events]):
                        raise ValueError(
                            'too many Fissions, all chromosomes are under lowerCap.')
                    events.append(idx2eventType.inv["chrTransloc"])

        # It is a fusion
        # FIXME This fusion chooses a random chr independently of the chromosome
        # length. The probability to choose small chromosomes is as high as the
        # probability to choose long chromosomes.
        # This may be the cause of the exponential distribution of the lengths
        # of chromosomes at the end of simulation, in extants genomes.
        elif idEvt == idx2eventType.inv['chrFusion']:
            if len(newGenome) < 2:
                # If there are not enough chromosomes, the fusion is
                # postponed, i.e. added at the end of the event vector.
                # Underlying Hypothesis: there may be fissions that will yield
                # more chromosomes before applying the last fission.
                if all([e == idx2eventType.inv['chrFusion']
                        or e == idx2eventType.inv['chrTransloc'] for e in events]):
                    print(events)
                    raise ValueError('too many fusions. Not enough chromosomes')
                if forceChrNumber:
                    events.append(idx2eventType.inv["chrFusion"])
                    reinsertEvtCnt += 1
                continue
            else:
                # Two chromosomes that will be fused
                (c1, c2) = random.sample(range(len(newGenome)), 2)
                print c1+1,
                cc1 = applyStrand(newGenome[c1], myEvolutionProbas.randomStrand())
                print "+", c2+1,
                cc2 = applyStrand(newGenome[c2], myEvolutionProbas.randomStrand())
                # both chromosomes are put end to end, they may also be
                # inversed
                newChromosome = cc1 + cc2
                sizesAreGood = (upperCap > len(newChromosome))
                if sizesAreGood:
                    newGenome[c1] = newChromosome
                    #newGenome[c2] = [] #FIXME Joseph commented this line
                    # WARNING the list newGenome is changed, indices won't be conserved.
                    # Here we can del newGenome[c2] because their is no need to conserve
                    # the indices signification in the loop on chromosomic rearrangements
                    # FIXME Joseph uncommented this line
                    del newGenome[c2]
                elif forceChrNumber:
                    if reinsertEvtCnt > numberOfTriedReinsertionsBeforeBreak and \
                            all([e == idx2eventType.inv['chrFusion']
                                 or e == idx2eventType.inv['chrTransloc'] for e in events]):
                        raise ValueError(
                            'too many Fusions, all chromosomes are over upperCap.')
                    events.append(idx2eventType.inv["chrFusion"])
                    reinsertEvtCnt += 1

        # It is a fission/break
        elif idEvt == idx2eventType.inv['chrFission']:

            # Choose chromosome ...
            if (fissionSamplingAlgo == "lengthPropSampling"):
                # ... based on its length
                (c, x) = myEvolutionProbas.randomIntergene(newGenome,
                                                           removedExtremityLength=1)
            elif (fissionSamplingAlgo == "simpleRandomSampling"):
                # ... independently of its length
                if all([len(e) < 2 for e in newGenome]):
                    raise ValueError('Not enough genes to do Fissions.')
                c = random.sample(range(len(newGenome)), 1)[0]
                if (len(newGenome[c])) > 1:
                    x = random.randint(1, len(newGenome[c]) - 1)
                else:
                    while len(newGenome[c]) < 2:
                        c = random.sample(range(len(newGenome)), 1)[0]
                        if (len(newGenome[c])) > 1:
                            x = random.randint(1, len(newGenome[c]) - 1)
                            break

            else:
                raise ValueError('Wrong sampling algorithm for Fissions. Please choose "lengthPropSampling" or "simpleRandomSampling".')

            sizesAreGood = (lowerCap <= x and lowerCap <= (len(newGenome[c]) - x))
            if sizesAreGood:
                # FIXME exclude the print outside of the function (cf the translocation and
                # fission functions)
                intervBreakStr = "%d:%d=(%s)" % (
                    c + 1, x, myMagSimusTools.getIntervStr(newGenome, c, x))
                newGenome.append(newGenome[c][:x])
                del newGenome[c][:x]

                print >> sys.stdout, intervBreakStr,
                print "to", len(newGenome),
            elif forceChrNumber:
                if reinsertEvtCnt > numberOfTriedReinsertionsBeforeBreak and \
                    all([e == 4 or e == 2 for e in events]):
                    raise ValueError('too many Fissions, all chromosomes are under lowerCap.')
                events.append(idx2eventType.inv["chrFission"])
                reinsertEvtCnt += 1

        print
        if printIntermGenomes:
            myMagSimusTools.printCurrentGenome(newGenome, prefix=prefix)

        nbEffEvents[idx2eventType[idEvt]] += 1

    # nbEffEvents = {}
    # for idEvt, nbEvt in idx2NbEffEvents.iteritems():
    #     nbEffEvents[idx2eventType[idEvt]] = nbEvt

    # DEBUG assertion
    # Conservation law of chromosomes
    #if len(newGenome) >= 2:
    #    assert len(newGenome)-niniC == nbBreak - nbFus, "%s - %s == %s - %s" % (len(newGenome), niniC, nbBreak, nbFus)

    return (newGenome, nbEffEvents)