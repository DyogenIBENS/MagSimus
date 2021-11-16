# -*- coding: utf-8 -*-

import os
import sys
import random
import collections
import myEvents

class ManagerOfEventsIterator:
    def __init__(self, moe):
        assert isinstance(moe, ManagerOfEvents)
        self.moe = moe

    def __iter__(self):
        return self

    def next(self):
        if len(self.moe) <= 0:
            raise StopIteration
        else:
            return self.moe.popleft()

class ManagerOfEvents(collections.deque):
    def __init__(self, paramValues, idx2eventType,
                 maxReinsertion=300):
        self.idx2eventType = idx2eventType
        lIdEvtNbEvts = [(idEvt, paramValues[evt]) for (idEvt, evt) in idx2eventType.iteritems()]
        maxSizeOfQueue = sum([nbEvts for (_, nbEvts) in lIdEvtNbEvts])
        listOfEvents = []
        # 'listOfEvents' is a list of ints in [0,3], the number of each int 'i' is equal to
        # paramValues[i].
        for (idEvt, nbEvts) in lIdEvtNbEvts:
            listOfEvents.extend([idEvt] * nbEvts)
        # random.shuffle is faster on lists
        random.shuffle(listOfEvents)
        collections.deque.__init__(self, maxlen=maxSizeOfQueue)
        for idEvt in listOfEvents:
            self.append(idEvt)
        del listOfEvents
        self.maxReinsertions = maxReinsertion
        self.reinsertEvtCpt = collections.defaultdict(int)
        self.chrInvertId = idx2eventType.inv['chrInvert']
        self.chrTranslocId = idx2eventType.inv['chrTransloc']
        self.chrFusionId = idx2eventType.inv['chrFusion']
        self.chrFissionId = idx2eventType.inv['chrFission']
        self.evtsIncreasingChromLengths = [self.chrFusionId]
        self.evtsDecreasingChromLengths = [self.chrFissionId]
        try:
            # For the usage in MagSimus1
            # In MagSimus 2, there are no geneBirth, geneDup, ...; hence this code will raise an expected exception
            # TODO
            # It would be cleaner to have an option "MagSimus 1 mode" instead of this default exception...
            self.geneBirthId = idx2eventType.inv['geneBirth']
            self.geneDupId = idx2eventType.inv['geneDup']
            self.geneLossId = idx2eventType.inv['geneLoss']
            # create two new ids for geneTandemDup and geneDispersedDup
            self.evtsIncreasingChromLengths.extend([self.geneDupId, self.geneBirth])
            self.evtsDecreasingChromLengths.extend([self.geneLossId])
        except:
            pass

    def __iter__(self):
        return ManagerOfEventsIterator(self)

    def _reinsertFiniteNbOfTimes(self, idEvt):
        # numberOfTriedReinsertionsBeforeBreak is necessary if, for
        # instance if all chromosomes are of size 99, with a upper
        # cap 100, finding the very few translocations that yield ok
        # chromsomes would take too long.
        self.reinsertEvtCpt[idEvt] += 1
        print >> sys.stderr, "event idEvt=%s already reinserted %s times" % (idEvt, self.reinsertEvtCpt[idEvt])
        if self.reinsertEvtCpt[idEvt] > self.maxReinsertions:
            if idEvt == self.chrTranslocId:
                raise ValueError(
                    'A translocation yielding two chromosomes with (lowerCap <= lengths < upperCap) is too unprobable.'
                    '\nPlease check the overall ratio fusions/fissions that may be unbalanced, leading to either,'
                    '\ntoo many small chromosomes or too many long chromosomes.')
            elif idEvt == self.chrFusionId:
                raise ValueError(
                    'A fusion yielding one chromosome with length <= upperCap is too unprobable.')
            elif idEvt == self.chrFissionId:
                raise ValueError(
                    'A fission yielding two chromosomes with lowerCap < lengths is too unprobable.')
        else:
            self.append(idEvt)

    def _handleExceptionEvt(self, idEvt, e, genome, cptRecursion=0, **kwargs):
        #cptRecursion = kwargs.get('cptRecursion', 0)
        flag = e.__class__.__name__ + ('' if cptRecursion == 0 else "/cptRecursion=%s" % cptRecursion)
        infos = None
        if isinstance(e, myEvents.LessThanTwoChromosomes):
            # If there are not enough chromosomes, the event is
            # postponed, i.e. added at the end of the event vector.
            # Underlying Hypothesis: there may be fissions that will yield
            # more chromosomes before applying this event.
            if any([evtId == self.chrFissionId for evtId in self]):
                flag = flag + '/ReinsertEventAfterCompensatingEvt'
                self.append(idEvt)
            else:
                raise ValueError(
                    "Unable to perform %s event, only one chromosome in the genome." % self.idx2eventType[idEvt]\
                    + "\nToo many fusions and not enough fissions on this branch.")
        elif isinstance(e, myEvents.NoChromosomeWithAtLeastTwoGenes):
            # If all chromosomes only contain one gene, the event (fission, inversion...) is
            # postponed, i.e. added at the end of the event vector.
            # Underlying Hypothesis: there may be (fusions, duplications, gene births, ...) that will yield
            # longer chromosomes before applying the last fission.
            if any([e in self.evtsIncreasingChromLengths for e in self]):
                flag = flag + '/ReinsertEventAfterCompensatingEvt'
                self.append(idEvt)
            else:
                print sys.stderr, [len(chrom) for chrom in genome.values()]
                raise ValueError("Unable to perform %s event, all chromosomes contain only one gene." % self.idx2eventType[idEvt]\
                                 + "\nToo many %s and not enough %s on this branch." %\
                                 ([self.idx2eventType[idEvt] for idEvt in self.evtsDecreasingChromLengths],
                                  [self.idx2eventType[idEvt] for idEvt in self.evtsIncreasingChromLengths]))
        elif isinstance(e, myEvents.ChromosomeTooLong):
            try:
                # try some other times... until cptRecursion == recursionMax
                (infos, flag) = self.handleEvent(idEvt, genome, cptRecursion=cptRecursion+1, **kwargs)
            except myEvents.MagSimusException as e:
                if isinstance(e, myEvents.TooMuchTime):
                    if any([(e in self.evtsDecreasingChromLengths) for e in self]):
                        flag = flag + '/ReinsertEventAfterCompensatingEvt'
                        self.append(idEvt)
                    else:
                        flag = flag + '/ReinsertEventAfterNoCompensatingEvt'
                        self._reinsertFiniteNbOfTimes(idEvt)
                else:
                    self._handleExceptionEvt(idEvt, e, genome, cptRecursion=cptRecursion+1, **kwargs)
        elif isinstance(e, myEvents.ChromosomeTooShort):
            try:
                # try some other times... until cptRecursion == recursionMax
                (infos, flag) = self.handleEvent(idEvt, genome, cptRecursion=cptRecursion+1, **kwargs)
            except myEvents.MagSimusException as e:
                if isinstance(e, myEvents.TooMuchTime):
                    if any([(e in self.evtsIncreasingChromLengths) for e in self]):
                        flag = flag + '/ReinsertEventAfterCompensatingEvt'
                        self.append(idEvt)
                    else:
                        flag = flag + '/ReinsertEventAfterNoCompensatingEvt'
                        self._reinsertFiniteNbOfTimes(idEvt)
                else:
                    self._handleExceptionEvt(idEvt, e, genome, cptRecursion=cptRecursion+1, **kwargs)
        elif isinstance(e, myEvents.TooMuchTime):
            raise myEvents.TooMuchTime("Events %s takes to much time to be prepared" % self.idx2eventType[idEvt])
        else:
            raise e
        return (infos, flag)

    def handleEvent(self, idEvt, genome, cptRecursion=0, **kwargs):
        if cptRecursion == 100:
            raise myEvents.TooMuchTime
        try:
            if idEvt == self.geneBirthId:
                (infos, flag) = myEvents.prepareGeneBirth(genome, **kwargs)
            elif idEvt == self.geneDupId:
                (infos, flag) = myEvents.prepareGeneDuplication(genome, **kwargs)
                assert flag in ['TandemDuplication', 'DispersedDuplication']
            elif idEvt == self.geneLossId:
                (infos, flag) = myEvents.prepareGeneLoss(genome, **kwargs)
                assert flag in ['OK', 'RemoveEmptyChromosome']
            elif idEvt == self.chrFissionId:
                (infos, flag) = myEvents.prepareFission(genome, **kwargs)
            elif idEvt == self.chrFusionId:
                (infos, flag) = myEvents.prepareFusion(genome, **kwargs)
            elif idEvt == self.chrInvertId:
                (infos, flag) = myEvents.prepareInversion(genome, **kwargs)
            elif idEvt == self.chrTranslocId:
                (infos, flag) = \
                    myEvents.prepareReciprocalTranslocation(genome, **kwargs)
            else:
                raise ValueError("idEvt %s unknown. idEvt should be in %s" % (idEvt, self.idx2eventType.inv.keys()))
        except myEvents.MagSimusException as e:
            # kwargs['cptRecursion'] = cptRecursion
            # print >> sys.stderr, type(e)
            (infos, flag) = self._handleExceptionEvt(idEvt, e, genome, cptRecursion=cptRecursion, **kwargs)
        return (infos, flag)
