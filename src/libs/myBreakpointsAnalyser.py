import collections
import copy
import sys
from libs import myEvents
import itertools
from utils import myTools, myMaths
from utils import myLightGenomes
from utils import myFile
import utils.myMapping as myMapping
from utils.myLightGenomes import LightGenome, Families
import utils.myIntervals as myIntervals
from utils.myIntervals import OAdjacency


# These breakpoint statistics must consider the possibility that two events separated in times may create smaller
# synteny blocks if two of there breakpoints fall on the same chromosome.
class BreakpointsAnalyser():
    # when an old adjecency is recovered and also deleting genes when a gene loss occur.
    def __init__(self, iniAncGenome, isLazy=False):
        # if isLazy, the breakpoint analyser is 'lazy' and do nothing when it is called
        self.lazy = isLazy
        if self.lazy:
            return
        isinstance(iniAncGenome, LightGenome)
        # This copy of the remainingAncestralGenome will record breaks and recoveredAdjs for each event
        # FIXME check that the copy is enough deep to not disturb the simulation
        self.conservedAncestralBlocks = copy.copy(iniAncGenome)
        # This copy won't evolve and is used to keep a record of the remainingAncestralGenome
        self.remainingAncestralGenome = copy.copy(iniAncGenome)
        self.setGenesAtAncChromExtremities = set()
        self.setAncChromExtremities = set()
        for (c, chrom) in iniAncGenome.iteritems():
            self.conservedAncestralBlocks[c] = copy.copy(iniAncGenome[c])
            self.remainingAncestralGenome[c] = copy.copy(iniAncGenome[c])
            g_beg = chrom[0]
            g_end = chrom[-1]
            self.setGenesAtAncChromExtremities.update({g_beg.n, g_end.n})
            ht_beg = 't' if g_beg.s == +1 else 'h'
            self.setAncChromExtremities.add((g_beg.n, ht_beg))
            ht_end = 'h' if g_end.s == +1 else 't'
            self.setAncChromExtremities.add((g_end.n, ht_end))

        self.remainingAncestralGenome.computeDictG2P()
        # they have the same dict g2p at the beginning
        self.conservedAncestralBlocks.g2p = copy.copy(self.remainingAncestralGenome.g2p)
        self.conservedAncestralBlocks.withDict = True
        # self.extremitiesOfConservedAncestralBlocks = set([gExtr for chrom in self.conservedAncestralBlocks
        #                                                  for gExtr in (myIntervals.geneExtremityFromGene(chrom[0], -1),
        #                                                                myIntervals.geneExtremityFromGene(chrom[-1], +1))])
        self.setGeneExtrOfConservedAncBlocksWithoutAncChromExtr = set()

        self.nbBreakpointsAtChromExtremity = 0
        self.setBrokenGeneExtrs = set()
        self.nbBreakpoints = 0
        self.nbBreakpointsReuse = 0
        self.nbInnerBreakpointsReuse = 0
        self.setOnceDisruptedAncAdjs = set()
        self.nbStandardInnerBreakpoints = 0
        # nb of recovered previously disrupted adjacencies
        self.nbRecoveredAncAdjs = 0
        self.setRecoveredAncAdjacencies = set()
        self.nbBreakpointsWithinRecoveredAdjacencies = 0
        self.nbAncChromsRemovedBecauseOfGeneLoss = 0


    #######################
    # Analyse breakpoints #
    #######################
    def recordBreakpoint(self, (c, x), genomeBeforeBreakpoint, ignoredGeneNames=set()):
        if self.lazy:
            return
        isBreakpointReuse = False

        assert isinstance(genomeBeforeBreakpoint, LightGenome)
        chrom = genomeBeforeBreakpoint[c]
        assert 0 <= x <= len(chrom)
        (idxGl, idxGr) = myIntervals.findConsideredFlankingGenesIdxs(chrom, x, ignoredGeneNames=ignoredGeneNames)
        if idxGl is None and idxGr is None:
            # this breakpoint only impact a chromosome full of ignoredGeneNames
            raise ValueError('This should not occur in our specific case since each rearrangement '
                             'should at least rearrange one ancestral gene')
        elif idxGl is None or idxGr is None:
            # 1) the breakpoint is at one extremity (of the chrom filtered to remove ignoredGeneNames)
            self.nbBreakpoints += 1
            self.nbBreakpointsAtChromExtremity += 1
            if idxGl is not None:
                assert 0 <= idxGl
                gl = chrom[idxGl]
                glht = 'h' if gl.s == +1 else 't'
                brokenIntergene = (gl.n, glht)
            else:
                assert idxGr is not None and idxGr <= len(chrom)
                gr = chrom[idxGr]
                grht = 't' if gr.s == +1 else 'h'
                brokenIntergene = (gr.n, grht)
            assert brokenIntergene[0] not in ignoredGeneNames
            if brokenIntergene in self.setBrokenGeneExtrs:
                self.nbBreakpointsReuse += 1
                isBreakpointReuse = True
            else:
                self.setBrokenGeneExtrs.add(brokenIntergene)
        else:
            # the breakpoint falls within an inner intergene
            self.nbBreakpoints += 1
            # gene left
            gl = chrom[idxGl]
            # gene right
            gr = chrom[idxGr]
            assert gl.n not in ignoredGeneNames and gr.n not in ignoredGeneNames
            adjacency = OAdjacency(gl, gr)
            # brokenGeneExtr1 = (g1n, g1ht)
            # brokenGeneExtr2 = (g2n, g2ht)
            (brokenGeneExtr1, brokenGeneExtr2) = myIntervals.geneExtremitiesFromAdjacency(adjacency)
            # Even if the two are equivalent in genomeBeforeBreakpoint, once the rearrangement is made we need to keep
            # recorded the two broken gene extremities for next analysis for next breakpoints (breakpoints reuse...)
            if brokenGeneExtr1 in self.setBrokenGeneExtrs or brokenGeneExtr2 in self.setBrokenGeneExtrs:
                # 2.a) The breakpoint falls in a previously broken inner intergene. (inner <=>  not a chrom extremity)
                self.nbBreakpointsReuse += 1
                isBreakpointReuse = True
                self.nbInnerBreakpointsReuse += 1
                if adjacency in self.setRecoveredAncAdjacencies:
                    # print >> sys.stderr, "breakpoint within recovered OAdjacency = %s" % str(adjacency)
                    self.nbBreakpointsWithinRecoveredAdjacencies += 1
                    self.setRecoveredAncAdjacencies.remove(adjacency)
                    self._tearApartConservedAncestralBlocks(adjacency)
            elif brokenGeneExtr1 in self.setAncChromExtremities or brokenGeneExtr2 in self.setAncChromExtremities:
                # 2.b) The breakpoint falls in an intergene newly created without being due to a breakpoint,
                # this can happen after a fusion for instance. The intergene is then flanked by at least one intergene that
                # was at a chrom extremity and is now an inner chrom intergene. In the case of a reciprocal translocation, a
                # chromosome extremity might become an inner intergene but it will be adjacent to a broken intergene,
                # thus this last case will be catched in 2.a).
                # set |= {item1, item2} <=> set.update({item1, item2})
                # FIXME
                # This is counted as a breakpoint reuse even if this might be debatable
                self.nbBreakpointsReuse += 1
                isBreakpointReuse = True
                self.nbInnerBreakpointsReuse += 1
                self.setBrokenGeneExtrs |= {brokenGeneExtr1, brokenGeneExtr2}
            else:
                # 3) The breakpoint falls within a not yet broken (or not newly created, cf 2.b)) inner intergene
                # A standard inner breakpoint is neither a breakpoint reuse or a breakpoint at a chromosomal extremity
                self.nbStandardInnerBreakpoints += 1
                self.setBrokenGeneExtrs |= {brokenGeneExtr1, brokenGeneExtr2}
                # Ancestral adjacencies that have been disrupted at least once (they may be recovered afterwards)
                self.setOnceDisruptedAncAdjs.add(adjacency)
                self._tearApartConservedAncestralBlocks(adjacency)
        return isBreakpointReuse

    # single_leading_underscore: weak "internal use" indicator. E.g. "from M import *" does not import objects whose
    # name starts with an underscore.
    def _tearApartConservedAncestralBlocks(self, adjacency):
        isinstance(self.conservedAncestralBlocks, LightGenome)
        brokenIntergene = myIntervals.intergeneFromAdjacency(adjacency, self.conservedAncestralBlocks, default=None)
        if brokenIntergene is not None:
            (c, x) = brokenIntergene
            # c -> c=c[:x] & newC=c[x:]
            self.conservedAncestralBlocks = myEvents.performFission(self.conservedAncestralBlocks, (c, x), updateGenomeDict=True)
            (ogl, ogr) = adjacency
            glExtr = myIntervals.geneExtremityFromGene(ogl, +1)
            grExtr = myIntervals.geneExtremityFromGene(ogr, -1)
            self.setGeneExtrOfConservedAncBlocksWithoutAncChromExtr |= {glExtr, grExtr}
        else:
            return

    def _fuseConservedAncestralBlocks(self, newAdjacency):
        isinstance(self.conservedAncestralBlocks, LightGenome)
        (ogl, ogr) = newAdjacency
        ag1pos = self.conservedAncestralBlocks.getPosition(ogl.n)
        # check that g1 is at one block extremity
        assert ag1pos.idx in {0, len(self.conservedAncestralBlocks[ag1pos.c]) - 1}, str(ag1pos.idx)
        c1 = ag1pos.c
        ag2pos = self.conservedAncestralBlocks.g2p[ogr.n]
        # check that g2 is at one block extremity
        assert ag2pos.idx in {0, len(self.conservedAncestralBlocks[ag2pos.c]) - 1}, str(ag2pos.idx)
        c2 = ag2pos.c

        ag1 = self.conservedAncestralBlocks[c1][ag1pos.idx]
        ag2 = self.conservedAncestralBlocks[c2][ag2pos.idx]

        if len(self.conservedAncestralBlocks[c1]) == 1:
            # only one gene in the conservedAncestralBlock containing g1
            s1 = +1 if ag1.s == ogl.s else -1
        elif ag1pos.idx == 0:
            assert (ag1.s is None and ogl.s is None) or ag1.s == - ogl.s
            s1 = -1
        else:
            assert ag1pos.idx == len(self.conservedAncestralBlocks[c1]) - 1
            assert (ag1.s is None and ogl.s is None) or ag1.s == ogl.s
            s1 = +1

        if len(self.conservedAncestralBlocks[c2]) == 1:
            # only one gene in the conservedAncestralBlock containing g1
            s2 = +1 if ag2.s == ogr.s else -1
        elif ag2pos.idx == 0:
            assert (ag2.s is None and ogr.s is None) or ag2.s == ogr.s
            s2 = +1
        else:
            assert ag2pos.idx == len(self.conservedAncestralBlocks[c2]) - 1
            assert (ag2.s is None and ogr.s is None) or ag2.s == - ogr.s
            s2 = -1
        # c1 & c2 -> c1 = [s1*c1 + s2*c2]
        self.conservedAncestralBlocks = myEvents.performFusion(self.conservedAncestralBlocks, ((c1, s1), (c2, s2)), updateGenomeDict=True)
        glExtr = myIntervals.geneExtremityFromGene(ogl, +1)
        grExtr = myIntervals.geneExtremityFromGene(ogr, -1)
        self.setGeneExtrOfConservedAncBlocksWithoutAncChromExtr -= {glExtr, grExtr}

        ###########################
    # Analyse new adjacencies #
    ###########################
    # replace lost gene 'ogLost' by its neighbours 'ogl' and 'ogr' (oriented gene left and oriented gene right resp.)
    def _replaceGLostByGNeighboursInAdjs(self, setAdjs, (ogl, ogLost, ogr)):
        newSetAdjs = set()
        for adjacency in setAdjs:
            assert isinstance(adjacency, OAdjacency)
            (og1, og2) = adjacency
            # bool(a) ^ bool(b) <=> bool(a) xor bool(b) (exclusive or)
            # <=> bool(a) or bool(b) and not (bool(a) and bool(b))
            assert og1.n != og2.n
            # print >> sys.stderr, "[%s, %s, %s]" % (str(ogl), str(ogLost), str(ogr)),  "adjacency = %s" % str(adjacency)
            if og1.n == ogLost.n:
                if og1.s == ogLost.s:
                    # og1.s == ogLost.s means that the adjacency is coherent with the orientation of the deleted gene
                    # before deletion.
                    # ogl is None <=> ogLost was at one chrom extremity
                    newSetAdjs.update({OAdjacency(ogl, og2)} if ogl is not None and og2.n != ogl.n else {})
                else:
                    assert og1.s == - ogLost.s
                    # og1.s == - ogLost.s means that the adjacency is reversed compared to the orientation of the
                    # deleted gene before deletion.
                    # ogr is None <=> ogLost was at one chrom extremity
                    newSetAdjs.update({OAdjacency(og2.reversed(), ogr)} if ogr is not None and og2.n != ogr.n else {})
            elif og2.n == ogLost.n:
                if og2.s == ogLost.s:
                    # same principle as before
                    newSetAdjs.update({OAdjacency(og1, ogr)} if ogr is not None and og1.n != ogr.n else {})
                else:
                    # same principle as before
                    assert og2.s == - ogLost.s
                    newSetAdjs.update({OAdjacency(ogl, og1.reversed())} if ogl is not None and og1.n != ogl.n else {})
            else:
                # No need to modify this adjacency, we keep it as it is
                newSetAdjs.add(adjacency)
        return newSetAdjs

    # replace lost gene 'ogLost' by its neighbours 'ogl' and 'ogr' (oriented gene left and oriented gene right resp.)
    def _replaceGLostByGNeighboursInGeneExtrs(self, setGeneExtremities, (ogl, ogLost, ogr)):
        newSetGeneExtrs = set()
        for geneExtr in setGeneExtremities:
            (gn, ht) = geneExtr
            if gn == ogLost.n:
                if ogLost.s == +1:
                    if ht == 'h':
                        # the removal of gLost translates the region from the head to the tail of gLost,
                        # towards the left, since ogLost.s == +1.
                        newSetGeneExtrs.update({myIntervals.geneExtremityFromGene(ogl, +1)} if ogl is not None else {})
                    else:
                        assert ht == 't'
                        # the removal of gLost translates the region from the tail to the head of gLost,
                        # towards the right, ogLost.s == +1.
                        newSetGeneExtrs.update({myIntervals.geneExtremityFromGene(ogr, -1)} if ogr is not None else {})
                elif ogLost.s == -1:
                    if ht == 'h':
                        # the removal of gLost translates the region from the head to the tail of gLost,
                        # towards the right, since ogLost.s == -1.
                        newSetGeneExtrs.update({myIntervals.geneExtremityFromGene(ogr, -1)} if ogr is not None else {})
                    else:
                        assert ht == 't'
                        # the removal of gLost translates the region from the tail to the head of gLost,
                        # towards the left, ogLost.s == -1.
                        newSetGeneExtrs.update({myIntervals.geneExtremityFromGene(ogl, +1)} if ogl is not None else {})
            else:
                newSetGeneExtrs.add(geneExtr)
        return newSetGeneExtrs


    def updateGeneContentAfterGeneLoss(self, gLostN):
        if self.lazy:
            return
        # update gene content in self.conservedAncestralBlocks (remove the deleted anc gene)
        gLostPosCAB = self.conservedAncestralBlocks.getPosition(gLostN, default=None)
        assert gLostPosCAB is not None
        if len(self.conservedAncestralBlocks[gLostPosCAB.c]) == 1:
            self.nbAncChromsRemovedBecauseOfGeneLoss += 1
        self.conservedAncestralBlocks = myEvents.performGeneLoss(self.conservedAncestralBlocks, gLostPosCAB, updateGenomeDict=True)

        # update gene content in self.remainingAncestralGenome (remove the deleted anc gene)
        gLostPosRAG = self.remainingAncestralGenome.getPosition(gLostN, default=None)
        assert gLostPosRAG is not None
        #print >> sys.stderr, "len(g[%s])=%s, g[%s][%s]" % (gLostPosRAG.c, len(self.remainingAncestralGenome[gLostPosRAG.c]), gLostPosRAG.c, gLostPosRAG.idx)
        ogLost = self.remainingAncestralGenome[gLostPosRAG.c][gLostPosRAG.idx]
        self.remainingAncestralGenome = myEvents.performGeneLoss(self.remainingAncestralGenome, gLostPosRAG, updateGenomeDict=True)
        if gLostPosRAG.c in self.remainingAncestralGenome:
            self.recordNewAdjacencyFromNewIntergene(self.remainingAncestralGenome, gLostPosRAG, ignoredGeneNames=set())

        chrom = self.remainingAncestralGenome[gLostPosRAG.c]
        # index of the new intergene in self.remainingAncestralGenome after the gene deletion
        x = gLostPosRAG.idx
        # no need of ignoredGeneNames since self.remainingAncestralGenome contains only considered anc genes
        (idxGl, idxGr) = myIntervals.findConsideredFlankingGenesIdxs(chrom, x, ignoredGeneNames=set())
        (ogl, ogr) = (None, None)
        if idxGl is not None:
            ogl = chrom[idxGl]
        if idxGr is not None:
            ogr = chrom[idxGr]
        # if ogl is not None and ogr is not None:
        #     newAdj = OAdjacency(ogl, ogr)
        #     self.recordNewAdjacency(newAdj)

        # listOfSetsOfGeneExtrs = [self.setAncChromExtremities, self.setBrokenGeneExtrs]
        # listOfSetsOfAdjs = [self.setRecoveredAncAdjacencies, self.setOnceDisruptedAncAdjs]
        self.setBrokenGeneExtrs = self._replaceGLostByGNeighboursInGeneExtrs(self.setBrokenGeneExtrs, (ogl, ogLost, ogr))
        self.setAncChromExtremities = self._replaceGLostByGNeighboursInGeneExtrs(self.setAncChromExtremities, (ogl, ogLost, ogr))
        self.setGeneExtrOfConservedAncBlocksWithoutAncChromExtr = self._replaceGLostByGNeighboursInGeneExtrs(self.setGeneExtrOfConservedAncBlocksWithoutAncChromExtr, (ogl, ogLost, ogr))

        self.setRecoveredAncAdjacencies = self._replaceGLostByGNeighboursInAdjs(self.setRecoveredAncAdjacencies, (ogl, ogLost, ogr))
        self.setOnceDisruptedAncAdjs = self._replaceGLostByGNeighboursInAdjs(self.setOnceDisruptedAncAdjs, (ogl, ogLost, ogr))
        # FIXME genes at extremities of chromosomes should not be here
        # FIXME error here, put some adjacencies that corresponds to breakpoints at chrom extremities

    def recordNewAdjacencyFromNewIntergene(self, genomeWithNewAdj, (c, x), ignoredGeneNames=set()):
        if self.lazy:
            return
        assert isinstance(genomeWithNewAdj, LightGenome)
        if c in genomeWithNewAdj:
            chrom = genomeWithNewAdj[c]
        else:
            # FIXME I had this error once:  "chr 0 not in genomeWithNewAdj"
            raise ValueError("chr %s not in genomeWithNewAdj = %s " % (c, [(c, len(chrom)) for (c, chrom) in genomeWithNewAdj.iteritems()]))
        (idxGl, idxGr) = myIntervals.findConsideredFlankingGenesIdxs(chrom, x, ignoredGeneNames=ignoredGeneNames)
        if idxGl is None or idxGr is None:
            # no new adjacency, this intergene corresponds to a chrom extremity
            return
        else:
            newAdj = OAdjacency(chrom[idxGl], chrom[idxGr])
            self.recordNewAdjacency(newAdj)

    def recordNewAdjacency(self, newAdj):
        #assert newAdj not in self.setRecoveredAncAdjacencies # not always the case since a recovery can be due to a gene deletion
        #(og1, og2) = newAdj
        # if myIntervals.geneExtremityFromGene(og1, +1) in self.setGeneExtrOfConservedAncBlocksWithoutAncChromExtr and\
        #     myIntervals.geneExtremityFromGene(og2, -1) in self.setGeneExtrOfConservedAncBlocksWithoutAncChromExtr:
        # This previous condition is not working since it also consider fusions of two previously different synteny blocks
        if newAdj in self.setOnceDisruptedAncAdjs:
            self.nbRecoveredAncAdjs += 1
            self.setRecoveredAncAdjacencies.add(newAdj)
            self._fuseConservedAncestralBlocks(newAdj)
        else:
            pass

    # This function is useless, user should directly use self.conservedAncestralBlocks, that will give him faster exactly
    # the same result
    @myTools.deprecated
    def computeSbs(self, ancestralGenome, removedAncGeneNames=set()):
        assert isinstance(ancestralGenome, LightGenome)
        # This set of broken gene extremities define synteny blocks as blocks of genes free of breakpoints
        # setGEaEoBFoB = 'setGeneExtrsAtExtremitiesOfBlocksFreeOfBreakoints'
        # This set of broken gene extremities define synteny blocks as blocks of genes that appear as free of breakpoints
        # (which are detectable)
        # sGEaEoCAB = setGeneExtrsAtExtremitiesOfConservedAncestralBlocks

        # The "conserved ancestral block" may have undergo a rearrangement followed by another that repaired it
        if not self.setBrokenGeneExtrs and not self.setOnceDisruptedAncAdjs:
            print >> sys.stderr, ("Warning: be sure that breakPointAnalyser has recorded breakpoints.",
                                  "\nIt does not contain any recorded breakpoint for the moment.")

        (ancestralGenome, _, (nbChrLoss, nbGeneLoss)) = myMapping.remapFilterGeneContent(ancestralGenome,
                                                                                         removedNames=removedAncGeneNames)
        assert nbGeneLoss == len(removedAncGeneNames), "%s %s %s" % (nbGeneLoss, len(removedAncGeneNames), nbChrLoss)

        # This next assertion is just to be sure, it is not necessary that ancestralGenome has no internal dict of
        # gene positions
        # function intergeneFromGeneExtremity needs a genome with an  dict g2p={..., geneName: genePosision, ...} up to date
        ancestralGenome.computeDictG2P()
        assert ancestralGenome.withDict is True

        brokenInnerIntergenesOnChr = collections.defaultdict(set)
        for brokenInnerGeneExtrs in self.setBrokenGeneExtrs - self.setAncChromExtremities:
            assert brokenInnerGeneExtrs[0] not in removedAncGeneNames
            (c, x) = myIntervals.intergeneFromGeneExtremity(brokenInnerGeneExtrs, ancestralGenome)
            brokenInnerIntergenesOnChr[c].add(x)

        self.setDirsuptedAncAdjacencies = self.setOnceDisruptedAncAdjs - self.setRecoveredAncAdjacencies
        self.setDetectableBrokenGeneExtrs = set()
        for adjacency in self.setDirsuptedAncAdjacencies:
            # we add all broken gene extremities corresponding to the disrupted adjacencies
            self.setDetectableBrokenGeneExtrs |= set(myIntervals.geneExtremitiesFromAdjacency(adjacency))
        # extremities of detectable contiguous ancestral regions, blocks of anc genes free of apparent breakpoints
        # detectableBreakpoints is not in self since detectableBreakpoints has a meanning only if it is with
        # ancestralGenome, which is not in self.
        detectableBrokenIntergenesOnChr = collections.defaultdict(set)
        for detectableBrokenGeneExtr in self.setDetectableBrokenGeneExtrs:
            # FIXME does it work well on chromosome extrimities
            (c, x) = myIntervals.intergeneFromGeneExtremity(detectableBrokenGeneExtr, ancestralGenome)
            detectableBrokenIntergenesOnChr[c].add(x)

        # Just for information:
        setDetectableBrokenIntergenes = set()
        for (c, setIntergenes) in detectableBrokenIntergenesOnChr.iteritems():
            for x in setIntergenes:
                setDetectableBrokenIntergenes.add((c, x))
        setBrokenIntergenes = set()
        for (c, setIntergenes) in brokenInnerIntergenesOnChr.iteritems():
            for x in setIntergenes:
                setBrokenIntergenes.add((c, x))

        # # set1 <= set 2 <=> set1.issubset(set2)
        # FIXME next assertion sent me AssertionError: ['16:0'] once... ['19:0'] another time...
        # FIXME once I had ['10:0 (len(c)=443)', '29:991 (len(c)=991)', '0:24 (len(c)=24)']
        assert setDetectableBrokenIntergenes <= setBrokenIntergenes, "%s" % \
                                                                     ["%s:%s (len(c)=%s)" % (c, x, len(ancestralGenome[c])) for (c, x) in setDetectableBrokenIntergenes - setBrokenIntergenes]
        if setDetectableBrokenIntergenes == setBrokenIntergenes:
            assert len(self.setRecoveredAncAdjacencies) == 0
            print >> sys.stderr, 'No recovered Ancestral Adjacency'
        else:
            setUndetectableBrokenGeneExtrs = setBrokenIntergenes - setDetectableBrokenIntergenes
            print >> sys.stderr, "There are %s undetectable broken intergenes: %s" % \
                                 (len(setUndetectableBrokenGeneExtrs), sorted(list(setUndetectableBrokenGeneExtrs)))

        self.detectableSyntenyBlocks = collections.defaultdict(list)
        # For each chromosome that did not undergo any rearrangement, the synteny block is equal to the old chromosome
        for c in set(ancestralGenome.keys()) - set(detectableBrokenIntergenesOnChr.keys()):
            self.detectableSyntenyBlocks[c] = [ancestralGenome[c]]
        for (c, xs) in detectableBrokenIntergenesOnChr.iteritems():
            assert len(xs) > 0
            for (x1, x2) in myTools.myIterator.slidingTuple(sorted(list(xs | {0, len(ancestralGenome[c])}))):
                assert x1 < x2
                self.detectableSyntenyBlocks[c].append(ancestralGenome[c][x1:x2])

        # DEBUG assertion, conservation of genes
        nbGenesInSbs = sum(len(sb) for (c, listOfSbs) in self.detectableSyntenyBlocks.iteritems() for sb in listOfSbs)
        nbGenesInAnc = sum(len(chrom) for chrom in ancestralGenome.values())
        assert nbGenesInSbs == nbGenesInAnc, "%s == %s" % (nbGenesInSbs, nbGenesInAnc)

        self.nbdetectableSyntenyBlocks = sum(len(listOfSbs) for (c, listOfSbs) in self.detectableSyntenyBlocks.iteritems())
        self.nbMonoGenicSbs = sum((1 if len(sb) == 1 else 0) for (c, listOfSbs) in self.detectableSyntenyBlocks.iteritems() for sb in listOfSbs)

    def printStats(self, stream=sys.stderr):
        print >> stream, "number of bra-standard inner breakpoints = %s" % self.nbStandardInnerBreakpoints
        print >> stream, "number of bra-inner breakpoints reused = %s" % self.nbInnerBreakpointsReuse
        if self.nbBreakpoints != 0:
            percentageBreakpointsReuse = "%.2f%%" % (100 * (float(self.nbBreakpointsReuse) / float(self.nbBreakpoints)))
        else:
            percentageBreakpointsReuse = "No breakpoint"
        print >> stream, "... number of breakpoints = %s, number of breakpoints reused = %s (%s)" % \
                             (self.nbBreakpoints, self.nbBreakpointsReuse, percentageBreakpointsReuse)
        print >> stream, "number of bra-breakpoints at chromosomal extremities = %s" % self.nbBreakpointsAtChromExtremity
        print >> stream, "number of bra-undetectable inner breakpoints (adjacency broken then recovered) = %s" % (self.nbRecoveredAncAdjs - self.nbBreakpointsWithinRecoveredAdjacencies)
        #print >> sys.stderr, "number of bra-recovered ancestral adjacencies = %s" % self.nbRecoveredAncAdjs
        # print >> sys.stderr, "number of bra-breakpoints within recovered adjacencies = %s" % self.nbBreakpointsWithinRecoveredAdjacencies
        print >> stream, "number of bra-synteny-blocks removed because of gene loss = %s" % self.nbAncChromsRemovedBecauseOfGeneLoss

        # Method 1
        #print >> sys.stderr, "# First method 1)"
        print >> stream, "number of bra-empirical detectable synteny blocks = %s" % len(self.conservedAncestralBlocks.keys())
        sbLens = sorted([len(sb) for (c, sb) in self.conservedAncestralBlocks.iteritems()])
        nbSbsOfLen = collections.Counter(sbLens)

        print >> stream, "distribution of bra-sb lengths (nbSbs:length) = %s" % myMaths.myStats.distribSummary(nbSbsOfLen)

        # debug
        # print >> stream, [sb for (c, sb) in self.conservedAncestralBlocks.iteritems() if len(sb) == 2]


        # # Method 2, returns exactly the same results
        # print >> sys.stderr, "# Second method 2)"
        # bra.computeSbs(iniGenome, removedAncGeneNames=setAncGenesLost)
        # print >> sys.stderr, "2) number of self-empirical detectable synteny blocks = %s" % self.nbdetectableSyntenyBlocks
        # sbLens = sorted([len(sb) for (c, listOfSbs) in self.detectableSyntenyBlocks.iteritems() for sb in listOfSbs])
        # nbSbsOfLen = collections.Counter(sbLens)
        # distribSbLens = [" %s:%s" % (nbSbsOfLen[sbLen], sbLen) for sbLen in sorted(nbSbsOfLen.keys())]
        # distribSbLens = distribSbLens[:5] + ["..."] + distribSbLens[-3:]
        # print >> sys.stderr, "2) distribution of bra-sb lengths (nbSbs:length) = %s" % ",".join(distribSbLens)
@myTools.verbose
@myTools.deprecated
def rewriteSbsWithOrthologs(sbs_A1_D1, sbs_A2_D2,
                  f_A1_Ds=None, f_A2_Ds=None,
                  f_LCA_Ds=None, verbose=False):
    """
    Combine two empirical sbs into one synthetic empirical sb

    For instance let us consider the evolution from Euarchontoglire to Homo.sapiens and Mus.musculus.
    MagSimus will return two empirical sbs:
        1) sbs.genes.Euarchontoglire.Homo.sapiens.list.bz2
        2) sbs.genes.Euarchontoglire.Mus.musculus.list.bz2
    combineTwoSbs(empiricalSb1, empiricalSb2) returns sbs.genes.Homo.sapiens.Mus.musculus.list.bz2
    This last file corresponds to the synteny blocks found by PhylDiag when Homo.sapiens is compared with
    Mus.musculus with the families of Euarchontoglire

    Warning, please ensure that LCA is LCA(A1, A2) and also LCA(D1, D2)

    :param sbs_A1_D1: LightGenome (with gene names in f_A1_Ds, if f_A1_Ds is not None)
    :param sbs_A2_D2: LightGenome (with gene names in f_A2_Ds, if f_A2_Ds is not None)
    :param f_A1_Ds families from A1 to Ds (Descendants)
    :param f_A2_Ds families from A2 to Ds (compatible Descendants)
    :param f_LCA_Ds families from the LCA of A1 and A2 to Descendants
    :return: LightGenome empiricalSb_A1_A2 with gene names of f_LCA_Ds
    """
    assert isinstance(sbs_A1_D1, LightGenome) and isinstance(sbs_A2_D2, LightGenome)
    if f_A1_Ds and f_A2_Ds and f_LCA_Ds:
        assert isinstance(f_A1_Ds, Families) and isinstance(f_A2_Ds, Families) and isinstance(f_LCA_Ds, Families)
    ###########################################################################################
    # 1) rewrite sbs_A1_D1, sbs_A2_D2 with family names of families_LCA_A1_A2
    ###########################################################################################
    if (f_A1_Ds is None and f_A2_Ds is None) or f_A1_Ds == f_A2_Ds:
        sbs_A1_D1_fLCA = copy.deepcopy(sbs_A1_D1)
        sbs_A2_D2_fLCA = copy.deepcopy(sbs_A2_D2)
        nbSbsLoss1 = 0
        nbSbsLoss2 = 0
        nbGeneLoss1 = 0
        nbGeneLoss2 = 0
    else:
        assert f_LCA_Ds is not None
        f_LCA_A1 = myLightGenomes.f_A0_A1_from_f_A0_Ds_and_f_A1_Ds(f_LCA_Ds, f_A1_Ds)
        assert all(len(sb) > 0 for sb in sbs_A1_D1.values())
        sbs_A1_D1 = myMapping.labelWithOrthologNames(sbs_A1_D1, f_LCA_A1)
        (sbs_A1_D1_fLCA, _, (nbSbsLoss1, nbGeneLoss1)) = myMapping.remapFilterGeneContent(sbs_A1_D1, {None})

        f_LCA_A2 = myLightGenomes.f_A0_A1_from_f_A0_Ds_and_f_A1_Ds(f_LCA_Ds, f_A2_Ds)
        assert all(len(sb) > 0 for sb in sbs_A2_D2.values())
        sbs_A2_D2 = myMapping.labelWithOrthologNames(sbs_A2_D2, f_LCA_A2)
        (sbs_A2_D2_fLCA, _, (nbSbsLoss2, nbGeneLoss2)) = myMapping.remapFilterGeneContent(sbs_A2_D2, {None})

    sg1 = sbs_A1_D1_fLCA.getGeneNames()
    sg2 = sbs_A2_D2_fLCA.getGeneNames()

    if f_LCA_Ds:
        familyNamesInLCA = set(fn for (fn, dns) in f_LCA_Ds)
        assert sg1 <= familyNamesInLCA, "%s" % (sg1 - familyNamesInLCA)
        assert sg2 <= familyNamesInLCA, "%s" % (sg2 - familyNamesInLCA)

    return (sbs_A1_D1, sbs_A2_D2)

def combineTwoSbs(sbs_1, sbs_2, verbose=False):
    """
    Combine two empirical sbs into one synthetic empirical sb

    :param sbs_1: sbs with positional orthologous gene names
    :param sbs_2: sbs with positional orthologous gene names
    :param verbose:
    :return: synthetic sbs
    """
    assert isinstance(sbs_1, myLightGenomes.LightGenome) and isinstance(sbs_2, myLightGenomes.LightGenome)
    sbs_1Gns = sbs_1.getGeneNames()
    sbs_2Gns = sbs_2.getGeneNames()

    ######################################
    # 1) reduce to the same set of genes #
    ######################################
    setGenesInBothSpecies = sbs_1Gns & sbs_2Gns
    setGenesToRemove = (sbs_1Gns | sbs_2Gns) - setGenesInBothSpecies
    (sbs_1, _, (nbRemovedSbsAnc, nbRemovedGenesAnc)) = myMapping.remapFilterGeneContent(sbs_1, setGenesToRemove)
    (sbs_2, _, (nbRemovedSbsDes, nbRemovedGenesDes)) = myMapping.remapFilterGeneContent(sbs_2, setGenesToRemove)
    assert sbs_1.getGeneNames() == sbs_2.getGeneNames()
    nbSbsRmed = nbRemovedSbsAnc + nbRemovedSbsDes

    ########################
    # 2) share breakpoints #
    ########################
    sbsExtremities1 = myIntervals.analyseGenomeIntoChromExtremities(sbs_1)
    sbsExtremities2 = myIntervals.analyseGenomeIntoChromExtremities(sbs_2)
    sbsExtremities = sbsExtremities1 | sbsExtremities2
    syntheticSbs = copy.deepcopy(sbs_1)
    assert all(len(sb) > 0 for sb in syntheticSbs.values())
    syntheticSbs.computeDictG2P()
    for breakpointGeneExtr in sbsExtremities:
        (c, x) = myIntervals.intergeneFromGeneExtremity(breakpointGeneExtr, syntheticSbs)
        syntheticSbs = myEvents.performFission(syntheticSbs, (c, x), updateGenomeDict=True)
    assert all(len(sb) > 0 for sb in syntheticSbs.values())
    assert isinstance(syntheticSbs, LightGenome)

    return (syntheticSbs, nbSbsRmed)

def integrateAllBreakpoints(evolutivePathFromLCAtoD, empiricalSbs):
    """
    :param evolutivePathFromLCAtoD: example ['Amniota', 'Theria', 'Boreoeutheria', 'Euarchontoglires', 'Homo sapiens']
    :param empiricalSbs: dict of sbs with keys (Amniota, Theria, ..., 'Homo sapiens', ...)
    :param ancGenes: dict of ancGenes Families with keys (Amniota, Theria, ..., ...)
    :return: sbs
    """
    # ancestor -> descendant chain
    a_d_chain = myTools.myIterator.slidingTuple(evolutivePathFromLCAtoD, width=2)
    (lca, d) = a_d_chain.next()
    sbs = empiricalSbs[(lca, d)]
    assert all(len(sb) > 0 for sb in sbs.values())
    old_d = d
    # a: ancestor
    # d: descendant
    nbSbsRmed = 0
    for (a, d) in a_d_chain:
        assert a == old_d
        new_sbs = empiricalSbs[(a, d)]
        assert all(len(sb) > 0 for sb in new_sbs.values())
        # print >> sys.stdout, "in integrateAllBreakpoints %s->%s" % (a, d)
        (sbs, tmp_nbSbsRmed) = combineTwoSbs(sbs, new_sbs, verbose=True)
        nbSbsRmed += tmp_nbSbsRmed
        assert all(len(sb) > 0 for sb in sbs.values())
        old_d = d
    return (sbs, nbSbsRmed)

def computePairwiseSyntenyBlocksOfMagSimus(speciesTree, pSimGenomes, pAncGenes, pReferenceSbs):
    # load families
    ancGenesOf = {}
    for s in speciesTree.listAncestr:
        ancGenesOf[s] = myLightGenomes.Families(pAncGenes % s)
    # load genomes
    genomeOf = {}
    for s in speciesTree.allNames:
        print >> sys.stderr, pSimGenomes % str(s)
        genomeOf[s] = myLightGenomes.LightGenome(pSimGenomes % str(s))
    # load sbs
    simSbsOf = {}

    def loadSbs(speciesTree, pSbs, parent):
        for (child, bLength) in speciesTree.items[parent]:
            simSbsOf[(parent, child)] = myLightGenomes.LightGenome(pSbs % (parent, str(child)))
            if child in speciesTree.items:
                # if child is an internal node of the species tree
                loadSbs(speciesTree, pSbs, child)
    loadSbs(speciesTree, pReferenceSbs, speciesTree.root)

    # for each pairwise comparison of two extant species compute the corresponding sbs
    for (sp1, sp2) in itertools.combinations(speciesTree.listSpecies, 2):
        lca = speciesTree.lastCommonAncestor([sp1, sp2])
        # print lca
        lca_genome = genomeOf[lca]
        nbGs_ini = len(lca_genome.getGeneNames(checkNoDuplicates=True))

        # sbs in the lineage from lca to sp1
        (sbs1, nbSbsRmed1) = integrateAllBreakpoints(speciesTree.dicLinks[lca][sp1], simSbsOf)
        # sbs in the lineage from lca to sp2
        (sbs2, nbSbsRmed2) = integrateAllBreakpoints(speciesTree.dicLinks[lca][sp2], simSbsOf)
        assert all(len(chrom) > 0 for chrom in sbs1.values())
        assert all(len(chrom) > 0 for chrom in sbs2.values())
        assert len(sbs1.getGeneNames(asA=list)) == len(sbs1.getGeneNames(asA=set))
        assert len(sbs2.getGeneNames(asA=list)) == len(sbs2.getGeneNames(asA=set))
        (sbs, nbSbsRmed3)  = combineTwoSbs(sbs1, sbs2, verbose=False)
        nbGsRmed = nbGs_ini - len(sbs.getGeneNames(checkNoDuplicates=True))
        nbSbsRmed = nbSbsRmed1 + nbSbsRmed2 + nbSbsRmed3
        # to have the species combination in alphabetical order
        speciesNames = sorted([sp1, sp2])
        print >> sys.stderr, "Computation of %s: %s sbs removed, %s genes removed" % \
                             (pReferenceSbs % (speciesNames[0], speciesNames[1]), nbSbsRmed, nbGsRmed)
        # sort sbs by decreasing sizes
        sbs.sort()
        with myFile.openFile(pReferenceSbs % (speciesNames[0], speciesNames[1]), 'w') as f:
            sbs.printIn(f)
