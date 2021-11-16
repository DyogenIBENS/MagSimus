# -*- coding: utf-8 -*-

import sys
import collections
import random

import utils.myTools as myTools
import utils.myFile as myFile

import myEvolutionProbas


class Clusterator:

    def __init__(self, geneTimeline, families=None):
        allvect = {}
        byname = {}
        byid = {}
        print >> sys.stderr,\
            "Loading vectors of gene evolution (genetimeline) from %s " % geneTimeline
        f = myFile.openFile(geneTimeline, 'r')
        for l in f:
            # FIXME : For the moment all 'None' values are replaced by 0 values. A
            # better consideration of 'None's in the coevolution score computation
            # could be done.
            t = l.replace('None', '0').split("\t")
            gname = t[0]
            ev = []
            for x in t[1:]:
                ev.extend(eval(x))
            # ev = (0,0,1,0,0,1,0,....), is the events vector, a concatenation of the
            # duplications vector with the deletions vector. It has 2*nbBranches
            # components.
            ev = tuple(ev)
            # len(ev) # returns 2*8=16

            # Collapse identical genes to save memory
            # Simultaneous build of:
            # - 'allvect': {ev1 : 1, ev2: 2, ...}
            # - 'byname': {g1name : index1, g2name : index2, ...}
            #  - 'byid':
            # Since there are few possible vectors, working with a list of all
            # possible events vectors saves more memory.
            # A given events vector <=> an index in 'allvect'
            # Now a surjection from the set of genes to the set of indices of events
            # vectors in 'allvect'
            # Instead of recording events vectors for each gene, the corresponding
            # index in 'allvect' is recorded.
            # If 2 genes have the same events vector they have the same index in
            # allvect.
            # For instance:
            # Given byname = {g0name: 1, g1name: 1} and allvect = {v0: 0, v1: 1}
            # Then g0name and g1name have the same events vector v1.

            # 'allvect' is a dict rather than a list because searches are faster
            # using dicts than using lists.
            # Either 'allvect' already contains 'ev' and thus its value is already
            # in 'allvect',
            # Or, 'allvect' does not contain 'ev' and thus we add 'len(allvect)' as
            # a value of 'ev'
            # Rq: setdefault(key[, default]) : if key is in the dictionary, return
            # its value. If not, insert key with a value of default and return
            # default.
            if families is None:
                byname[gname] = allvect.setdefault(ev, len(allvect))
            else:
                byname[families.getFamID(gname)] = allvect.setdefault(ev, len(allvect))
            # byid is a dict {idx in allvect -> corresponding events vector}
            byid[byname[gname]] = ev
        f.close()
        self.allvect = allvect
        self.byname = byname
        self.byid = byid

    # FIXME : (v1, v2, lstindices, byid) is not cacheable:
    # > [coevoluTionScore: 0 values cached, 1 calls, 0.00x speedup]/<function coevoluTionScore at 0x965a28>)
    # FIXME consider 'None's and values >1 in the events vectors
    # 'lstindices' is the list of the branches (the indices of the events
    # vectors 'v1' and 'v2') that are considered.
    # For the initial genome all branches are considered, then only the indices of
    # the subTrees below the current ancestor are considered.
    @myTools.memoizeMethod
    def coevoluTionScore(self, v1, v2, lstindices):
        # For optimisation, 'v1' are 'v2' the indices of the corresponding event
        # vectors in 'allvect'
        v1 = self.byid[v1]
        v2 = self.byid[v2]
        # v1 and v2 are now event vector e.g. (0,0,1,0,0,1, ....) with
        # 2*(nbBranches) elements
        s = 0  # Number of coevolutions
        n = 0  # Number of singular evolutions
        # s = nb of times that v1[i] == v2[i]
        # n = nb of times that:
        #       either: v1[i] != 0
        #       or:     v2[i] != 0
        #       (this include (v1[i] and v2[i] are != 0))
        for i in lstindices:
            # Warning, for optimisation purposes we do not test (v1[i]>0)
            if v1[i]:  # v1[i] == 1
                n += 1
                s += v2[i]
            else:
                n += v2[i]
        # Nb of times that both genes undergo the same event over the total number
        # of events undergone by one gene or the other
        coevolScore = s/float(n) if n != 0 else 0
        # coevoluTionScore is always >= 0
        return (coevolScore, n, s)

    # Find the best block to insert at the location 'x' of the chromosome 'chrom'
    # Version where the inserted block is chosen among a list of blocks
    # - A position (chrom, x) of the genome is fixed
    # - 'lblocks' is the list of blocks that are candidate for the insertion. Each
    #       block is a list of genes
    # - The block that maximises the sum of the coevolutionScores with at its both
    #       extremities is chosen
    #   Each coevolutionScore calculation is performed between a gene at the
    #       extremity of the block and a gene on the chromosome at the location of
    #       the insertion
    # - Every block is involved twice in the calculation, once with a positive
    #       orientation and once turned over, in a negative orientation
    # Look at the coevolutionScores of the genes at the extremities of a candidate
    # block with the genes around the location 'x' (at the location x-1 and x+1)
    # This function is also used to know if it is a good thing to substitute a gene
    # for a block of its duplicates. This would correspond to a gene tandem
    # duplications.
    # 'B' stands for blocks, here the insertion position is fixed.
    # Returns the list of blocks that have the highest coevolScores for an insertion
    # at the location 'x' of 'chrom', with their best orientation or a random
    # orientation when it is impossible to determine the best orientation.
    def scoreWithNeighboursB(self, lblocks, chrom, x, lstindices, byname):
        # 'name left': events vector of the left gene
        nl = byname[chrom[max(x-1, 0)][0]]
        # 'name right': events vector of the right gene
        # FIXME : J. Lucas added +1
        nr = byname[chrom[min(x+1, len(chrom)-1)][0]]

        bestS = 0
        # list of the chosen block (or blocks if equality)
        bestB = []
        for block in lblocks:
            # 1st gene of the block, 'l': left
            bl = byname[block[0]]
            # last gene of the block, 'r': right
            br = byname[block[-1]]
            # score of the insertion in the 'f': forward orientation
            sf = self.coevoluTionScore(nl, bl, lstindices)[0] +\
                self.coevoluTionScore(nr, br, lstindices)[0]
            # score of the insertion in the 'r': reverse orientation
            sr = self.coevoluTionScore(nl, br, lstindices)[0] +\
                self.coevoluTionScore(nr, bl, lstindices)[0]
            if sf == sr:
                (score, strand) = (sf, myEvolutionProbas.randomStrand())
            else:
                if sf > sr:
                    (score, strand) = (sf, 1)
                else:
                    (score, strand) = (sr, -1)
            if score > bestS:
                bestS = score
                bestB = [(strand, block)]
            elif score == bestS:
                bestB.append((strand, block))
        # bestB contains at least one block if lblocks is not empty
        return bestB

    # Find the best place to insert the 'block' among the locations contained in
    # 'lpos'.
    # - A bock of genes is fixed.
    # It is tried to bie inserted at all the locations in all orientations. The
    # location and orientation with the highest score of coevolution is returned.
    # If many positions-orientations are equivalent a list of these positions are
    # returned.
    # 'lpos' has the form of : {chrom : [xpos1, xpos2, ...], ...}
    # 'P' stands for position
    #
    # Warning, to take care of extremities, bgenome has the form :
    # (c.f. insertAtTheBestPlace)
    # for (c,chrom) in enumerate(genome):
    #     bgenome[c] = [byname[x[0]] for x in chrom]
    #     bgenome[c].append(bgenome[c][-1])
    #     bgenome[c].append(bgenome[c][0])
    # If bchrom = genome[c]:
    # bchrom[i] = chrom[i] for i in [0,len(chrom)-1] and
    # bchrom[-1] = chrom[0] = bchrom[len(chrom)+1]
    # bchrom[-2] = chrom[-1] = bchrom[len(chrom)]
    # and bchrom[-3] = chrom[-1] = bchrom[len(chrom)-1]
    # soit bchrom=[chrom + [chrom[-1]] + [chrom[0]] ]
    def scoreWithNeighboursP(self, block, genome, allInsertPos, indices, clusterator):
        # evolution-vector-genome
        eVgenome = {}
        # as allInsertPos, allInsertPosInChrom contains all the positions where it is possible to
        # insert a block, but only in chrom this time
        # allInsertPosInChrom is only used in this if condition
        allInsertPosInChrom = collections.defaultdict(list)
        for (chrom,x) in allInsertPos:
            allInsertPosInChrom[chrom].append(x)
        for (c, chrom) in genome.iteritems():
            # assert len(chrom) > 0 # DEBUG assertion
            eVgenome[c] = [clusterator.byname[x[0]] for x in chrom]
            # Here there is a trick to consider extremities,
            # c.f. 'scoreWithNeighboursP'
            eVgenome[c].append(eVgenome[c][-1])  # (cf 'scoreWithNeighboursP')
            eVgenome[c].append(eVgenome[c][0])   # (cf 'scoreWithNeighboursP')
            # FIXME : Joseph added +1. This corresponds to a possible insertion
            # at the end of the chromosome.
            # DEBUG assertion
            #assert len(allInsertPosInChrom[c]) == len([(ch,tata) for (ch,tata) in allInsertPos if ch == c])

        bl = clusterator.byname[block[0][0]]
        br = clusterator.byname[block[-1][0]]
        bestS = 0
        bestInsertions = []
        for c in allInsertPosInChrom:
            bchrom = eVgenome[c]
            assert len(allInsertPosInChrom[c]) <= len(eVgenome[c])+1
            for x in allInsertPosInChrom[c]:
                # x is in [0,len(lpos[c])-1] = [0,len(chrom)]
                # assert x < len(bchrom) # DEBUG
                nl = bchrom[x - 1]
                nr = bchrom[x]
                # if x = 0: bchrom[x-1] = bchrom[0] by construction
                # and nl = nr
                # if x = len(chrom): bchrom[x] = bchrom[x-1] by construction
                # and nl = nr
                sf = self.coevoluTionScore(nl, bl, indices)[0] +\
                    self.coevoluTionScore(nr, br, indices)[0]
                sr = self.coevoluTionScore(nl, br, indices)[0] +\
                    self.coevoluTionScore(nr, bl, indices)[0]
                if sf == sr:
                    (score, strand) = (sf, myEvolutionProbas.randomStrand())
                else:
                    if sf > sr:
                        (score, strand) = (sf, 1)
                    else:
                        (score, strand) = (sr, -1)
                if score > bestS:
                    bestS = score
                    bestInsertions = [(strand, c, x)]
                elif score == bestS:
                    bestInsertions.append((strand, c, x))
        return bestInsertions

    # Segment a list of genes and clusterise them thanks to the coevolutionScore
    # Returns allblocks, format : [Bloc1, Bloc2, ...] avec Bloc = [geneName1, geneName2, ...]
    @myTools.verbose
    def mkClusterByScore(self, lstGenes, indices, forceClustering, clusteredGenesRatio,
                         prefix='', verbose=False):
        # lstGenes : list of the gene names
        # indicies : frozenset of the indices of eventsVectors t consider
        #   usually only the branches below the current node are considered
        print prefix, "clustering",

        # if len(indices) == 0 <=> means that we are on a leaf
        if (len(indices) == 0) and (not forceClustering):
            print "singletons", len(lstGenes), "|".join(lstGenes)
            # return lists with only one gene, no clustering
            return [[g] for g in lstGenes]

        # ng: number of genes
        ng = len(lstGenes)
        if ng == 0:
            print "empty", len(lstGenes)
            return []

        # assert len(lstGenes) > 0 # DEBUG
        print "ok", len(lstGenes),

        #FIXME warning after that !!

        # list of all blocks,
        # allblocks format : [Bloc1, Bloc2, ...] avec Bloc = [geneName1, geneName2, ...]
        allblocks = []
        # boolean True when it is requested to build a new block.
        # It happens when the preceding block cannot be extended any more or if the
        # random consideration of the clusterRatio yields a new block.
        # If forceclustering = True, newBlock is True for the first gene and false
        # for all the others and THERE IS ONLY ONE BLOCK!
        newBlock = True
        # 'i' is a random integer in [0, ..., ng-1]
        i = random.choice(range(ng))
        # Deep copy because we will modify it
        lstGenes2 = list(lstGenes)  # FIXME
        reasons = []
        # coevolution score
        coevS = None

        nbGenes = len(lstGenes2)
        progressBar = myTools.ProgressBar(nbGenes, step=1)
        while True:
            progressBar.printProgressIn(sys.stderr, nbGenes - ng)
            if newBlock:
                allblocks.append([])
            # refG : reference gene
            # refG is extracted from list 'lstGenes2' at the index 'i'
            refG = lstGenes2.pop(i)
            ng -= 1
            #FIXME : for optimisation purposes it is not requested to return 'm'
            # refGene is inserted in the last/current block, format:
            # (coevolutionScore, gene)
            allblocks[-1].append((coevS, refG))
            if ng == 0:
                break

            # Start a new block with a complete random gene
            # It happens more when clusteredGenesRatio is low
            if random.random() > clusteredGenesRatio:
                coevS = None
                # i, index of a random gene
                i = random.choice(range(ng))
                newBlock = True
                # TODO R as 'r'andom ?
                reasons.append("R")

            # Genes with best coevS are added to the current block
            else:
                # Here we continue the clustering in the current block
                if len(indices) > 0:
                    # refV is the index of the events vector of the gene 'refG' in
                    # 'allvect'
                    refV = self.byname[refG]
                    # 'l': list of all coevoluTion scores for refG with all the genes
                    # of 'lstGenes2'.
                    # 'l' format: [(coevS, indexOfGenesInlstGenes2), ...]
                    l = []
                    for (i, g) in enumerate(lstGenes2):
                        # coevolution score between gene 'refG' and 'g'
                        (coevolScore, n, s) = self.coevoluTionScore(refV, self.byname[g], indices)
                        l.append((coevolScore, n, s, i))
                else:
                    # Case when len(indices) == 0, True for extant species
                    # If it is True, it also means that 'forceClustering' = True
                    # (c.f. 'if' at the beginning of the function)
                    # TODO understand: score = None, means no gene to clusterise
                    l = [(None,)]
                # Search max(l) owing to the first component of each element,
                # i.e. coevoluTionScore(refG, g) and take the first component,i.e.
                # the best score. It is a float.
                coevS = max(l, key=lambda x: x[0])[0]
                # If all genes have a coevoluTionScore = 0 elif len(indices) == 0
                # (extant species)
                if coevS == 0 or coevS is None:
                    # A random gene is chosen
                    i = random.choice(range(ng))
                    # If clustering is forced, all genes will be clusterised, else a
                    # new block is started.
                    newBlock = not forceClustering
                    # TODO 'l' and 'L', meaning of this letter ?
                    reasons.append("l" if forceClustering else "L")
                else:
                    newBlock = False
                    bestCoevolGenes =\
                        [(_n, _s, _i) for (_cS, _n, _s, _i) in l if _cS == coevS]
                    # Random choice among all genes that have a coevoluTionScore = m
                    (n, s, i) = random.choice(bestCoevolGenes)
                    #print >> sys.stderr, "B%.1f/%d/%d" % (coevS,n,s)
                    # TODO : B as 'b'locks ?
                    reasons.append("B%.1f/%d/%d" % (coevS, n, s))

        assert len(lstGenes2) == 0
        print "|".join(str(len(x)) for x in allblocks), "-".join(reasons),\
            "|".join("+".join([g for (_, g) in x]) for x in allblocks)

        if forceClustering and clusteredGenesRatio == 1:
            assert len(allblocks) == 1

        return allblocks
