# -*- coding: utf-8 -*-

#import myClustering
import random
import collections
import sys
import copy
import itertools
import operator

import myEvolutionProbas
import myChromosomalEvents
import myMagSimusTools

import utils.myLightGenomes as myLightGenomes
from utils.myLightGenomes import OGene as Gene
import utils.myTools as myTools


# Read the geneEvent file and return a list of genes that are 'renamed' and a
# list of genes that are 'gained'.
# 'renamed' genes include genes that are kept in the child species and also
# gene duplicates. 'gained' genes are the de novo gene births.
# Remember that in MagSimus2 the evolution of the gene set matches the real
# data.
def fixGeneEvents(genome, child, dist, speciesTree, geneEventsReader,
                  prefix=''):
    # Rq: in geneEvents
    # The 1st line only contains a set(['name1','name2',...]) corresponding to
    # the set of genes in the first ancestral genome.
    # Other lines contains:
    # {'n':['nt1, 'nt2']} set(['nd1', 'nd2']) set(['ng1', 'ng2'])
    # with t:translate, d:delete, g:gained
    # translated = name changes + duplications.
    # Excepted the 1st line, each line corresponds to a branch of the species
    # tree
    l = geneEventsReader.readline().split("\t")

    # The eval function changes a str into a python expression, it does as if
    # the str value was written without the " around, in the code
    # Here l[0] = "{'n':['nt1','nt2',...],...}"
    #                   => eval(l[0]) = {'n':['nt1','nt2',...],...}
    rename = eval(l[0])
    deleted = eval(l[1])
    gained = eval(l[2])
    assert deleted.isdisjoint(rename)
    assert deleted.isdisjoint(gained)
    assert gained.isdisjoint(rename.keys())

    # Add identical genes
    for chrom in genome.values():
        for (g, _) in chrom:
            if (g not in deleted) and (g not in rename):
                # These are genes that won't undergo any event, not even a
                # duplication, and that won't change there name neither
                rename[g] = [g]

    # FIXME
    # Losses of genes in 2X species are not well taken into account, 1/3rd of
    # the genome is removed later.
    # It is not because certain descendant genes are not seen in the child
    # genome (2X) that they are not here in real data. Thus all genes that are
    # supposed to be removed are kept.
    if child in speciesTree.lstEsp2X:
        # The 3 first letters of the 1st part of the child species, then the 3
        # first letters of the 2nd part of the name of the child species, these
        # last 3 letters are written upcases.
        # In other words, instead of removing genes on a branch that goes to a
        # 2X species, these genes are kept and are re-named SIMUG*****.
        prefix = (child[:3] + child.split()[1][:3]).upper()
        for (i, g) in enumerate(deleted):
            rename.add(OGene(g, ["SIMUG%s%06d" % (prefix, i)]))
        # these genes are removed from deleted
        deleted.clear()
        # FIXME rename.update( (g,["SIMUG%s%06d" % (prefix,i)]) for (i,(g,l)) in enumerate(trans.iteritems()) if len(l) == 0)

    return (rename, list(gained))


def renameAndTandemBlocks(genome, rename, indices, clusterator,
                          forceClustering, clusteredGenesRatio, tandemDupRatio, prefix=''):
    # nb of duplications
    # FIXME Will not be exactly the nb of tandem duplications since all gene events
    # are followed by all chromosomal rearrangements, thus breakpoints may
    # disrupt tandem blocks.
    nbGeneEvents = collections.defaultdict(int)
    # list of all blocks of duplicates that are not duplicated in tandem
    # format: set([Block1, Block2, ...]) with Block1 = ((g1,s1), (g2,s2), ...)
    dispersedTandemBlocks = set()
    for (c, chrom) in genome.iteritems():
        lastWasDeletion = False
        # 'sIdxDel' conserves the index of a 1st block of deleted genes
        sIdxDel = None
        nbDelOnThisChrom = 0
        for (idx, (g, _)) in enumerate(chrom):
            if g not in rename:
                # Consequence it means that it is a deleted gene
                print prefix,\
                    "geneevent genedeletion", g, "at %s:%d" % (c, idx+1)
                if not lastWasDeletion:
                    # if the previous gene was not a deletion, the index of the
                    # gene is recorded in 'sIdxDel'
                    sIdxDel = idx
                lastWasDeletion = True
                nbGeneEvents['geneDel'] += 1
                nbDelOnThisChrom += 1
                if nbDelOnThisChrom == len(chrom):
                    print >> sys.stderr, "Warning: gene deletions lead to an empty chromosome"
            # 'g' is in rename, so it is not deleted, it is kept.
            else:
                if lastWasDeletion:
                    # getIntervStr -> gene names around the intergene : g1-g2
                    print prefix, "geneevent segmentdeletion from %s:%d=(%s) to %s:%d=(%s)"\
                        % (c, sIdxDel, genome.getIntervStr(c, sIdxDel),
                           c, idx, genome.getIntervStr(c, idx))
                lastWasDeletion = False
        # Consider the case where the last gene is removed
        if lastWasDeletion:
            print prefix, "geneevent segmentdeletion from %s:%d=(%s) to %s:%d=(%s)" %\
                (c, sIdxDel, genome.getIntervStr(c, sIdxDel),
                 c, idx+1, genome.getIntervStr(c, idx+1))  # correction by Joseph
        # defines 'tmpChrom' in order to not disrupt 'chrom', because we are in
        # a for loop over the elements of 'chrom'.
        # WARNING : here tmpChrom may be null !
        tmpChrom = [OGene(g, strand) for (g, strand) in chrom if g in rename]
        # newChrom = buffer chromosome. It will replace 'chrom' at the end of the
        # loop
        newChrom = []
        # 'r' as rename
        for (idx, (g, strand)) in enumerate(tmpChrom):
            r = rename[g]
            if len(r) == 1:
                if r[0] != g:
                    print prefix, "geneevent rename", g, "to", r[0]
                newChrom.append(OGene(r[0], strand))
            elif len(r) >= 2:
                nbGeneEvents['geneDup'] += len(r)-1
                # 'g' has several child genes (duplicas) in its child genome.
                # They are contained in 'r'. 'r' also contains the new name of
                # the initial gene.
                # Identification of gene to be duplicated and distribution into
                # blocks
                # Consider 'r' as a list of genes and clusterise this list of
                # genes depending on the indices
                tandemBlocks = clusterator.mkClusterByScore(r, indices,
                                                            forceClustering,
                                                            tandemDupRatio,
                                                            prefix=prefix,
                                                            verbose=False)
                # tandemBlocks format: [(strand, block), ...]
                tandemBlocks = [[g for (m, g) in tb] for tb in tandemBlocks]
                # DEBUG assertion
                foo = len(tandemBlocks)
                if len(indices) > 0:
                    # 'lBestTbsForG' is the list of the best blocks(s) of 'tandemBlocks' to
                    # insert at the location of 'g'
                    lBestTbsForG = clusterator.scoreWithNeighboursB(tandemBlocks, tmpChrom, idx, indices, clusterator.byname)
                    if len(lBestTbsForG) > 0:
                        # list of clustered genes, considered segmental tandem
                        # duplicates of 'g'
                        (bestOrder, bestTbForG) = random.choice(lBestTbsForG)
                        tandemBlocks.remove(bestTbForG)  # FIXME Joseph correction
                        if bestOrder == -1:
                            bestTbForG.reverse()
                    else:
                        bestTbForG = random.choice(tandemBlocks)
                        tandemBlocks.remove(bestTbForG)  # FIXME Joseph correction
                else:
                    bestTbForG = random.choice(tandemBlocks)
                    tandemBlocks.remove(bestTbForG)  # FIXME Joseph correction
                # DEBUG assertion
                assert len(tandemBlocks) == foo - 1

                # 'g' is replaced by a list segmental tandem duplicates: bestTbForG.
                # Other duplicates (listed in 'tandemBlocks' for blocks != of bestTbForG[1])
                # are dispersed in the genome.
                print prefix, "geneevent duplication", g, " renamed to",\
                    "+".join(bestTbForG), ", new genes:",\
                    "|".join("+".join(x) for x in tandemBlocks if x != bestTbForG)
                # FIXME : same orientation for all genes duplicated in tandem or
                # not ? what is the biological reality ?
                # Previous 'g' is replaced by the best block of the duplicated
                # genes. That is similar to a segment of genes duplicated in
                # tandem.
                idxConservedStrand = random.choice(range(len(bestTbForG)))
                for (idx, g) in enumerate(bestTbForG):  # FIXME correction Joseph
                    # At least one gene with the initial orientation
                    if idx == idxConservedStrand:
                        newChrom.append(OGene(g, strand))
                    else:
                        newChrom.append(OGene(g, myEvolutionProbas.randomStrand()))
                for tb in tandemBlocks:
                    #assert tb != bestTbForG[1] # DEBUG assertion
                    # Rq: it is necessary tu put tuple otherwise it is a generator
                    tb = tuple(OGene(g, myEvolutionProbas.randomStrand()) for g in tb)
                    dispersedTandemBlocks.add(tb)  # FIXME : changed by Joseph
                    #dispersedTbs.update( tuple((g,randomStrand()) for g in b) for b in tandemBlocks if b != bestTbForG[1] ) # ici, creation des blocs de duplica a inserer un peu partout dans le genome # FIXME remove ?
            else:
                # 'r' cannot be empty
                raise

        # Replacement of the initial genome by the new one after evolution
        if len(newChrom) > 0:
            # newChrom overwrites chrom
            genome[c] = newChrom
        else:
            genome[c] = []
            #del genome[c] #WARNING, this could perturbate a loop on genome

    # Empty chromosomes are removed
    listOfCs = genome.keys()
    for c in listOfCs:
        if len(genome[c]) == 0:
            print >> sys.stderr, "Warning: Because of gene losses, one empty chromosome has been removed. This may lead to an unexpected nb of chromosomes in the child species."
            del genome[c]

    return (genome, nbGeneEvents, dispersedTandemBlocks)


# Insert some genes at the best place and keep the remaining genes in a list
def insertAtTheBestPlace(blocksToBeInserted, genome, indices, allInsertPos,
                         clusteredGenesRatio, clusterator, prefix=''):
    todo = []

    remainingBlocksToBeInserted = copy.deepcopy(blocksToBeInserted)

    # A part of the blocks are inserted at the best places
    if len(indices) > 0:

        # In 'newBlock' only a part of the blocks are conserved. Removed blocks
        # won't be clusterised.
        # 'nextnb' is the number of blocks that will be clusterised.

        nbRemainingExpectedInsertions = int(len(blocksToBeInserted) * clusteredGenesRatio)
        # 0 <= nextnb <= len(newBlocks) because 0 <= clusteredGenesRatio <= 1

        # This loop is the longest processing of the simulation:
        # For each cluster, all the insertions of remaining genes are tried.
        print >> sys.stderr, ("For every cluster of new genes we find its best "
                              "location thanks to the coevolution scores of "
                              "extremities with genes around the insertion position")
        iniRemainingExpectedInsertions = nbRemainingExpectedInsertions
        progressBar = myTools.ProgressBar(iniRemainingExpectedInsertions, step=1)
        while nbRemainingExpectedInsertions > 0:
            # Force Python's print function to output to the screen.
            # "to flush" = nettoyer a grandes eaux
            # Force python to write its buffer into sys.stdout in order to avoid
            # too long buffer and free memory.
            sys.stdout.flush()
            # Return a 'nextnb' length list of unique elements chosen from
            # the list 'blocksToBeInserted', the returned list is in selection
            # order.
            blocksThatWillBeInserted = random.sample(blocksToBeInserted, nbRemainingExpectedInsertions)
            # Remove from 'newBlocks' all the elements of 'blocksThatWillBeInserted'
            remainingBlocksToBeInserted.difference_update(blocksThatWillBeInserted)
            for b in blocksThatWillBeInserted:
                bestInsertions = clusterator.scoreWithNeighboursP(b, genome,
                                                                  allInsertPos,
                                                                  indices,
                                                                  clusterator)
                if len(bestInsertions) == 0:
                    # If the current block could not be inserted, it will be
                    # placed at a random position later.
                    remainingBlocksToBeInserted.add(b)
                    continue
                # A location for the insertion of the block is chosen
                # 's' = strand, 'c' = chrom, 'x' = pos
                (s, c, x) = random.choice(bestInsertions)
                # DEBUG assertions
                #assert x <= max([tata for (ch,tata) in allInsertPos if ch == c] ), "%s > max([ tata for (ch,tata) in allInsertPos if ch == %s] ) = %s" % (x, c, max([ tata for (ch,tata) in allInsertPos if ch == c] ) )

                print prefix,\
                    "cluster insertion optimized", "%s->%s in %s:%d=(%s)" %\
                    (b[0][0], b[-1][0], c, x, genome.getIntervStr(c, x)),
                # 'applyStrrand(b,s)' returns 'b' if 's' == +1 and 'b' in the
                # opposite direction (including genes' orientations) if 's' == -1
                todo.append((myChromosomalEvents.applyStrand(b, s), (c, x)))
                print  # FIXME ?
                # FIXME : why do we remove the position for forthcoming insertions ?
                allInsertPos.remove((c, x))
                nbRemainingExpectedInsertions -= 1
                progressBar.printProgressIn(sys.stderr, (iniRemainingExpectedInsertions - nbRemainingExpectedInsertions))
            # Used in the case where it was not possible to insert some blocks,
            # If ever there are some blocks that cannot be inserted we may
            # exhaust the stock of available blocks.
            # In this case we stop.
            nbRemainingExpectedInsertions = min(nbRemainingExpectedInsertions, len(remainingBlocksToBeInserted))
    return (genome, remainingBlocksToBeInserted, allInsertPos, todo)


# Insert genes in 'newGenes' and in 'dispersedDuplicates' in genome.
def insertNewGenes(newGenes, dispersedTandemBlocks, genome, indices,
                   forceClustering, clusteredGenesRatio, clusterator, nbGeneEvents, prefix=''):

    # FIXME just a name change or need a deep copy ?
    # blocksToBeInserted = copy.deepcopy(dispersedDuplicates)
    blocksToBeInserted = dispersedTandemBlocks
    for g in newGenes:
        print prefix, "geneevent newgene", g
        nbGeneEvents['geneGain'] += 1

    # Blocks of new genes
    # New genes are clusterised and they receive a random orientation.
    blocksOfNewGenes =\
        clusterator.mkClusterByScore(newGenes, indices,
                                     forceClustering,
                                     clusteredGenesRatio,
                                     prefix=prefix,
                                     verbose=False)
    blocksOfNewGenes = [[g for (m, g) in block] for block in blocksOfNewGenes]
    for b in blocksOfNewGenes:
        # changed by Joseph.
        # It is necessary tu put a tuple otherwise it is a generator.
        tmpBlock = tuple(OGene(g, myEvolutionProbas.randomStrand()) for g in b)
        # FIXME changed by Joseph
        blocksToBeInserted.add(tmpBlock)
        # 'blocksToBeInserted' format: set([Block1, Block2, ...])
        # with Block = [(g1,s1), g2,s2), ...]

    # DEBUG assertion
    #assert sum(len(block) for block in blocksToBeInserted) == len(newGenes) + sum(len(block) for block in dispersedDuplicates), "%s == %s + %s" % ( sum(len(block) for block in blocksToBeInserted) , len(newGenes) , sum(len(block) for block in dispersedDuplicates) )

    # 'allInsertPos' contains all the positions where it is possible to insert a block.
    allInsertPos = set()
    for (c, chromosome) in genome.iteritems():
        # assert len(chromosome) > 0 # DEBUG assertion
        # Joseph added +1 that corresponds to a possible insertion at the end of
        # the chromosome
        for x in xrange(len(chromosome)+1):
            allInsertPos.add((c, x))
    # set of possible pairs (c = chrom, x = position)

    # list of insertions to perform
    # 'todo' format: [((block),(chrom,x)), ...]
    # with block = set([(g1,s1),(g2,s2),...])
    (genome, remainingBlocksToBeInserted, allInsertPos, todo) =\
        insertAtTheBestPlace(blocksToBeInserted, genome, indices, allInsertPos,
                             clusteredGenesRatio, clusterator, prefix=prefix)

    # DEBUG assertion
    #assert sum(len(block) for block in remainingBlocksToBeInserted) + sum(len(block) for (block,_) in todo) == sum(len(block) for block in blocksToBeInserted), "%s +%s == %s" % (sum(len(block) for block in remainingBlocksToBeInserted) , sum(len(block) for (block,_) in todo) , sum(len(block) for block in blocksToBeInserted))

    # Remaining blocks are inserted at random locations
    # len(randomPosBlocks) locations are sampled in the genome to insert the
    # blocks in 'randomPosBlocks'.
    # 'itertools.izip' : makes an iterator that aggregates elements from each of
    # the iterables. Like zip() except that it returns an iterator instead of a
    # list.
    foo = itertools.izip(remainingBlocksToBeInserted,
                         random.sample(allInsertPos, len(remainingBlocksToBeInserted)))
    for (b, (c, x)) in foo:
        print prefix, "cluster insertion random", "%s->%s in %s:%d=(%s)" %\
            (b[0][0], b[-1][0], c, x, genome.getIntervStr(c, x)),
        strandedBlock = myChromosomalEvents.applyStrand(b, myEvolutionProbas.randomStrand())
        todo.append((strandedBlock, (c, x)))
        allInsertPos.remove((c, x))
        print

    # DEBUG assertion
    #assert sum(len(block) for (block,_) in todo) == sum(len(block) for block in blocksToBeInserted) , "%s == %s " % (sum(len(block) for (block,_) in todo), sum(len(block) for block in blocksToBeInserted))

    # Insertion of new blocks
    # Sort the list 'todo' depending of (c,x) by decreasing order.
    # This allows us to insert blocks at the end of chromosome without
    # disrupting the index locations that are at the beginning of the chromosomes.
    todo.sort(reverse=True, key=operator.itemgetter(1))
    for (l, (c, x)) in todo:
        # 'l' is inserted at the index 'x' on chromosome 'c'
        genome[c][x:x] = l
        #Rq l[x:x]='o' insert 'o' at the position :
        #               x                       if 0 <= x < len(l)
        #               len(l)+x                if -len(l) <= x < 0
        #               at the end of l         if x >= len(l)
        #               at the beginning of l   if x ==0 or if x < -len(l)

    return (genome, nbGeneEvents)


# Only genic events
def performGeneEvents(rename, newGenes, genome, indices, forceClustering,
                      clusteredGenesRatio, tandemDupRatio, clusterator, printIntermGenomes=False, prefix=''):
    # nb of genes in the genome (Number INItial)
    nbIniGenes = sum(len(chrom) for chrom in genome.values())

    if printIntermGenomes:
        famId = myMagSimusTools.buildFamId(genome)
        myMagSimusTools.printCurrentGenome(genome, famId, prefix=prefix)

    # Rename, remove and insert genes duplicated in tandem
    print >> sys.stderr, ("Renaming genes: conservation and duplications (clustering of tandem duplicated genes)")
    debugSet1 = set(rename.keys())
    debugSet2 = set([g.n for chrom in genome.values() for g in chrom])
    #if not debugSet1.issubset(debugSet2):
    #    raise ValueError('')
    assert debugSet1.issubset(debugSet2), "%s and rename|pblematic genes = %s" %\
        (debugSet1 - debugSet2, dict([(gn, rename[gn]) for gn in debugSet1 - debugSet2]))
    (genome, nbGeneEvents, dispersedTandemBlocks) =\
        renameAndTandemBlocks(genome, rename, indices, clusterator,
                              forceClustering, clusteredGenesRatio, tandemDupRatio,
                              prefix=prefix)
    nbDispersedTbs = len(dispersedTandemBlocks)
    nbTandemDup = nbGeneEvents['geneDup'] - nbDispersedTbs

    # DEBUG assertions
    nbDupTmp = sum([len(l)-1 for l in rename.values()])
    assert nbGeneEvents['geneDup'] == nbDupTmp, "%s == %s" % (nbGeneEvents['geneDup'], nbDupTmp)

    if printIntermGenomes:
        famId = myMagSimusTools.buildFamId(genome)
        myMagSimusTools.printCurrentGenome(genome, famId, prefix=prefix)

    # Insertion of blocks
    print >> sys.stderr, "Insertion of new genes: clustering and insertion of blocks"
    (genome, nbGeneEvents) = insertNewGenes(newGenes, dispersedTandemBlocks, genome, indices,
                                            forceClustering, clusteredGenesRatio, clusterator,
                                            nbGeneEvents, prefix=prefix)

    # DEBUG assertions
    # Conservation law of genes
    assert sum(len(chrom) for chrom in genome.values()) - nbIniGenes == nbGeneEvents['geneGain'] + nbGeneEvents['geneDup'] - nbGeneEvents['geneDel'], "%s - %s == %s + %s - %s" % (sum(len(chrom) for chrom in genome), nbIniGenes, nbGeneEvents['geneGain'], nbGeneEvents['geneDup'], nbGeneEvents['geneDel'])

    # DEBUG verification
    # Empty chromosomes areremoved
    for chrom in genome:
        if len(chrom) == 0:
            # FIXME: this should never happen, since this case should already be
            # fixed in renameAndTandemBlocks.
            raise

    listOfCs = genome.keys()
    for c in listOfCs:
        if len(genome[c]) == 0:
            del genome[c]

    nbGenes = sum(len(chrom) for chrom in genome.values())
    nbGeneEvents['geneLoss'] = nbIniGenes + nbGeneEvents['geneDup'] + nbGeneEvents['geneGain'] - nbGenes

    return (genome, nbGeneEvents, nbTandemDup)
