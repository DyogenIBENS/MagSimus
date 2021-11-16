#! /usr/bin/python
# -*- coding: utf-8 -*-

__doc__ = """
Simulate the evolution of an initial and artificial ancestral genome through
a species tree, knowing the gene events.
Events of evolution are:
    - genic events: de novo births, duplications (tandem and dispersed), losses
    - chromosomal rearrangements: fusions, fissions, translocations and inversions
Also simulates:
    - incomplete sequencing
"""

# TODO :
# 1) Make a choice between 'cluster' and 'block'. These two words are used for
# the same thing. This is because of the use of the verb 'to clusterise' when
# genes are gathered into blocks.
# 2) Remove 'changed by Joseph'
# 3) Remove 'M. Muffato'

import sys
import copy
import random

import bidict as bidict

import utils.myFile as myFile
import utils.myTools as myTools
import utils.myPhylTree as myPhylTree
import utils.myLightGenomes as myLightGenomes
from utils.myLightGenomes import OGene as Gene

import libs.myClustering as myClustering
import libs.myFixedGenicEvents as myFixedGenicEvents
import libs.myMagSimusTools as myMagSimusTools
import libs.myEvolutionProbas as myEvolutionProbas
import libs.myManagerOfEvents as myManagerOfEvents


# Performs additional operations on genome to simulate the non-perfect
# quality of 2X and 6X coverage sequenced genomes in real data.
#def applyCoverage(genome, coverage):
#    global prefix
#    if coverage == "2X":
#        # nbScaff = Number of scaffolds in function of the mean observed lengths
#        # N such as (1.1<=N<=2.9) is the mean length of scaffolds for 2X coverage
#        # nbScaff is equal to the nb of genes divided by N
#        nbScaff = int(sum([len(chrom) for chrom in genome.values()]) / random.uniform(1.1, 2.9))
#        # nbScaff - len(genome.keys()) = nb of scaffolds to create
#        for _ in xrange(nbScaff - len(genome.keys())):
#            print prefix, "2X-coverage", "fragmentation",
#            intervBreakStr = myChromosomalEvents.doChrBreak(genome)
#            print >> sys.stdout, intervBreakStr,
#            print
#        # There is only a part of the genome that is sequenced
#        # To remove chromosome randomly we first shuffle the genome
#        # random.shuffle shuffles a list and returns nothing
#        cs = genome.keys()
#        random.shuffle(cs)
#        newGenome = []
#        for c in cs:
#            newGenome.append(genome[c])
#        for (c, chrom) in zip(cs, newGenome):
#            genome[c] = chrom
#        # Only around 2/3 of the genome is sequenced
#        # Remove scaffolds that are between random.uniform(0.55, 0.85)th index
#        # and the end of the list of the scaffolds of the genome.
#        # This is equivalent to remove between roughly 45 and 15% of the
#        # scaffolds.
#        del genome[int(len(genome.keys()) * random.uniform(0.55, 0.85)):]
#        print >> sys.stderr, "OK (%d chromosomes)" % len(genome)
#    elif coverage == "6X":
#        # FIXME Should it be rather around 2000 and 2500 fissions ?
#        # Around 850 breaks
#        for _ in xrange(random.randint(800, 900)):
#            print prefix, "6X-coverage", "fragmentation",
#            intervBreakStr = myChromosomalEvents.doChrBreak(genome)
#            print >> sys.stdout, intervBreakStr,
#            print
#    return genome


# Chromosomal rearrangements
def performChromEvents(genome,
                       paramValues,
                       idx2eventType,
                       ChromSegVM_Mean=0.01,
                       ChromSegVM_Kappa=3,
                       fissionChrSamplingAlgo="simpleRandomSampling",
                       fusionChrSamplingAlgo="simpleRandomSampling",
                       inversionChrSamplingAlgo="simpleRandomSampling",
                       translocChrSamplingAlgo="simpleRandomSampling",
                       lowerCap=0,
                       upperCap=sys.maxint,
                       breakpointAnalyzer=True):
    setEvts = {"chrInvert", "chrTransloc", "chrFusion", "chrFission"}
    assert setEvts == set(idx2eventType.inv.keys())
    setParams = setEvts | {'geneTandemDupProp'}
    assert setParams == set(paramValues.keys())
    nbChrsIni = len(genome.keys())
    nbGenesIni = sum([len(chrom) for chrom in genome.values()])
    nbEffEvents = dict((idEvt, 0) for idEvt in idx2eventType.keys())

    # nbThEvents: nb theoretical events
    # nbThEvents contains 4 int values, each value corresponds to a number of
    # events for a specific chromosomal rearrangement.
    moe = myManagerOfEvents.ManagerOfEvents(paramValues,
                                            idx2eventType)

    for idEvt in moe:
        # Check that there is at least one non-empty chromosome
        # Using the implicit booleanness of the empty list
        # bool([]) = False
        # bool([1,2]) = True
        if bool(genome.keys()):
            if bool(genome[genome.keys()[0]]):
                pass
                # at least one gene in the first chromosome
            else:
                raise ValueError('First chromosome contains no gene')
        else:
            raise ValueError('No chromosome in the genome')

        # It is an inversion of a chromosomal segment
        if idEvt == moe.chrInvertId:
            (genome, flag) = moe.handleInversionEvent(genome,
                                                      ChromSegVM_Mean,
                                                      ChromSegVM_Kappa,
                                                      stream=sys.stdout)
        # It is a reciprocal translocation
        elif idEvt == moe.chrTranslocId:
            (genome, flag) = moe.handleReciprocalTranslocationEvent(genome,
                                                                    translocSamplingAlgo=translocSamplingAlgo,
                                                                    lowerCap=lowerCap,
                                                                    upperCap=upperCap,
                                                                    stream=sys.stdout)
        # It is a fusion
        elif idEvt == moe.chrFusionId:
            (genome, flag) = moe.handleFusionEvent(genome,
                                                   upperCap=sys.maxint,
                                                   stream=sys.stdout)
        # It is a fission
        elif idEvt == moe.chrFissionId:
            (genome, flag) = moe.handleFissionEvent(genome,
                                                    fissionSamplingAlgo=fissionSamplingAlgo,
                                                    lowerCap=lowerCap,
                                                    stream=sys.stdout)
        else:
            raise ValueError("Unknown idEvt, idEvt should be in %s" % idx2eventType.keys())

        if flag == 'OK':
            #print >> sys.stderr, "performed %s" % idx2eventType[idEvt]
            nbEffEvents[idEvt] += 1
        else:
            # If the event has not been performed, print the reason and how it
            # has been managed.
            print >> sys.stderr, flag

    nbEffEvents = dict([(idx2eventType[idEvt], nbEvt) for (idEvt, nbEvt) in nbEffEvents.iteritems()])
    # DEBUG assertion
    # Conservation law of chromosomes
    nbFusion = nbEffEvents['chrFusion']
    nbFission = nbEffEvents['chrFission']
    nbChrs = len(genome.keys())
    assert nbChrs - nbChrsIni == nbFission - nbFusion, "%s - %s == %s - %s" % (nbChrs, nbChrsIni, nbFission, nbFusion)
    # Conservation law of genes
    nbGenes = sum([len(chrom) for chrom in genome.values()])
    assert nbGenes == nbGenesIni, "%s = %s" % (nbGenes, nbGenesIni)
    return (genome, nbEffEvents)

# Recursive evolution/build of genomes
def evolve(speciesTree,
           node,
           genome,
           geneEventsReader,
           idx2eventType,
           arguments,
           branchSpecificParameters):
    # Performs additional operations on genome to simulate the non-perfect
    # quality of 2X and 6X coverage sequenced genomes in real data.
    #coverage = None
    #if node in speciesTree.lstEsp2X:
    #    coverage = "2X"
    #elif node in speciesTree.lstEsp6X:
    #    coverage = "6X"
    #if coverage is not None:
    #    print >> sys.stderr, "Applying %s coverage on %s ..." % (coverage, node),
    #    genome = applyCoverage(genome, coverage)
    #    print >> sys.stderr, "OK (%d chromosomes)" % len(genome.keys())

    # Writte the 'genome' of 'node' in a genes.<species>.list.bz2 file
    # 'families' is a dict that contains the gene names of 'node' and their
    # descendant gene names in children species.
    # For extant species (leaves of the species tree), 'families' contains the
    # extant gene and itself considered as a child gene.
    # 'Families' is a dict that contains the gene names of genome and their
    # children gene names in the species that descent from this genome.
    families = myLightGenomes.Families()
    print >> sys.stderr, "Writing %s genome ..." % node,
    genome.printIn(myTools.myFile.openFile(arguments["out:genomeFiles"] % speciesTree.fileName[node], 'w'))
    print >> sys.stderr, "OK"
    for chrom in genome.values():
        for g in chrom:
            families.addFamily(myLightGenomes.Family(g.n, [g.n] if node in speciesTree.listSpecies else []))
    # Used only by 'buildFamId' and 'printCurrentGenome'
    # FIXME is it used ?
    # global famid

    # If 'node' is not a leaf <=> if it has children
    if node in speciesTree.items:

        # 'child': name of the child species of 'node'
        for (child, bLength) in speciesTree.items[node]:
            # utilise la function __repr__ du decorateur memoize pour afficher
            # le nb de valeurs enregistrees en cache, le nombre d'appels et
            # l'acceleration due a memoize.
            #print >> sys.stderr, coevoluTionScore

            # Deep copy of the genome
            # This copy will be modified and the genome of 'node' will not
            # change.
            newGenome = copy.deepcopy(genome)
            # The genome evolves along a branch with a given duration and
            # chromosomal rearrangement rates fixed by the user.
            print >> sys.stdout, "BRANCH %s %s" % (speciesTree.fileName[node],
                                                   speciesTree.fileName[child])
            speciesTree.printBranchName(child, stream=sys.stderr)
            paramValues = branchSpecificParameters.paramVal[child]


            ###################
            # 1) Genic Events #
            ###################
            # 1- re-name genes:
            #   1.a- terminal branches: pass from the ENSGT00 naming to
            #   the ENSG00 / ENSMUSG00 / etc naming.
            #   1.b- on intern branches (not leading to leaves): duplications
            # 2- remove lost genes (= not presents in 'rename')
            # 3- add new gene (in 'newGenes')
            (rename, newGenes) =\
                myFixedGenicEvents.fixGeneEvents(newGenome, child, bLength, speciesTree,
                                            geneEventsReader)
            # List of the indices of the branches for the coevolutionScore
            # computations.
            indicesDup = []
            indicesDel = []
            nbOfBranches = len(speciesTree.indBranches)
            for e in speciesTree.allNames:
                if speciesTree.isChildOf(e, child) and (child != e):
                    indicesDup.append(speciesTree.indBranches[e])
                    indicesDel.append(nbOfBranches + speciesTree.indBranches[e])
            indices = frozenset(indicesDup + indicesDel)
            # Rq: The gene events vector is equal to DupG + DelG, and has a total
            # of 8+8 = 16 indices.
            # Re-initialize the cache of the memoize decorator
            clusterator.coevoluTionScore.reinit_stats()

            # Perform genic events in a random order
            (newGenome, nbGeneEvents, nbDupTandem) =\
                myFixedGenicEvents.performGeneEvents(rename, newGenes, newGenome,
                                                indices,
                                                arguments['forceClustering'],
                                                arguments['geneClusteredRatio'],
                                                paramValues['geneTandemDupProp'],
                                                clusterator,
                                                printIntermGenomes=arguments['printIntermGenomes'])
            # DEBUG asserion
            #assert len(newGenome) > 0, "node %s -> child %s" % (node, child)
            #assert len(newGenome[0]) > 0, "node %s -> child %s" % (node, child)

            table = []
            for evt in ['geneGain', 'geneDup', 'geneLoss']:
                nbEvt = nbGeneEvents[evt]
                table.append(["#%s" % evt,
                              "= %s" % nbEvt,
                              "rate = %.2f" % (float(nbEvt)/bLength)])
            table.append(["propTandemDup",
                          "= %.2f%%" % (100 * ( float(nbDupTandem) / float(nbGeneEvents['geneDup']) )),
                          ""])
            table.append(["#genes",
                          "= %s" % sum(len(chrom) for chrom in newGenome.values()),
                          ""])
            myTools.printTable(table, sys.stderr)

            #################################
            # 2) Chromosomal rearrangements #
            #################################
            # Compute the nb of chromosomal rearrangement for each type of
            # chromosomal rearrangement.
            print >> sys.stderr, "Compute and perform chromosomic events"

            (newGenome, nbChromEvents) = \
                performChromEvents(
                    newGenome,
                    paramValues,
                    idx2eventType,
                    ChromSegVM_Mean=arguments['chr:invDistMean'],
                    ChromSegVM_Kappa=arguments['chr:invDistShape'],
                    lowerCap=arguments['chr:sizeLowerCap'],
                    upperCap=arguments['chr:sizeUpperCap'],
                    fissionChrSamplingAlgo=arguments['chr:fissionSamplingAlgo'],
                    fusionChrSamplingAlgo="simpleRandomSampling",
                    inversionChrSamplingAlgo="simpleRandomSampling",
                    translocChrSamplingAlgo=arguments['chr:translocSamplingAlgo'],
                    onlyDetectableEvents=arguments['passNonViableEvents'],
                    printIntermGenomes=arguments['printIntermGenomes'])

            #if [nbInv, nbTransloc, nbFus, nbBreak] != nbEvents:
            #    print >> sys.stderr,\
            #        ("Warning, the expected number of chromosomal events differ "
            #         "from the number of chromosomal rearrangements actually "
            #         "performed")

            table = []
            for (_, evt) in sorted(idx2eventType.iteritems(), key=lambda x:x[0]):
                nbEvts = nbChromEvents[evt]
                table.append(["#%s" % evt,
                              "= %s" % nbEvts,
                              "rate = %.2f" % (float(nbEvts)/bLength),
                              "randomFactor = %.2f" % branchSpecificParameters.randomFactors[child][evt]])
            table.append(["#chromosomes",
                          "= %s" % len(newGenome),
                          "",
                          ""])
            myTools.printTable(table, sys.stderr)

            # Recursive call, continue evolution
            subFam = evolve(speciesTree, child, newGenome, geneEventsReader, idx2eventType, arguments, branchSpecificParameters)

            # Descendant genes in child species are added in families.
            # For genes that are conserved from parent to child
            for (parentGn, childGns) in rename.iteritems():
                # For all genes of the parent (with the gene mames of the child
                # genes)
                for gn in childGns:
                    # If 'gn' is in a family of child
                    subfamily = subFam.getFamilyByName(gn, None)
                    if subfamily is not None:
                        family = families.getFamilyByName(parentGn, None)
                        if family is not None:
                            family = myLightGenomes.Family(family.fn, list(set(family.dns + subfamily.dns)))
    # Write the list of ancestral genes
    print >> sys.stderr, "Writing %s ancestral genes ..." % node,
    families.printIn(myFile.openFile(arguments["out:ancGenesFiles"] % speciesTree.fileName[node], 'w'))
    print >> sys.stderr, "OK"
    return families

# Design an initial ancestral genome.
# Since MagSimus2 know "a priori" the gene events thanks to the input gene trees,
# the location of genes is important and have a huge impact on the disruption
# level of gene order during simulation.
# For instance, if two genes have a coevolution (let's say a deletion) that
# occurs on the same branch. If they are next to each other, only one adjacency
# changes, otherwise two adjacencies are disrupted.
def randomGenomeWithFixedIDs(geneNames, forceClustering, clusteredGenesRatio,
                             clusterator, speciesTree):
    # DEBUG assertion
    # verify that geneNames is strictly included in byname
    #assert set(str(x) for x in geneNames).issubset(byname), set(str(x) for x in geneNames).difference(byname)

    # FIXME use only indices of deletions for coevolutionScore calculations.
    bIndices = speciesTree.indBranches.values()
    #indicesDup = bIndices
    indicesDel = [val + len(bIndices) for val in bIndices]
    # The events vector is equal to [DelG], and has a total of 8 indices
    #indices = frozenset(indicesDup + indicesDel)
    indices = frozenset(indicesDel)

    print >> sys.stderr, ("Gene clusterisation of the 1st ancestral genome,"
                          " this task is the longest of the simulation."
                          " It can take some time ...")
    # Clusterisation considering all branch indices of the species tree.
    # The first designed ancestral genome is the root of the species tree.
    blocks = clusterator.mkClusterByScore(geneNames, indices,
                                          forceClustering, clusteredGenesRatio,
                                          verbose=True)
    print >> sys.stderr, "End of the clusterisation!"
    print >> sys.stderr, "Nb of blocks of coevolving genes for the artificial design of the first ancestral genome = %s" % len(blocks)
    blocks = [[g for (m, g) in block] for block in blocks]
    random.shuffle(blocks)
    genes = [gn for block in blocks for gn in block]
    return genes

if __name__ == '__main__':
    argsFormatList = [
        ###############
        # Input Files #
        ###############
        # Species tree
        ("speciesTree.phylTree", file),
        # File with the genic events: rename, deleted, gained
        ("gene:geneEvents", file),
        # File with the vectors of gene evolution (duplications and deletions)
        ("gene:timeline", file)
        ]
    optsFormatList = [
        # File with listed parameters and values
        ("parameterFile", str, ""),
        # File with custom branch specific and chrom-rearrangement specific factors
        ("userRatesFile", str, ""),
        # Global factor of chrom-rearrangements
        # If it is 2 times higher there are 2 times more chrom-rearrangements
        ("globalFactor", float, 1.),
        #################################
        # initial genome configuration  #
        #################################
        ("startChrSizes", tuple, (0.141,0.095,0.081,0.074,0.064,0.054,0.036,0.035,0.033,0.03,0.029,0.028,0.025,0.024,0.024,0.023,0.023,0.022,0.021,0.02,0.017,0.017,0.016,0.016,0.016,0.013,0.012,0.008,0.003,0.001)),
        ###################################
        # Chromosomal rearrangement rates #
        ###################################
        # Factors specific to each type of chrom-rearrangement
        # Warning: these names should correspond with the names of events in
        # idx2evenType and also with the names of events in the userRatesFile
        ("chrInvert", float, 1.),
        ("chrTransloc", float, 1.),
        ("chrFusion", float, 1.),
        ("chrFission", float, 1.),
        ###################################################
        # Distribution of rearranged chromosomal segments #
        ###################################################
        # Choose sampling algorithm for translocations
        # Alternatives are "lengthPropSampling" or "simpleRandomSampling"
        ("chr:translocSamplingAlgo", str, "simpleRandomSampling"),
        # Choose sampling algorithm for Fissions
        # Alternatives are "lengthPropSampling" or "simpleRandomSampling"
        ("chr:fissionSamplingAlgo", str, "simpleRandomSampling"),
        # Non viable events are passed if True
        # If False, all events specified in specRates are effective events
        # recorded in the lineage and extant genomes have a constant nb of
        # chromosomes in all simulations.
        ("passNonViableEvents", bool, False),
        # Lower cap for chromosomes in genes
        # If chromosome size falls below due to fission or translocation,
        # the event is either postponed or ignored
        # (depending on chr:forceChrNumber)
        ("chr:sizeLowerCap", int, 1),
        # Upper cap for chromosomes in genes
        # If exceeded due to fusion or translocation, the event is either
        # postponed or ignored (depending on chr:forceChrNumber)
        ("chr:sizeUpperCap", int, 1000000),
        # Parameters for the shape of the distribution of the lengths of
        # rearranged segments (for inversions and translocations)
        # The distribution is defined in myMaths.randomValue.myVonMise.
        # It uses two parameters.
        # Pseudo-mean of the modified von-Mises
        # FIXME: solve the issue of the mean
        ("chr:invDistMean", float, 0.010),
        # Concentration of the bell-shape of the von Mises
        # The higher Kappa, the higher and thinner, the bell-shape.
        ("chr:invDistShape", float, 3.0),
        ###################
        # Gene clustering #
        ###################
        # Ratio of the desired number of clusterised genes
        # nb of clusterised genes / nb of genes candidate to the clustering
        # Rq, 'clusteredGenesRatio' is dominant over 'forceClustering'.
        # I.e. geneClusteredRatio=0.75 => ~75% of the genes are clusterised even
        # if 'forceClustering' = True
        ("geneClusteredRatio", float, 0.99),
        # Ratio of expected tandem duplications
        ("geneTandemDupProp", float, 1.0),
        # Force the clustering of genes that have a scoreSameEvent = 0 one after
        # the other, after genes already clusterised with the gene of referencee.
        # Clusterise also genes without history at the level of extant species.
        # Rq : if arguments["forceClustering"] and
        #       arguments["gene:clusteredGenesRatio"] == 1:
        #       mkClusterByScore returns blocks of len = 1
        #       In this case the clustering return one block with all the genes.
        ("forceClustering", bool, True),
        # Seed to initialise the random generator
        ("seed", str, ""),
        #############################
        # Random variation of rates #
        #############################
        ("b_randomAccel", bool, True),
        # Nominale value, difference from this value are caused by the random
        # factor.
        ("rate:eventMaxRandomFactor", float, 1.5),
        # Concentration of the bell-shape of the von Mises
        # The higher Kappa, the higher and thinner, the bell-shape.
        ("rate:eventRandomVonMisesKappa", float, 2.),

        # Parameters of the former gamma function.
        # The gamma function was historically used before the modified
        # von-Mises.
        #("chr:invertGammaAlpha",float,.2),
        #("chr:invertGammaBeta",float,2.),

        # Boolean to write intermediary genomes during the evolution along a
        # branch.
        ("printIntermGenomes", bool, False),

        ################
        # Output Files #
        ################
        # Files containing the simulated genomes at each node and leaf of the
        # species tree.
        ("out:genomeFiles", str, "simu/genes/genes.%s.list.bz2"),
        # File containing gene families at each node and leaf of the
        # species tree.
        ("out:ancGenesFiles", str, "simu/ancGenes/ancGenes.%s.list.bz2"),

        # File to load an initial ancestral genome instead of designing it.
        # Since the clustering process takes time, his can be used for a faster
        # execution after the design of the first ancestor.
        ("in:firstAncestralGenome", str, None)
        ]

    # do not output the first arguments before loading parameter file
    arguments = myTools.checkArgs(argsFormatList, optsFormatList, __doc__, showArgs=False)

    # Read parameter file
    arguments = myMagSimusTools.readParameterFile(arguments)
    myTools.printArguments(arguments, stream=sys.stderr)

    assert arguments["geneClusteredRatio"] <= 1.0, "geneClusteredRatio must be between 0 and 1"

    speciesTree = myPhylTree.PhylogeneticTree(arguments["speciesTree.phylTree"])

    print >> sys.stdout, "INITIAL ANCESTRAL GENOME %s" % (speciesTree.root)

    if len(arguments["seed"]) > 0:
        random.setstate(eval(arguments["seed"]))
        print >> sys.stderr, "Random number generator internal state initialised from arguments"
    # Information about the random seed
    #print >> sys.stderr, "Random number generator internal state:", random.getstate()

    geneEventsReader = myFile.openFile(arguments["gene:geneEvents"], "r")

    #This must match with arguments names
    idx2eventType = bidict.bidict({1: "chrInvert",
                                   2: "chrTransloc",
                                   3: "chrFusion",
                                   4: "chrFission"})

    clusterator = myClustering.Clusterator(arguments["gene:timeline"])

    if arguments["in:firstAncestralGenome"] is not None:
        print >> sys.stderr, "Load first ancestral genome"
        # 1st line of 'geneEvents': set(['gn1','gn2',...])
        # Do nothing with it, It is just to position the cursor on the second line.
        geneEventsReader.readline()
        genome = myLightGenomes.LightGenome(arguments["in:firstAncestralGenome"], withDict=False)
    else:
        print >> sys.stderr, "Loading %s genes" % speciesTree.root
        print >> sys.stderr, "Design the first ancestral genome"
        chromLengths = []
        startChrSizes = list(float(lengthChrom) for lengthChrom in arguments['startChrSizes'])
        # In any case that the user put them in sorted order
        # Since the most clustered genes come first in listOfClusteredGenes,
        # this could introduce a bias... By shuffling the distribution of
        # chromosomes we overcome this possible bias.
        random.shuffle(startChrSizes)
        # Read the first line of 'geneEvents', format: set(['gn1','gn2',...])
        iniGeneNames = eval(geneEventsReader.readline())
        nbGenesIni = len(iniGeneNames)
        listOfClusteredGenes = randomGenomeWithFixedIDs(iniGeneNames,
                                                        arguments['forceClustering'],
                                                        arguments['geneClusteredRatio'],
                                                        clusterator,
                                                        speciesTree)
        #from IPython import embed
        #embed()
        assert set(listOfClusteredGenes) == iniGeneNames
        print >> sys.stderr, "End of clusterization nb of vectors  = %s " % len(clusterator.allvect)
        # Random chromosome lengths
        # chromLengths = myEvolutionProbas.randomSizes(arguments["iniChrNumber"],
        #                                        len(listOfClusteredGenes), 1)
        for chromLength in startChrSizes:
            chromLengths.append(int(len(listOfClusteredGenes) * chromLength/sum(arguments['startChrSizes'])))
        # warning, the length of the last chromosome may vary a bit from the
        # user input specified in arguments['startChrSizes']
        chromLengths[-1] += len(listOfClusteredGenes) - sum(chromLengths)
        assert len(listOfClusteredGenes) == sum([chromLength for chromLength in chromLengths])

        # Genes are put into chromosomes
        genome = myLightGenomes.LightGenome(withDict=False)
        for (c, nbGenes) in enumerate(chromLengths):
            chromosome = listOfClusteredGenes[:nbGenes]
            chromosome = [OGene(gn, myEvolutionProbas.randomStrand()) for gn in chromosome]
            genome[str(c)] = chromosome
            del listOfClusteredGenes[:nbGenes]
        assert len(listOfClusteredGenes) == 0
        assert sum([len(chrom) for chrom in genome.values()]) == nbGenesIni

    print >> sys.stderr, "#genes in %s = %s" %\
        (speciesTree.root, sum([len(chrom) for chrom in genome.values()]))
    print >> sys.stderr, "#chromosomes in %s = %s" %\
        (speciesTree.root, len(genome.keys()))
    print >> sys.stderr, "Distribution of chromosome lengths in %s = %s" %\
        (speciesTree.root, sorted([len(chrom) for chrom in genome.values()]))

    # event specific parameter value
    evtSpeParamValue = {}
    for param in idx2eventType.values() + ['geneTandemDupProp']:
        evtSpeParamValue[param] = arguments[param]

    branchSpecificParameters = myMagSimusTools.BranchParameters(speciesTree,
                                                                evtSpeParamValue,
                                                                arguments['rate:eventRandomVonMisesKappa'],
                                                                paramBranchSpeParamFile=arguments['userRatesFile'],
                                                                b_randomAccel=arguments['b_randomAccel'],
                                                                eventMaxRandomFactor=arguments['rate:eventMaxRandomFactor'],
                                                                globalFactor=arguments['globalFactor'])

    print >> sys.stderr, "Start simulation of the evolution"
    # Start evolution in silico
    evolve(speciesTree,
           speciesTree.root,
           genome,
           geneEventsReader,
           idx2eventType,
           arguments,
           branchSpecificParameters)
