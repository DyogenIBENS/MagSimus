#!/usr/bin/python
# -*- coding: utf-8 -*-

__doc__ = """
        Return the future history of each gene along branches of the phylogenetic tree
"""

import sys

import utils.myFile as myFile
import utils.myTools as myTools
import utils.myPhylTree as myPhylTree
import utils.myLightGenomes as myLightGenomes

arguments = myTools.checkArgs([("speciesTree.phylTree", file)],
                              [("genesFile", str, "genesST.list.bz2"),
                               ("ancGenesFile", str, "ancGenes.list.bz2")],
                              __doc__)
# Load data
speciesTree = myPhylTree.PhylogeneticTree(arguments["speciesTree.phylTree"])
genomes = {}
for species in speciesTree.listSpecies:
    genomes[species] = myLightGenomes.LightGenome(arguments["genesFile"] % speciesTree.fileName[species], withDict=True)
for ancestor in speciesTree.listAncestr:
    genomes[ancestor] = myLightGenomes.Families(arguments["ancGenesFile"] % speciesTree.fileName[ancestor])


def genesOf(genome):
    res = set()
    if isinstance(genome, myLightGenomes.LightGenome):
        for c in genome:
            for g in genome[c]:
                res.add(g.n)
    elif isinstance(genome, myLightGenomes.Families):
        families = genome
        for family in families:
            res.add(family.fn)
    else:
        raise ValueError('genome should be either of type LightGenome or Families')
    return res


def do(node):
    for (child, _) in speciesTree.items.get(node, []):
        foo = "# Compute the %s -> %s branch #" % (node, child)
        bar = "#" * len(foo)
        print >> sys.stderr, bar
        print >> sys.stderr, foo
        print >> sys.stderr, bar
        trans = {}
        deleted = set()
        gained = genesOf(genomes[child])

        assert isinstance(genomes[node], myLightGenomes.Families)
        nbGenesParent = len([family for family in genomes[node]])
        nbDuplications = 0
        progressBar = myTools.ProgressBar(nbGenesParent)
        if isinstance(genomes[child], myLightGenomes.Families):
            for (i, family) in enumerate(genomes[node]):
                progressBar.printProgressIn(sys.stderr, i)
                # for each family name of 'node' we search descendant family names
                # in 'child'
                #assert genomes[child].getFamilyByName(family.fn, default=None) is None
                heritedG = set([])
                for n in family.dns:
                    # family name
                    f = genomes[child].getFamilyByName(n)
                    if f is not None:
                        heritedG.add(f.fn)
                heritedG = set(heritedG)
                if len(heritedG) == 0:
                    deleted.add(family.fn)
                else:
                    gained = gained - heritedG
                    trans[family.fn] = list(heritedG)
                    nbDuplications += len(heritedG)-1
            nbGenesSon = len([family for family in genomes[child]])
        elif isinstance(genomes[child], myLightGenomes.LightGenome):
            for (i, family) in enumerate(genomes[node]):
                progressBar.printProgressIn(sys.stderr, i)
                positions = []
                assert genomes[child].getPosition(family.fn, default=None) is None
                for gn in family.dns:
                    positions.append(genomes[child].getPosition(gn))
                heritedG = [genomes[child][p.c][p.idx].n for p in positions if p is not None]
                assert len(heritedG) == len(set(heritedG))
                # no doublon
                if len(heritedG) == 0:
                    deleted.add(family.fn)
                else:
                    gained = gained - set(heritedG)
                    trans[family.fn] = heritedG
                    nbDuplications += len(heritedG)-1
            nbGenesSon = sum(len(chrom) for chrom in genomes[child].values())
        else:
            raise ValueError('genomes[child] should either be a LightGenome or a Families')
        nbGenesHerited = sum(len(heritedG) for heritedG in trans.values())
        assert len(gained) + nbGenesHerited == nbGenesSon
        print >> sys.stderr, "%d genes -> %d genes" % (nbGenesParent, nbGenesSon)
        print >> sys.stderr, "%d herited (dup included), %d deleted, %d appeared" % (nbGenesHerited, len(deleted), len(gained))
        print >> sys.stderr, "%d Duplications (initial gene excluded) " % nbDuplications
        print >> sys.stdout, myFile.myTSV.printLine([trans, deleted, gained])
        do(child)
# main
assert isinstance(genomes[speciesTree.root], myLightGenomes.Families)
ini = set([family.fn for family in genomes[speciesTree.root]])
print >> sys.stderr, "%s: %d genes" % (speciesTree.root, len(ini))
print >> sys.stdout, ini
do(speciesTree.root)
