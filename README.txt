# -*- coding: utf-8 -*-
# python 2.7 at least is needed
#
#                   MagSimus v1.00 and MagSimus v2.00
#
# This code may be freely distributed and modified under the terms of the GNU General Public License version 3 or later and CeCILL. This should be distributed with the code. If you do not have a copy, see:
#
#      http://www.gnu.org/licenses/gpl-3.0-standalone.html
#      http://www.cecill.info/licences/Licence_CeCILL_V2-en.html
#
# Copyright for this code is held jointly by the Dyogen (DYnamic and Organisation of GENomes) team of the Institut de Biologie de l'Ecole Normale Supérieure and the individual authors. 
#
# Copyright © 2013 IBENS/Dyogen Joseph LUCAS, Matthieu MUFFATO and Hugues ROEST CROLLIUS
# mail : hrc@ens.fr or jlucas@ens.fr
# This is free software, you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 (GPL v3) or later and the CeCiLL v2 license in France
#
# This code is based on
# TODO <publication>
# authors :
# Joseph LUCAS (IBENS, Paris, France, jlucas@ens.fr)
# Matthieu MUFFATO (EBI, Cambridge, United Kingdom, muffato@ebi.ac.uk)
# Hugues ROEST CROLLIUS (IBENS, Paris, France, hrc@ens.fr)

#####################################
# MagSimus v1.00 and MagSimus v2.00 #
#####################################

###########
# INSTALL #
###########

Install python 2.7:
: sudo apt-get install python2.7

Define the folder that will contain all the project
: INSTALLDIR=<pathToInstallDir> # default ${HOME}/Libs

Download LibsDyogen and the PhylDiag projects that are used by magSimus scripts
: mkdir ${INSTALLDIR}
: cd ${INSTALLDIR}
: git clone git@depot.biologie.ens.fr:LibsDyogen
: echo "export PYTHONPATH=\"$PYTHONPATH:${INSTALLDIR}/LibsDyogen\"" >> ~/.bashrc
: git clone git@depot.biologie.ens.fr:PhylDiag
: chmod +x PhylDiag/src/*.py
: git clone git@depot.biologie.ens.fr:MagSimus
: chmod +x MagSimus/src/*.py
: chmod +x MagSimus/src/analys/*.py

From now on we consider that the user is in the root directory of the project
: cd ${INSTALLDIR}/MagSimus

################################################################
# Preprocessing step 1 : define gene families using gene trees #
################################################################

Convert nhx (or .nwk, newick) gene trees to our tabular format (phylTree):
: src/nhxGeneTrees2phylTreeGeneTrees.py data/geneTrees.example.nhx > res/geneTrees.protTree

Convert a newick species tree into a phylTree species tree:
: src/newickSpeciesTree2phylTreeSpeciesTree.py data/speciesTree.nwk > res/speciesTree.phylTree

Extract the ancestral gene content (ancGene) from the gene trees:
: src/ancGenesFromGeneTrees.py res/speciesTree.phylTree res/geneTrees.protTree -out:ancGenes=res/ancGenes.example.%s.list.bz2 > res/geneTrees.afterExtractingAncGenes.protTree

These ancGenes files can be used to define gene families.

########################
# Processing MagSimus1 #
########################

Access to the manual of MagSimus1
: src/magSimus1.py

This returns :
//////////////////////
Usage : ./src/magSimus1.py
	1: phylTree.conf <type 'file'>
	2: root <type 'str'>
	3: iniGeneNumber <type 'int'>
	4: iniChrNumber <type 'int'>
	+/-b_randomAccel (True)
	  -rate:eventMaxAccel <type 'float'> (1.5)
	  -rate:eventVonMisesKappa <type 'float'> (2.0)
	  -chr:CsegmentVonMisesMean <type 'float'> (0.01)
	  -chr:CsegmentVonMisesKappa <type 'float'> (3.0)
	  -globalFactor <type 'float'> (1.0)
	  -userRatesFile <type 'str'> ()
	  -geneLoss <type 'float'> (6.0)
	  -geneGain <type 'float'> (1.0)
	  -geneDuplication <type 'float'> (5.0)
	  -geneTandemDupProp <type 'float'> (0.95)
	  -chrInvert <type 'float'> (2.0)
	  -chrTransloc <type 'float'> (0.2)
	  -chrFusion <type 'float'> (0.3)
	  -chrBreak <type 'float'> (0.3)
	  -duplicationGeneSwap <type 'float'> (0.0)
	  -out:genomeFiles <type 'str'> (simu/genes/genes.%s.list.bz2)
	  -out:ancGenesFiles <type 'str'> (simu/ancGenes/ancGenes.%s.list.bz2)


Takes a species tree and simulate the evolution of an artificial ancestral genome
Yield files similar to Ensembl files
-> Take into account specific evolution/errors of sequencing:
    - chromosomal rearrangements constraints
    - incomplete sequencing
    - 2X coverage
//////////////////////

# wrong parameters, not enough chromosomes
#: src/magSimus1.py res/speciesTree.phylTree Amniota 20000 25 -out:genesFiles=res/simu1/genes.%s.list.bz2 -out:ancGenesFiles=res/simu1/ancGenes.%s.list.bz2 > res/simu1/logStd 2> res/simu2/logErr & tail -f --pid=$! -s0.1 res/simu1/logErr

With specified user rates of chromosomal rearrangements
: src/magSimus1.py res/speciesTree.phylTree Amniota 20000 25 -userRatesFile=data/specRates.v74 -out:genesFiles=res/simu1/genes.%s.list.bz2 -out:ancGenesFiles=res/simu1/ancGenes.%s.list.bz2 > res/simu1/logStd 2> res/simu2/logErr & tail -f --pid=$! -s0.1 res/simu1/logErr

#################################################
# Preprocessing step 2 specific to MagSimus2 :  #
#################################################

Calc the evolution of gene content along the species tree
: src/calcGeneEvolution.py res/speciesTree.phylTree -genesFile=data/genesST.%s.list.bz2 -ancGenesFile=data/ancGenes.%s.list.bz2 > res/simu2/geneEvents.txt 2> res/simu2/geneEvents.log

Calc the vectors of evolution of each gene
: src/getGeneTimeline.py res/speciesTree.phylTree res/simu2/geneEvents.txt > res/simu2/geneTimeline.txt 2> res/simu2/geneTimeline.log

########################
# Processing MagSimus2 #
########################

Access to the manual of MagSimus2
: src/magSimus2.py

This returns :
//////////////////////
 Usage : ./src/magSimus2.py
	1: phylTree.conf <type 'file'>
	2: iniChrNumber <type 'int'>
	3: gene:geneEvents <type 'file'>
	4: gene:timeline <type 'file'>
	  -userRatesFile <type 'str'> (res/specRates)
	  -globalFactor <type 'float'> (1.0)
	  -chrInvertRate <type 'float'> (2.0)
	  -chrTranslocRate <type 'float'> (0.2)
	  -chrFusionRate <type 'float'> (0.2)
	  -chrBreakRate <type 'float'> (0.2)
	  -gene:clusteredGenesRatio <type 'float'> (0.99)
	  -gene:tandemDupRatio <type 'float'> (0.6)
	+/-forceClustering (True)
	  -seed <type 'str'> ()
	+/-b_randomAccel (True)
	  -rate:eventMaxRandomFactor <type 'float'> (1.5)
	  -rate:eventRandomVonMisesKappa <type 'float'> (2.0)
	  -chr:invDistMean <type 'float'> (0.01)
	  -chr:invDistShape <type 'float'> (3.0)
	+/-printIntermGenomes (False)
	  -out:genomeFiles <type 'str'> (simu/genes/genes.%s.list.bz2)
	  -out:ancGenesFiles <type 'str'> (simu/ancGenes/ancGenes.%s.list.bz2)
	  -in:firstAncestralGenome <type 'str'> (None)


        Simulate chromosomal rearrangements of an artificial ancestral genome knowing the gene events)
//////////////////////

Classical example
: src/magSimus2.py res/speciesTree.phylTree 25 geneEvents.txt  geneTimeline.txt -userRatesFile=data/specRates.v74 -out:genomeFiles=res/simu2/genes.%s.list.bz2 -out:ancGenesFiles=res/simu2/ancGenes.%s.list.bz2 > res/simu2/logStd 2> res/simu2/logErr &tail -f --pid=$! -s0.1 res/simu2/logErr

With a given first ancestral genome
: src/magSimus2.py res/speciesTree.phylTree 20 res/simu2/geneEvents.txt res/simu2/geneTimeline.txt -userRatesFile=data/specRates.v74 -out:genomeFiles=res/simu2/genes.%s.list.bz2 -out:ancGenesFiles=res/simu2/ancGenes.%s.list.bz2 -in:firstAncestralGenome=res/simu2/genes.Amniota.list.bz2 > res/simu2/logStd_withGiven1stAncestralGenome 2> res/simu2/logErr & tail -f --pid=$! -s0.1 res/simu2/logErr


#TODO: update with the version of ./checkMagSimusIntegrity.sh once the pipeline is ok

###########################
## Launch 100 simulations #
###########################
#
#Define folders for results and binaries
#: MagSimus=${INSTALLDIR}/MagSimus/src
#: PhylDiag=${INSTALLDIR}/PhylDiag/src
#Be sure that the hard drive has enough space (typically around !!!TODO!!! Mo)
#: Res=<PathToSimus>
#
#Copy input data in the res folder to easily find them
#: cp ${INSTALLDIR}/MagSimus/data/speciesTreeCoverage.phylTree ${Res}/speciesTree.phylTree
#: cp ${INSTALLDIR}/MagSimus/data/specRates.v74 ${Res}/
#: ln -s /kingdoms/dyogen/workspace4/workspace4/jlucas/Workflow/Data/GenesWithOnlyGoodChrs/Data72 ${Res}/GenesST
#: ln -s /kingdoms/dyogen/workspace2/workspace2/alouis/GENOMICUS_SVN/data72/trees/0.30/ancGenes/all ${Res}/AncGenes
#
#Calc geneEvents.txt
#: ${MagSimus}/calcGeneEvolution.py ${Res}/speciesTree.phylTree -genesFile=${Res}/GenesST/genesST.%s.list.bz2 -ancGenesFile=${Res}/AncGenes/ancGenes.%s.list.bz2 > ${Res}/geneEvents.txt 2> ${Res}/geneEvents.log
#
#Calc geneTimeline.txt
#: ${MagSimus}/getGeneTimeline.py ${Res}/speciesTree.phylTree ${Res}/geneEvents.txt > ${Res}/geneTimeline.txt 2> ${Res}/geneTimeline.log
#
#Launch 100 simulations using condor
#: for i in `seq 0 99`;do mkdir -p ${Res}/$i && condor-submit.sh ${MagSimus}/magSimus2.py ${Res}/speciesTree.phylTree 20 ${Res}/geneEvents.txt ${Res}/geneTimeline.txt -userRatesFile=${Res}/specRates.v74 -out:genomeFiles=${Res}/${i}/genes.%s.list.bz2 -out:ancGenesFiles=${Res}/${i}/ancGenes.%s.list.bz2 +b_randomAccel / ${Res}/${i}/logStd // ${Res}/${i}/logErr /// ${Res}/logCondor =; done | $jl/Libs/Scripts/condor-clean.py > ${Res}/COND
#: condor_submit ${Res}/COND
#
#To monitor simulations:
#: for i in `seq 10000`; do date ; condor_q|grep -B4 running; sleep 1m; done > foo & tail -f foo
#
#Check that no simulation has been aborted, because 'too many fusions. Not enough chromosomes':
#: grep Error ${Res}/*/logErr
#
#Otherwise launch:
#: for i in `grep ValueErr ${Res}/*/logErr|cut -d/ -f5|sort|uniq`; do condor-submit.sh ${MagSimus}/magSimus2.py ${Res}/speciesTree.phylTree 20 ${Res}/geneEvents.txt ${Res}/geneTimeline.txt -userRatesFile=${Res}/specRates.v74 -out:genomeFiles=${Res}/${i}/genes.%s.list.bz2 -out:ancGenesFiles=${Res}/${i}/ancGenes.%s.list.bz2 +b_randomAccel / ${Res}/${i}/logStd // ${Res}/${i}/logErr /// ${Res}/logCondor =; done | $jl/Libs/Scripts/condor-clean.py > ${Res}/COND2
#: condor_submit ${Res}/COND2
#
#Repeat this last step untill all simulations finish.
#
##################
##	Post-analysis	#
##################
#
## Distribution of chromosome lengths of the Human:
###################################################
#: S=Homo.sapiens && ${MagSimus}/analysis/distributionChrLengths_Real_vs_Simus.py ${Res}/GenesST/genesST.${S}.list.bz2 -in:simulatedGenomes=${Res}/%s/genes.${S}.list.bz2 > ${Res}/distributionChrLengths_Simus_vs_Real_${S}.svg
#
## Extract synteny blocks
#########################
#Real case
#---------
#Extract synteny blocks in pairwise comparisons of real extant species using phylDiag and the cluster of the IBENS:
#: Libs/MagSimus/src/analysis/extractSbs_Real_S1_vs_S2.bsh -g 5 -d CD ${Res}/GenesST ${Res}/AncGenes ${MagSimus}/analysis/speciesCombi.txt ${Res}
#
#Simulations
#-----------
#Extract synteny blocks in pairwise comparisons of simulated extant species using phylDiag and the cluster of the IBENS
#: Libs/MagSimus/src/analysis/extractSbs_Simus_S1_vs_S2.bsh -g 5 -d CD ${Res} ${Res}/AncGenes ${MagSimus}/analysis/speciesCombi.txt ${Res}
#
#'-g 5' means that the gapMax is set to 5
#'-d CD' means that the distance metric used is the Chebyshev Distance Metric
#
#FIXME: verify that the cluster worked properly, often, some cores have troubles and some output files are lacking
#
## N50 of the pairwise comparisons (boxplot with matplotlib)
############################################################
#: Libs/MagSimus/src/analysis/graph_N50_Real_MagSimus.py ${MagSimus}/analysis/speciesCombi.txt -realDiagsLog=${Res}/Real_%s_vs_%s.sbs.log -simulatedDiagsLog=${Res}/%s/%s_vs_%s.sbs.log > ${Res}/graph_N50_Real_MagSimus.svg
#
## Distribution of the lengths of synteny blocks in the comparison Human-Mouse
##############################################################################
#: S1=Homo.sapiens && S2=Mus.musculus && Libs/MagSimus/src/analysis/distribSbsLengths_Real_vs_Simu.py ${Res}/Real_${S1}_vs_${S2}.sbs -in:Sim%s_S1_vs_S2.sbs=${Res}/%s/${S1}_vs_${S2}.sbs > ${Res}/distribSbsLengths_${S1}_vs_${S2}.svg
#
##Whole genome comparison
#########################
#Homology Matrix
#---------------
#Real case:
#: Libs/MagSimus/src/analysis/misc.compareGenomes.py ${Res}/GenesST/genesST.Homo.sapiens.list.bz2 ${Res}/GenesST/genesST.Mus.musculus.list.bz2 ${Res}/AncGenes/ancGenes.Euarchontoglires.list.bz2 +sortBySize -matrix:pointSize=0.05 > ${Res}/DotMatrix_Real_Human_vs_Mouse_pS0.05.ps
#
#Remark on how to convert ps files into pdf files:
#: for i in `ls ./Do*.ps`; do filename=$(basename "${i}"); extension="${filename##*.}"; filename="${filename%.*}"; ps2pdf -dAutoRotatePages=/None ${filename}.ps ${filename}.pdf; done
#The conversion often invert the x-axis from bottom to top... :( to avoid this add the -dAutoRotatePages=/None option.
#
#Simulation:
#: i=0 && Libs/MagSimus/src/analysis/misc.compareGenomes.py ${Res}/${i}/genes.Mus.musculus.list.bz2 ${Res}/${i}/genes.Homo.sapiens.list.bz2  ${Res}/AncGenes/ancGenes.Euarchontoglires.list.bz2 +sortBySize -matrix:pointSize=0.05 > ${Res}/DotMatrix_Human_vs_Mouse_pS0.05_tandemDupProp0.60.ps
#
#Karyotype:
#----------
#Real:
#: S1=Homo.sapiens && S2=Mus.musculus && i=SimuTest_tandemDupProp0.8 && Libs/MagSimus/src/analysis/misc.compareGenomes.py ${Res}/GenesST/genesST.${S1}.list.bz2 ${Res}/GenesST/genesST.${S2}.list.bz2 ${Res}/AncGenes/ancGenes.Euarchontoglires.list.bz2 -mode=drawKaryotype +sortBySize > ${Res}/Karyotype_Real_${S1}_${S2}_tandemDupProp0.8.ps
#
#Simulation:
#: S1=Homo.sapiens && S2=Mus.musculus && i=SimuTest_tandemDupProp0.8 && Libs/MagSimus/src/analysis/misc.compareGenomes.py ${Res}/${i}/genes.${S1}.list.bz2 ${Res}/${i}/genes.${S2}.list.bz2 ${Res}/AncGenes/ancGenes.Euarchontoglires.list.bz2 -mode=drawKaryotype +sortBySize > ${Res}/Karyotype_Simu_${S1}_${S2}_tandemDupProp0.8.ps
#
## Assess Coevolution Score relevance
######################################
#Compute the coevolution scores of adjacent genes in a real genome and compare if coevolution scores are higher than if the genome was completely shuffled.
#: ${MagSimus}/analysis/quantifyCoevolutionScoreBetweenAdjacentGenes.py ${Res}/GenesST/genesST.Homo.sapiens.list.bz2 ${Res}/AncGenes/ancGenes.Euarchontoglires.list.bz2 ${Res}/geneTimeline.txt ${Res}/speciesTree.phylTree > ${Res}/CoevolutionScoreBetweenAdjacentGenes.svg
#
#Compute distances separating pairs of coevolving genes in a real genome and compare if distances are smaller than expected between two randomly chosen genes.
#: ${MagSimus}/analysis/quantifyDistanceBetweenClusterisedGenes.py ${Res}/GenesST/genesST.Homo.sapiens.list.bz2 ${Res}/AncGenes/ancGenes.Amniota.list.bz2 ${Res}/geneEvents.txt ${Res}/geneTimeline.txt ${Res}/speciesTree.phylTree -in:AllBlocks=${Res}/serialisedBlocks.pickle -minCoevolScore=0.2 > ${Res}/logClustering_MagSimus
#



Update DB:
DB=data78
for name in Mus.musculus Canis.lupus.familiaris Homo.sapiens Gallus.gallus Monodelphis.domestica ; do cp $gen/DB/genes/genesST.${name}.list.bz2 data/; done
for name in Amniota Theria Euarchontoglires Sauria Boreoeutheria; do cp $gen/DB/trees/0.25/ancGenes/all/ancGenes.${name}.list.bz2 data/; done
