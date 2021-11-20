#!/bin/bash
#Launch all the commands in the README file and stops on errors if any
set -e
# Any subsequent commands which fail will cause the shell script to exit
# immediately
red='\e[0;31m'
green='\e[0;32m'
NC='\e[0m' # No Color


############################################
#	Check integrity of pre-processing scripts #
############################################
#preProcessCommandLines=(
## convet a .nhx tree into a protTree (forest of gene trees)
#"src/nhxGeneTrees2phylTreeGeneTrees.py data/geneTreesPrunedEnsembl.nhx > res/geneTreesPrunedEnsembl.protTree"
## convet a .nwk tree into a phylTree
#"src/newickSpeciesTree2phylTreeSpeciesTree.py data/speciesTree.nwk > res/speciesTree.phylTree"
## extract ancGenes (family)  from the species tree and the forest of gene trees
#"src/ancGenesFromGeneTrees.py res/speciesTree.phylTree res/geneTreesPrunedEnsembl.protTree -out:ancGenes=res/ancGenes.example.%s.list.bz2 > res/geneTreesPrunedEnsembl.afterExtractingAncGenes.protTree"
#)
#for line in "${preProcessCommandLines[@]}"
#	do
#		echo -e "${green}${line}${NC}"
#		eval ${line}
#done

#############################################
##	Check integrity of magSimus1  					 #
#############################################
magSimus1CommandLines=(
"src/magSimus1.py data/speciesTree.phylTree -out:genomeFiles=res/simu1/genes.%s.list.bz2 -out:ancGenesFiles=res/simu1/ancGenes.%s.list.bz2 -parameterFile=data/parametersG_ABC.v83 -userRatesFile=data/specRates_MS1.v84 > res/simu1/logStd"
)
for line in "${magSimus1CommandLines[@]}"
	do
		echo -e "${green}${line}${NC}"
		eval ${line}
done

########################################################
##	Check the quality of PhylDiag sbs using magSimus 1 #
########################################################
phylDiagUsingMagSimus1CommandLines=(
"src/magSimus1.py data/speciesTree.phylTree -out:genomeFiles=res/simu1/genes.%s.list.bz2 -out:ancGenesFiles=res/simu1/ancGenes.%s.list.bz2 -userRatesFile=data/specRates_MS1.v84 +breakpointAnalyzer -out:empiricalSbsAsGenomes=res/simu1/sbs.genes.%s.%s.list.bz2 > res/simu1/logStdout"
"src/analysis/comparePhylDiagSbsToSimulatedSbs.py data/speciesTree.phylTree Mus.musculus Gallus.gallus Amniota"
)
for line in "${phylDiagUsingMagSimus1CommandLines[@]}"
	do
		echo -e "${green}${line}${NC}"
		eval ${line}
done

#############################################
#	Check integrity of magSimus2  					  #
#############################################
magSimus2CommandLines=(
"src/calcGeneEvolution.py data/speciesTree.phylTree -genesFile=data/genesST.%s.list.bz2 -ancGenesFile=data/ancGenes.%s.list.bz2 > res/simu2/geneEvents.txt"
"src/getGeneTimeline.py data/speciesTree.phylTree res/simu2/geneEvents.txt > res/simu2/geneTimeline.txt"
"src/magSimus2.py data/speciesTree.phylTree res/simu2/geneEvents.txt  res/simu2/geneTimeline.txt -userRatesFile=data/specRates_MS2.v80 -out:genomeFiles=res/simu2/genes.%s.list.bz2 -out:ancGenesFiles=res/simu2/ancGenes.%s.list.bz2 -b_randomAccel > res/simu2/logStd"
"src/magSimus2.py data/speciesTree.phylTree res/simu2/geneEvents.txt res/simu2/geneTimeline.txt -userRatesFile=data/specRates_MS2.v80 -out:genomeFiles=res/simu2/genes.%s.list.bz2 -out:ancGenesFiles=res/simu2/ancGenes.%s.list.bz2 -in:firstAncestralGenome=res/simu2/genes.Amniota.list.bz2 -b_randomAccel > res/simu2/logStd_withGiven1stAncestralGenome"
"src/magSimus2.py data/speciesTree.phylTree res/simu2/geneEvents.txt res/simu2/geneTimeline.txt -userRatesFile=data/specRates_MS2.v80 -out:genomeFiles=res/simu2/genes.%s.list.bz2 -out:ancGenesFiles=res/simu2/ancGenes.%s.list.bz2 -parameterFile=data/parameters.v80 > res/simu2/logStd"
)
for line in "${magSimus2CommandLines[@]}"
	do
		echo -e "${green}${line}${NC}"
		eval ${line}
done

#################################################
# Check integrity of analysis/scripts  				  #
#################################################
W4=/kingdoms/dyogen/workspace4/workspace4/jlucas
INSTALLDIR=${W4}/Libs
MagSimus=${INSTALLDIR}/MagSimus
Res=${W4}/Workflow/MagSimus/Data72/100SimusSpecRates74
#############################################
# Distribution of chromosome lengths
#############################################
S=Homo.sapiens

analysisCommandLines=(
"${MagSimus}/src/analysis/distribChrLengths_Real_vs_Simus.py ${Res}/GenesST/genesST.${S}.list.bz2 -in:simulatedGenomes=${Res}/%s/genes.${S}.list.bz2 > ${Res}/distribChrLengths_${S}.svg"
)
#############################################
tandemGapMax=2
gapMax=5
DM=CD
filterType=InCommonAncestor
if [ "${filterType}" == "InCommonAncestor" ]
then
	fM="Ica"
elif [ "${filterType}" == "InBothSpecies" ]
then
	fM="IbS"
elif
	[ "${filterType}" == "None" ]
then
	fM="None"
	echo 'ERROR: The filterType is unknown. It should be either "InCommonAncestor", "InBothSpecies" or "None"' >&2
fi

#############################################
# Parallel execution of the extractions of sbs for real genomes and simulated
# genomes.
#############################################

# FIXME Launch by hand the two next lines
# For real genomes
#  ${MagSimus}/src/analysis/extractSbs_S1_vs_S2.bsh -f InCommonAncestor -t ${tandemGapMax} -g ${gapMax} -d ${DM} -r real ${Res}/GenesST ${Res}/AncGenes ${MagSimus}/src/analysis/speciesCombi.txt ${Res}
# For simulated genomes
#  ${MagSimus}/src/analysis/extractSbs_S1_vs_S2.bsh -f InCommonAncestor -t ${tandemGapMax} -g ${gapMax} -d -r simu ${DM} ${Res} ${Res}/AncGenes ${MagSimus}/src/analysis/speciesCombi.txt ${Res}
# TODO monitor the condor execution
# wait for end of computations to launch the next computations

###############################################
# N50 of pairwise comparisons
###############################################
S1="%s" && S2="%s" && CA="%s"
sbsFileName=${S1}-${S2}_${fM}${CA}_Tgm${tandemGapMax}Gm${gapMax}${DM}
analysisCommandLines+=(
"${MagSimus}/src/analysis/graph_N50_Real_MagSimus.py ${MagSimus}/src/analysis/speciesCombi.txt -realDiags=${Res}/real_${sbsFileName}.sbs -simuDiags=${Res}/%s/simu_${sbsFileName}.sbs > ${Res}/N50_pairwiseComp_fLCA_Tgm${tandemGapMax}Gm${gapMax}${DM}.svg"
)
###############################################
# Distribution of the lengths of synteny blocks for a pairwise comparison of
# two species
###############################################
S1=Homo.sapiens && S2=Mus.musculus && CA=Euarchontoglires
sbsFileName=${S1}-${S2}_${fM}${CA}_Tgm${tandemGapMax}Gm${gapMax}${DM}

analysisCommandLines+=(
"${MagSimus}/src/analysis/distribSbsLengths_Real_vs_Simu.py ${Res}/real_${sbsFileName}.sbs -in:Sim%s_S1_vs_S2.sbs=${Res}/%s/simu_${sbsFileName}.sbs > ${Res}/distribSbsLengths_${sbsFileName}.svg"
)
###############################################
# Dot Matrices of the comparison of two genomes
###############################################
S1=Homo.sapiens && S2=Mus.musculus && CA=Euarchontoglires
dotMatrixName="DotMatrix_${S1}-${S2}_f${CA}"

analysisCommandLines+=(
"${MagSimus}/src/analysis/misc.compareGenomes.py ${Res}/GenesST/genesST.${S1}.list.bz2 ${Res}/GenesST/genesST.${S2}.list.bz2 ${Res}/AncGenes/ancGenes.${CA}.list.bz2 +sortBySize -matrix:pointSize=0.05 > ${Res}/real${dotMatrixName}.ps"
)

simu=0

analysisCommandLines+=(
"${MagSimus}/src/analysis/misc.compareGenomes.py ${Res}/${simu}/genes.${S1}.list.bz2 ${Res}/${simu}/genes.${S2}.list.bz2  ${Res}/AncGenes/ancGenes.${CA}.list.bz2 +sortBySize -matrix:pointSize=0.05 > ${Res}/simu${simu}_${dotMatrixName}.ps"
)
###############################################
# Karyotype of the comparison of two genomes
###############################################
S1=Homo.sapiens && S2=Mus.musculus && CA=Euarchontoglires
karyotypeName=Karyotype_${S1}-${S2}_f${CA}

analysisCommandLines+=(
"${MagSimus}/src/analysis/misc.compareGenomes.py ${Res}/GenesST/genesST.${S1}.list.bz2 ${Res}/GenesST/genesST.${S2}.list.bz2 ${Res}/AncGenes/ancGenes.${CA}.list.bz2 -mode=drawKaryotype +sortBySize > ${Res}/real${karyotypeName}.ps"
)

simu=0

analysisCommandLines+=(
"${MagSimus}/src/analysis/misc.compareGenomes.py ${Res}/${simu}/genes.${S1}.list.bz2 ${Res}/${simu}/genes.${S2}.list.bz2 ${Res}/AncGenes/ancGenes.${CA}.list.bz2 -mode=drawKaryotype +sortBySize > ${Res}/simu${simu}${karyotypeName}.ps"
)
###############################################
# Convert from .ps to .pdf
###############################################
WD=$(pwd)
cd ${Res}
for i in `ls *.ps`; do filename=$(basename "${i}"); extension="${filename##*.}"; filename="${filename%.*}"; ps2pdf -dAutoRotatePages=/None ${filename}.ps ${filename}.pdf; done
cd ${WD}
###############################################
# Asses the coevolution score relevance
###############################################
analysisCommandLines+=(
"${MagSimus}/src/analysis/quantifyCoevolutionScoreBetweenAdjacentGenes.py ${Res}/GenesST/genesST.${S1}.list.bz2 ${Res}/AncGenes/ancGenes.${CA}.list.bz2 ${Res}/geneTimeline.txt ${Res}/speciesTree.phylTree > ${Res}/CoevolutionScoreBetweenAdjacentGenes.svg"
)

analysisCommandLines+=(
"${MagSimus}/src/analysis/quantifyDistanceBetweenClusterisedGenes.py ${Res}/GenesST/genesST.${S1}.list.bz2 ${Res}/AncGenes/ancGenes.${CA}.list.bz2 ${Res}/geneEvents.txt ${Res}/geneTimeline.txt ${Res}/speciesTree.phylTree -in:AllBlocks=${Res}/serialisedBlocks.pickle -minCoevolScore=0.2 > ${Res}/logClustering_MagSimus"
)
###############################################
# Tandem duplicates analysis
###############################################
S=Homo.sapiens
CA=Euarchontoglires
tandemGapMax=5
analysisCommandLines+=(
"${MagSimus}/src/analysis/countTandemDuplicates.py ${MagSimus}/data/genesST.${S}.list.bz2 ${MagSimus}/data/ancGenes.${CA}.list.bz2 -tandemGapMax=${tandemGapMax}"
)

for line in "${analysisCommandLines[@]}"
	do
		echo -e "${green}${line}${NC}"
		eval ${line}
done
