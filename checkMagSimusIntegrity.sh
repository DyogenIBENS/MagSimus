#!/bin/bash
#Launch all the commands in the README file and stops on errors if any
set -e
set -x
# Any subsequent commands which fail will cause the shell script to exit
# immediately
red='\e[0;31m'
green='\e[0;32m'
NC='\e[0m' # No Color

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
