# -*- coding: utf-8 -*-
# Python 2.7 at least is needed, but not Python 3.*
#
#                   MagSimus v1.00
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

##################
# MagSimus v1.00 #
##################

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
: git clone https://github.com/DyogenIBENS/LibsDyogen.git
: echo "export PYTHONPATH=\"$PYTHONPATH:${INSTALLDIR}/LibsDyogen\"" >> ~/.bashrc

From now on we consider that the user is in the root directory of the project
: cd ${INSTALLDIR}/MagSimus

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

