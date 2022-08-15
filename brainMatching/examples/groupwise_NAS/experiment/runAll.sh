#!/bin/bash
# 

plotExtension=png #png or svg

###run experiments
# ./runMe.sh expRun

##process connectomes
./runMe.sh expProc

## connectome level group difference (-- aka boxplot)
./runMe.sh grpDiff $plotExtension

## linear mixed effects model for calculating change of NAS with days since injury -- (done to calc correl)
# ./runMe.sh mixedModel_time $plotExtension