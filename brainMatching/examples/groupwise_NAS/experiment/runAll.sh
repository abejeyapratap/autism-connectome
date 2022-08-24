#!/bin/bash
# 

plotExtension=png #png or svg

###run experiments
# ./runMe.sh expRun

##process connectomes at connectome-level
# ./runMe.sh expProc

##process connectomes at Yeo system-level
# ./runMe.sh sysProc

## connectome level group difference (-- aka boxplot)
# ./runMe.sh grpDiff $plotExtension

## system level group difference (-- aka boxplot)
# ./runMe.sh sysDiff $plotExtension

## correlations between sub-groups
# need to provide correct NAS.res path
./runMe.sh correl



## linear mixed effects model for calculating change of NAS with days since injury -- (done to calc correl)
# ./runMe.sh mixedModel_time $plotExtension