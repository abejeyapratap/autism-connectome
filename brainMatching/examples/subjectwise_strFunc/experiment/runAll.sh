#!/bin/bash
# 

plotExtension=png #png or svg

###run experiments
./runMe.sh expRun  

###process connectomes
./runMe.sh expProc $plotExtension



