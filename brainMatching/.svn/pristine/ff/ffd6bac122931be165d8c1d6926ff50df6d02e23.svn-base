#!/bin/bash
#usage: ./runMeExperiments.sh {samplesFilename} {taskPerJob} {taskCount} {beta} {funcConn} {solver} {pathType} {pathLength} {preprocessGraphs} {assignmentCost}
#sample: ./runMeExperiments.sh allSamples_clean.txt 10 1000 0.5 positive LinAss direct 0 logScaleEdgesStructure_traffic edgesIncludeDiagRandDiag
#for sign in positive negative; do for path in 'direct' 'wShortest' 'uShortest' 'searchInformation' 'pathTransitivity' 'uCommunicability' 'wCommunicability';do ./runStrFuncCouplingExperiment_PNC.sh allSamples_Lausanne129_full.txt 200 999 0 $sign LinAss $path 2 logScaleStructuralConnectome_traffic_normalizeEdges edgesIncludeDiagRandDiag streamline_full Lausanne129; done; done


pathType=${1} #('direct' 'shortest' 'strongest' 'wCommunicability' 'uCommunicability' 'pathTransitivity' 'searchInformation') #${8}
pathLength=${2} #used for wCommunicability and uCommunicability
preprocessGraphs=${3} #traffic_logScaleEdgesStructure_normalizeEdges #traffic_logScaleEdgesStructure_normalizeEdges
assignmentCost=${4} #edgesIncludeDiagZeroDiag #edgesIgnoreDiag, edgesIncludeDiagRandDiag, edgesIncludeDiagZeroDiag
atlas=${5} #Lausanne234
functionType=${6} #full,lasso
funcConn=${7} #'positive' or 'negative'
taskPerJob=${8} #how many tasks per job (needed for matching if randDiag is chosen, and needed for permutation test for shuffling graphs)
taskCount=${9} #total task count
memory_match=${10} #2G
memory_permutation=${11} #2G

currTime=$(date +%s)
mkdir -p jobOutputs_$currTime
for seed in $(seq 0 $taskPerJob $taskCount); do
	echo "qsub -l h_vmem=$memory_match qsub_nimg.sh " $pathType $pathLength $preprocessGraphs $assignmentCost $atlas $functionType $funcConn $taskPerJob $seed matching >> jobOutputs_$currTime/matching_"$pathLength".jobs
   	qsub -l h_vmem=$memory_match qsub_nimg.sh $pathType $pathLength $preprocessGraphs $assignmentCost $atlas $functionType $funcConn $taskPerJob $seed matching
	
	echo "qsub -l h_vmem=$memory_match qsub_nimg.sh " $pathType $pathLength $preprocessGraphs $assignmentCost $atlas $functionType $funcConn $taskPerJob $seed permutation >> jobOutputs_$currTime/permutation_"$pathLength".jobs
   	qsub -l h_vmem=$memory_permutation qsub_nimg.sh $pathType $pathLength $preprocessGraphs $assignmentCost $atlas $functionType $funcConn $taskPerJob $seed permutation
done 
