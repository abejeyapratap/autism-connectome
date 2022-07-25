#!/bin/bash

#$ -S /bin/bash
#$ -cwd

#$ -o /cbica/home/osmanliy/sge_job_output/$JOB_ID.stdout
#$ -e /cbica/home/osmanliy/sge_job_output/$JOB_ID.stderr


# Check if the lab-specific variable is defined...if not, fall back to "/tmp"
if [ ! -z ${SBIA_TMPDIR} ] ; then
    tmpDataDIR=${SBIA_TMPDIR}

    # If there's no subdirectory named for the current user, try to make one
    if [ ! -d ${tmpDataDIR}/${USER} ] ; then
        mkdir ${tmpDataDIR}/${USER}
        if [ $? != 0 ] ; then
            echo "Failed to make directory ${tmpDataDIR}/${USER}"
        fi
    fi

    # use mktemp to create a unique random subdirectory, with a
    # name based on the programname
    tmpDataDIR=`mktemp -d ${tmpDataDIR}/${USER}/programname.XXXXXX`
    if [ $? != 0 ] ; then
        echo "Failed to make a temp directory in ${tmpDataDIR}/${USER}"
        exit
    fi
else
    tmpDataDIR=${outDataDIR}/tmp
fi

#parameters coming from the runMe.sh

pathType=${1} #('direct' 'shortest' 'strongest' 'wCommunicability' 'uCommunicability' 'pathTransitivity' 'searchInformation') #${8}
pathLength=${2} #used for wCommunicability and uCommunicability
preprocessGraphs=${3} #traffic_logScaleEdgesStructure_normalizeEdges #traffic_logScaleEdgesStructure_normalizeEdges
assignmentCost=${4} #edgesIgnoreDiag, edgesIncludeDiagRandDiag, edgesIncludeDiagZeroDiag
atlas=${5} #Lausanne234
functionType=${6} #full,lasso
funcConn=${7} #'positive' or 'negative'
taskPerJob=${8} #how many tasks per job (needed for matching if randDiag is chosen, and needed for permutation test for shuffling graphs)
seed=${9} #seed to start permutation test
task=${10} #matching or permutation

solver=LinAss #('Diagonal' 'MLPD' 'QAPD' 'LinAss') #let this be fixed for this study
#beta=${4} #this parameter will be used if the solver was QAPD

#exact paths for the input files
codePath=/cbica/home/osmanliy/code/tbiStrAndFunc/dist/Debug_Special/GNU-Linux-x86
connectomes_str=/cbica/home/osmanliy/data/PNC/dti/${atlas}
connectomes_func=/cbica/home/osmanliy/data/PNC/fmri/${atlas}"_"${functionType}
samples=/cbica/home/osmanliy/data/PNC/src/samples/list_jacob_fmri.txt #allSamples_reduced.txt

strFuncCouplingResultsFolder=/cbica/home/osmanliy/comp_space/$task
outputFileFolder=$atlas"__"$functionType"__"$pathType"__"$funcConn"__"$preprocessGraphs"__"$assignmentCost
outputFileName=$pathLength"__"$seed"__"
outputFilePath=$strFuncCouplingResultsFolder/$outputFileFolder/$outputFileName

module load gcc/5.2.0

#if the outfile folder was not generated before, generate it now
mkdir -p $strFuncCouplingResultsFolder/$outputFileFolder

#to run for a single task: #./brainMatch -experiment strFuncCoupling -permutation 100 -graphs ~/data/PNC/graphs/structureFunctionMappingProject/notNormalized/ -samples ~/data/PNC/src/samples/allSamples.txt -printMatches -outFile QAPD_matching -solver QAPD -beta 0.5
if( [ "$task" == "matching" ] ); then
	echo "${codePath}/brainMatch -experiment subjectwise match -permutation $taskPerJob -modality str_func -funcConn $funcConn -data matrix 1 -dti $connectomes_str/ -fmri $connectomes_func/ -samples $samples -printMatches -outputPath $outputFilePath -solver $solver -pathType $pathType -pathLength $pathLength -assCost $assignmentCost -preprocessGraphs $preprocessGraphs"
	
	${codePath}/brainMatch -experiment subjectwise match -permutation $taskPerJob -modality str_func -funcConn $funcConn -data matrix 1 -dti $connectomes_str/ -fmri $connectomes_func/ -samples $samples -printMatches -outputPath $outputFilePath -solver $solver -pathType $pathType -pathLength $pathLength -assCost $assignmentCost -preprocessGraphs $preprocessGraphs
fi

if( [ "$task" == "permutation" ] ); then
	echo "${codePath}/brainMatch -experiment subjectwise match -permutation $taskPerJob -shuffle structure 1 -seed $seed -modality str_func -funcConn $funcConn -data matrix 1 -dti $connectomes_str/ -fmri $connectomes_func/ -samples $samples -printMatches -outputPath $outputFilePath -solver $solver -pathType $pathType -pathLength $pathLength -assCost $assignmentCost -preprocessGraphs $preprocessGraphs"
	
	${codePath}/brainMatch -experiment subjectwise match -permutation $taskPerJob -shuffle structure 1 -seed $seed -modality str_func -funcConn $funcConn -data matrix 1 -dti $connectomes_str/ -fmri $connectomes_func/ -samples $samples -printMatches -outputPath $outputFilePath -solver $solver -pathType $pathType -pathLength $pathLength -assCost $assignmentCost -preprocessGraphs $preprocessGraphs
fi
