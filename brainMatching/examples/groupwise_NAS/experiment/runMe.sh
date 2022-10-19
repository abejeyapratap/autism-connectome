#!/bin/bash
# 

job=${1} #expRun, expProc, grpDiff, grpDiffLong, corr, complete (to do everything)
plotExtension=${2} #png or svg

experimentFolder=.

### path to codes that do brain match, python3 code for processing matching results and calculating group differences
brainMatch=$experimentFolder/../../../bin/brainMatch
scriptFolderPath=$experimentFolder/../scripts
processResults_py=$scriptFolderPath/processResults.py
groupDifference_similarity_py=$scriptFolderPath/groupDifference_similarity.py
listScores_py=$scriptFolderPath/listScores.py
mixedModel_similarity_time_R=$scriptFolderPath/mixedModel_similarity_time.R
mixedModel_plot_similarity_time_py=$scriptFolderPath/mixedModel_plot_similarity_time.py


### parameters to the matching algorithm: assignmentCost can be one of the following <edgesIncludeDiagZeroDiag,edgesIgnoreDiag, edgesIncludeDiagRandDiag> 
### with the hardness of the problem to be solved in increasing order
assignmentCost=" -assCost edgesIgnoreDiag"
pathType="-pathType direct" # could have been "-pathType shortestPath" or "-pathType wCommunicability -pathLength 2"

# preprocessGraphs="-preprocessGraphs none"
preprocessGraphs="-preprocessGraphs normalizeEdges"
# preprocessGraphs="-preprocessGraphs logScaleEdgesStructure"
# preprocessGraphs="-preprocessGraphs logScaleEdgesStructure_normalizeEdges"

connectomeName="schaefer" # for correlation

### path to connectomes, list of samples to be used in the experiment, and cognitive scores of samples
# subjectsInfoPath='../data/tbi_longitudinal_dtiQAPass_20210121.csv'
# samples=$experimentFolder/../data/subjects_qa.txt
# connectomes=$experimentFolder/../data/connectomes_desikan/norm_desikan_2
# samples=$experimentFolder/../data/subjectNames/desikan_filtered.txt
# samples=$experimentFolder/../data/smallD.txt
# samples=$experimentFolder/../data/desikan_male.txt
connectomes=$experimentFolder/../data/connectomes_schaefer/schaefer_200 ## 200 ROIs instead of 220
# samples=$experimentFolder/../data/smallS.txt
# samples=$experimentFolder/../data/schaefer_below12.txt
# samples=$experimentFolder/../data/subjectNames/schaefer_male.txt
samples=$experimentFolder/../data/subjectNames/schaefer_filtered.txt
# connectomes=$experimentFolder/../data/connectomes_schaefer/norm_schaefer_200_100 # MY manual normalization
# connectomes=$experimentFolder/../data/connectomes_schaefer
subjectsInfoPath='../data/tobacco_demographics.csv'

### output path for the results and plots
# results=$experimentFolder/results
results=$experimentFolder/results_schaefer/norm_2_schaefer220
plotsRoot=$experimentFolder/plots
numSys=7

########### run matching experiment ##############
##calculate distance between healthy controls using graph matching
if( [ "$job" == "expRun" ] || [ "$job" == "complete" ] );then
	echo -e "\tRunning groupwise brain matching experiment..."
	mkdir -p $results
	$brainMatch -experiment groupwise match -groups all all -data matrix 1 -dti $connectomes/ -samples $samples -printMatches -printSimilarity -modality str_str -outputPath $results/matching_raw $pathType $assignmentCost $preprocessGraphs
fi

########### system & connectome level analysis ##############
sysFileOut=$results/sys_level
sysResults_py=$scriptFolderPath/processSystems.py
sysYeo=$experimentFolder/../data/yeo_7systems_schaefer.txt
if( [ "$job" == "sysProc" ] || [ "$job" == "complete" ] );then
	echo -e "\tSystem-level processing raw experiment results..."
	echo -e "\tNumber of sub-systems: $numSys"
	mkdir -p $sysFileOut

	python3 $sysResults_py -r $results/"matching_raw_0.res" -o $sysFileOut -s $samples -mt accuracy -al group --relativeTo healthy --sysMap $sysYeo --numSys $numSys
fi

resultFileAccuracy=$results/NAS_rt_healthy.res
if( [ "$job" == "expProc" ] || [ "$job" == "complete" ] );then
	echo -e "\tProcessing raw experiment results..."

	python3 $processResults_py -r $results/"matching_raw_0.res" -o $resultFileAccuracy -s $samples -mt accuracy -al group --relativeTo healthy
fi

########### group difference calcs ##############
sysDifference_similarity_py=$scriptFolderPath/groupDiff_system.py
if( [ "$job" == "sysDiff" ] || [ "$job" == "complete" ] );then
	echo -e "\tSystem-level group difference and generating "$plotType" plots..."
	echo -e "\tNumber of sub-systems: $numSys"
	outpath=$results/sys_level/plots
	mkdir -p $outpath

	plotType=box #box or violin

	python3 $sysDifference_similarity_py -r $resultFileAccuracy -o $outpath/ -s $samples --scoreType dist  --plotExtension $plotExtension --plotType $plotType --timePoints any --sysFiles $sysFileOut --numSys $numSys
	python3 $sysDifference_similarity_py -r $resultFileAccuracy -o $outpath/ -s $samples --scoreType dist  --plotExtension $plotExtension --plotType $plotType --timePoints any --noHealthyPlot --sysFiles $sysFileOut --numSys $numSys
fi

if( [ "$job" == "grpDiff" ] || [ "$job" == "complete" ] );then
	echo -e "\tCalculating group difference and generating "$plotType" plots..."
	outpath=$plotsRoot/groupDifference
	mkdir -p $outpath/histograms

	plotType=box #box or violin

	python3 $groupDifference_similarity_py -r $resultFileAccuracy -o $outpath/ -s $samples --scoreType dist  --plotExtension $plotExtension --plotType $plotType --timePoints any
	python3 $groupDifference_similarity_py -r $resultFileAccuracy -o $outpath/ -s $samples --scoreType dist  --plotExtension $plotExtension --plotType $plotType --timePoints any --noHealthyPlot
fi

correl_py=$scriptFolderPath/correlation.py
if( [ "$job" == "correl" ] );then
	echo -e "\tCalculating correlations at connectome-level..."
	correlOut=$plotsRoot/correl
	mkdir -p $correlOut
	correlConnectome=$correlOut/$connectomeName
	mkdir -p $correlConnectome

	python3 $correl_py --subjectsInfoPath $subjectsInfoPath -r $resultFileAccuracy -o $correlConnectome --ados
	python3 $correl_py --subjectsInfoPath $subjectsInfoPath -r $resultFileAccuracy -o $correlConnectome --no-ados
fi

correlSys_py=$scriptFolderPath/correlationSystem.py
if( [ "$job" == "correlSys" ] );then
	echo -e "\tCalculating correlations at system-level..."
	correlOut=$plotsRoot/correl
	correlConnectome=$correlOut/$connectomeName/sys_level
	mkdir -p $correlConnectome
	resPaths=$results/sys_level

	python3 $correlSys_py --subjectsInfoPath $subjectsInfoPath -r $resPaths -o $correlConnectome
fi


##### Linear model ##### 
if( [ "$job" == "mixedModel_time" ] || [ "$job" == "complete" ] );then
	echo -e "\tCalculating mixed effects model for predicting similarity score from DSI..."
	#generate scores file where cognitive socres and matching accuracy scores are listed for each subject
	mkdir -p $results
	python3 $listScores_py -r $resultFileAccuracy -sl $samples  --subjectsInfoPath $subjectsInfoPath -tOut $results/ --scoreType dist --timePoints any

	#run mixed effects model for proSpeed, executive, and verbal
	#and then generate plots for the model
	resultsPath=$results/mixedModel
	plotPath=$plotsRoot/mixedModel
	mkdir -p $plotPath $resultsPath

	### mixed effect model to predict similarity from days since injury: similarityScore ~ DSI + age + (1|subjectID)
	echo -e "\t\tCalculating mixed effects model: similarityScore ~ DSI + age + (1|subjectID)"
	for tp in none; do # <none> using all time points, <s2> removing second time point
		Rscript $mixedModel_similarity_time_R  --scoreFile $results/allScores.txt --summaryOutputFile $resultsPath/"similarityTime_model_summary_"$tp"Discarded.txt" --interceptOutputFile $resultsPath/"similarityTime_intercept_"$tp"Discarded_" --discardTimePoint $tp > $resultsPath/"similarityTime_model_details_"$tp"Discarded.txt" 2> $resultsPath/"errors_similarityTime_"$tp"Discarded.txt"

		coloringScheme="--coloringScheme group --grayLines mixedModel --lineCoef $resultsPath/similarityTime_intercept_" #alternatively --coloringScheme <gender,group> --dimColorWithAge --noSubjectName --grayLines <mixedModel,spaghetti>

		python3 $mixedModel_plot_similarity_time_py -m $resultsPath/"similarityTime_model_summary_" -s $results/allScores.txt -pOut $plotPath/ --plotExtension $plotExtension --discardTimePoint $tp  --coloringScheme group --grayLines spaghetti --noSubjectName --lineCoef $resultsPath/similarityScore_intercept_"$tp"Discarded_similarityScore0.txt 
	done
fi

