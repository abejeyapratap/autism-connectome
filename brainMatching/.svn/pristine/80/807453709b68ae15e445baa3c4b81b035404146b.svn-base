/* 
 * File:   experiment.h
 * Author: yusuf
 *
 * Generated on September 21, 2016, 5:33 PM
 * Modified on January 21st, 2017 
 */

#ifndef EXPERIMENT_H
#define EXPERIMENT_H

#include <iostream> //cin,cout
#include <cstdlib> //atoi,stof, etc.
#include <string> //std::stof Note: stof requires C++11 as the compiler standard
#include <vector>
#include "dataset.h"
#include "edge.h"
#include "parameterSet.h"
#include "initializerParameters.h"


///this class includes the parameters for Experiments: several other classes are derived from this base class
//This class also includes some basic utility experiments such as calculating average connectomes etc.
class Experiment
{
    public:
		Experiment(InitializerParameters _initializerParameters, std::string _commandLine):initializerParameters(_initializerParameters), commandLine(_commandLine){}
		~Experiment(){};

		void listSubjectIDsOfSamplesUsedInExperiment(std::string dataFolder);

		////////utility functions
		//applies consistency thresholding for a given connectivity matrix
		static void thresholdConsistency(std::string controlsSubjectList,std::string allSubjectList, std::string streamlineFolder, float density, std::string outputFolder);
		//load group information for the subjects, for a 2 group dataset
		static void identifyGroups(std::string samplesFilename,std::string controlText, std::string patientText, std::vector<int> &controlsOrder, std::vector<int> &patientsOrder,std::vector<std::string> &controlIDs, std::vector<std::string> &patientIDs);

		//mutators for parameters
		void setFiles(std::string samplesFile, std::string outputPath, std::string normalizerFile);
		void setDataset(Dataset *_dataset1,Dataset *_dataset2=NULL);
		void setModality(std::string modality);
		void setPrint(bool _printMatches,bool _printSimilarityScores,bool _printAverageSimilarityScores, bool _printMatchingAccuracies, bool _printAverageMatchingAccuracies);
		
    protected:
		std::string samplesFile, groundTruthFile, outputPath, normalizerFile;
		Dataset *dataset,*dataset2;
		Edge::Feature edgeType1, edgeType2;
		std::string commandLine;//to keep record of the command line while generating a certain result file
		InitializerParameters initializerParameters;
		bool printMatches,printSimilarityScores,printAverageSimilarityScores,printMatchingAccuracies,printAverageMatchingAccuracies;
};

#endif /* EXPERIMENT_H */

