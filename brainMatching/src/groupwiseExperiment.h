/* 
 * File:   groupwiseExperiment.h
 * Author: yusuf
 *
 * Generated on February 22, 2018, 3:55 PM
 */

#ifndef GROUPWISEEXPERIMENT_H
#define GROUPWISEEXPERIMENT_H


#include <iostream> //cin,cout
#include <cstdlib> //atoi,stof, etc.
#include <string> //std::stof Note: stof requires C++11 as the compiler standard
#include <vector>
#include "experiment.h"

//this class includes experiments where two groups of subjects are compared with each other for a single modality.
//That is, each subject in the first group is compared with each subject of the second group.
class GroupwiseExperiment : public Experiment
{
    public:
		GroupwiseExperiment(InitializerParameters _initializerParameters, std::string _commandLine):Experiment(_initializerParameters,_commandLine){}
		~GroupwiseExperiment(){};

		void groupwiseMatchingExperiment(int numPermutation=1,std::string controlsLabel="all",std::string patientsLabel="all");
			void groupwiseCorrelationExperiment(std::string controlsLabel="all",std::string patientsLabel="all");
		void groupwiseMatrixDistanceExperiment(std::string distanceMeasure,std::string controlsLabel="all",std::string patientsLabel="all");
		void groupwiseSubjectSpecificMatrixSubtractionExperiment(std::string controlsLabel="all",std::string patientsLabel="all");
		void groupwiseSubjectSpecificMatrixZScoreExperiment(std::string controlsLabel="all",std::string patientsLabel="all");

		void calculateAverageResults(std::string resultsFolder,std::string controlText, std::string patientText);

    private:
		void saveMatchingResults(int numOfGroup1, int numOfGroup2, int numOfNodes, int ***matchesMatrix, std::vector<float> *averageSimilarityScores, float **similarityScores, std::vector<float> *averageMatchingAccuracies, float **matchingAccuracies, std::string outfile, bool printInformationLines=true);
        void saveSimilarityResults(int numOfGroup1, int numOfGroup2, int numOfNodes, std::vector<float> *averageSimilarities, float **similarityScores, std::string outfile, bool printInformationLines=true);
};


#endif /* GROUPWISEEXPERIMENT_H */

