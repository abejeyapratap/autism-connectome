/* 
 * File:   subjectwiseExperiment.h
 * Author: yusuf
 *
 * Generated on May 1, 2018, 2:00 PM
 */

#ifndef SUBJECTWISEEXPERIMENT_H
#define SUBJECTWISEEXPERIMENT_H

#include <iostream> //cin,cout
#include <cstdlib> //atoi,stof, etc.
#include <string> //std::stof Note: stof requires C++11 as the compiler standard
#include <vector>
#include "experiment.h"


//this class includes experiments for structure-function coupling for each subject of a dataset. 
//Specifically, str-func, str-str or func-func of each subject is matched. No inter-subject matching is performed in this class.
class SubjectwiseExperiment : public Experiment
{
    public:
        SubjectwiseExperiment(InitializerParameters _initializerParameters, std::string _commandLine):Experiment(_initializerParameters,_commandLine){}
        ~SubjectwiseExperiment(){};

        void subjectwiseMatchingExperiment(int numPermutation=1, bool saveAvgConnectome=false, int seedSupplement=0, std::string shuffleType="",int shuffleDataset=1);
        void subjectwiseDistanceExperiment(std::string distanceMeasure);

        void saveStructureFunctionCorrelationOfSingleSubject(int subjectId);
        void saveAverageStructureFunctionCorrelationOfAllSubjects();

    private:
        void saveSubjectwiseMatchingResults(int numNodes, int numSubjects, int **matchesMatrix, std::vector<float>& similarityScores, std::string outfile="");
};

#endif /* SUBJECTWISEEXPERIMENT_H */

