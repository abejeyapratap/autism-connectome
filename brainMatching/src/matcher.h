/* 
 * File:   matcher.h
 * Author: yusuf
 *
 * Generated on September 15, 2016, 4:29 PM
 */

#ifndef MATCHER_H
#define MATCHER_H

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <algorithm> // std::find
#include <dirent.h>  //for readdir, opendir, closedir
#include "dataset.h"
#include "graph.h"
#include "edge.h"
#include "normalizer.h"
#include "parameterSet.h"
#include "initializerParameters.h"
#include "matrixInitializer.h"

class Matcher
{
	public:
		Matcher(const Matcher &other){dataset= other.dataset;dataset2=other.dataset2;};
		Matcher(Dataset *_dataset,Dataset *_dataset2=NULL){dataset= _dataset;dataset2=_dataset2;};
		~Matcher(){}
		
		void matchGroups(float** similarityScores, int*** matchesMatrix, std::vector<int> &group1IDs, std::vector<int> &group2IDs, Edge::Feature edgeType1, Edge::Feature edgeType2);
		void matchEveryoneWithItself(std::vector<float> &similarityScores, int **matchesMatrix, Edge::Feature edgeType1, Edge::Feature edgeType2);
      
		//legacy code for the grail paper, works with leaveOneOutExperiment.cpp
		void matchEveryoneWithEveryone(std::string outputFile, int rowStart, int rowEnd, int columnStart, int columnEnd);
		void matchEveryoneWithEveryone(std::string outputFile, bool saveMatches=false);
		void matchEveryoneWithEveryone(float **resultMatrix, int ***matchesMatrix);
		void matchOneWithEveryone(int testGraph,float **resultMatrix, int **matchesMatrix=NULL);
		void matchTwoGraphs(int g1, int g2, bool printMatches);

	private:     
		Dataset *dataset,*dataset2;
      
		struct TableRow
		{
			int obj1, obj2;
			int timeInit;
			int timeRun;
			float objValue;
			std::vector<int> matches;
		};
		
		void saveResults(std::string outputFile, int taskSize, TableRow *results, bool saveMatches=false);
		MatrixInitializer *allocateMatrixInitializer(Edge::Feature edgeType1, Edge::Feature edgeType2);
};

#endif /* MATCHER_H */

