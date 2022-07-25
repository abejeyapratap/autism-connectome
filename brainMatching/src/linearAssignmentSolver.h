/* 
 * File:   linearAssignmentSolver.h
 * Author: yusuf
 *
 * Generated on July 7, 2017, 11:05 AM
 */

#ifndef LINEARASSIGNMENTSOLVER_H
#define LINEARASSIGNMENTSOLVER_H

#include <cstdlib>
#include <math.h>
#include <iostream>
#include <fstream>
#include <map>
#include "graph.h"
#include "utility.h"
#include "matrixInitializer.h"

class LinearAssignmentSolver
{
	public:
		LinearAssignmentSolver();
		LinearAssignmentSolver(const Graph & objectGraph, const Graph & labelGraph);
		virtual ~LinearAssignmentSolver();

		float doLabeling(MatrixInitializer &initializer);
		void  incrementMatchesMatrix(int **matchesMatrix);
		void  saveMatchesToArray(int *matchesMatrix);
		void  saveMatchesToVector(std::vector<int> &matchesVector);

		inline std::map<int,int>::iterator getMatchesIteratorBegin(){return matches.begin();}
		inline std::map<int,int>::iterator getMatchesIteratorEnd(){return matches.end();}

		//functions for DEBUG
		void printMatches();
		void printObjectsWeightMatrix();
		void printAssignmentCostMatrix();
		void printLabelsDistanceMatrix();
		float  calculateIdentityMatchesRatio();

	protected:
	   virtual float solveOptimizationProblem(double ** pairingProbabilities);
	   float calculateObjectiveValue(std::map<int, int> &_matches);

	   float **objectsWeightMatrix;
	   float **labelsDistanceMatrix;
	   float **assignmentCostMatrix;

	   std::map<int, int> matches;

	   Graph objects;
	   Graph labels;
};

#endif /* LINEARASSIGNMENTSOLVER_H */

