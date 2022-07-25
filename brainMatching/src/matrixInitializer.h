/* 
 * File:   matrixInitializer.h
 * Author: yusuf
 *
 * Generated on September 19, 2016, 2:34 PM
 */

#ifndef MATRIXINITIALIZER_H
#define MATRIXINITIALIZER_H

#include <cstdlib>
#include <vector>
#include <map>
#include "graph.h"
#include "normalizer.h"
#include "initializerParameters.h"

class MatrixInitializer
{
	public:
		MatrixInitializer(InitializerParameters &_initializerParameters):initializerParameters(_initializerParameters){}
		void normalizeGraphs(std::map<int,Graph> &graphs){normalizer.normalizeData(graphs);}
		void initializeNormalizer(Normalizer &_normalizer){normalizer=_normalizer;}
		
		//initializers
		virtual void initializeAssignmentCostMatrix(float **assignmentCostMatrix, Graph &objects, Graph &labels)=0;
		virtual void initializeObjectsWeightMatrix(float **objectsWeightMatrix, Graph &objects)=0;
		virtual void initializeLabelsDistanceMatrix(float **labelsDistanceMatrix, Graph &labels)=0;
		
	protected:
		InitializerParameters initializerParameters;
		Normalizer normalizer;
};

#endif /* MATRIXINITIALIZER_H */

