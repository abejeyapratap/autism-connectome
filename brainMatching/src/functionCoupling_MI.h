/* 
 * File:   functionCoupling_MI.h
 * Author: yusuf
 *
 * Generated on February 22, 2018, 11:30 AM
 */

#ifndef FUNCTIONCOUPLING_MI_H
#define FUNCTIONCOUPLING_MI_H

#include "matrixInitializer.h"

class FunctionCoupling_MI : public MatrixInitializer
{
	public:
		FunctionCoupling_MI(InitializerParameters &_initializerParameters):MatrixInitializer(_initializerParameters){}
		
		//initializers
		void initializeAssignmentCostMatrix(float **assignmentCostMatrix, Graph &objects, Graph &labels) override;
		void initializeObjectsWeightMatrix(float **objectsWeightMatrix, Graph &objects) override;
		void initializeLabelsDistanceMatrix(float **labelsDistanceMatrix, Graph &labels) override;
};

#endif /* FUNCTIONCOUPLING_MI_H */

