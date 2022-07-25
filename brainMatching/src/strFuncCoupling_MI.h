/* 
 * File:   strFuncCoupling_MI.h
 * Author: yusuf
 *
 * Generated on February 22, 2018, 11:16 AM
 */

#ifndef STRFUNCCOUPLING_MI_H
#define STRFUNCCOUPLING_MI_H

#include "matrixInitializer.h"

class StrFuncCoupling_MI : public MatrixInitializer
{
	public:
		StrFuncCoupling_MI(InitializerParameters &_initializerParameters):MatrixInitializer(_initializerParameters){}
		
		//initializers
		void initializeAssignmentCostMatrix(float **assignmentCostMatrix, Graph &objects, Graph &labels) override;
		void initializeObjectsWeightMatrix(float **objectsWeightMatrix, Graph &objects) override;
		void initializeLabelsDistanceMatrix(float **labelsDistanceMatrix, Graph &labels) override;
};


#endif /* STRFUNCCOUPLING_MI_H */

