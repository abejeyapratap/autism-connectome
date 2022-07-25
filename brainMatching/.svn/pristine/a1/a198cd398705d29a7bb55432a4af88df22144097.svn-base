/* 
 * File:   strCoupling_MI.h
 * Author: yusuf
 *
 * Generated on February 22, 2018, 11:29 AM
 */

#ifndef STRUCTURECOUPLING_MI_H
#define STRUCTURECOUPLING_MI_H

#include "matrixInitializer.h"

class StructureCoupling_MI : public MatrixInitializer
{
	public:
		StructureCoupling_MI(InitializerParameters &_initializerParameters):MatrixInitializer(_initializerParameters){}
		
		//initializers
		void initializeAssignmentCostMatrix(float **assignmentCostMatrix, Graph &objects, Graph &labels) override;
		void initializeObjectsWeightMatrix(float **objectsWeightMatrix, Graph &objects) override;
		void initializeLabelsDistanceMatrix(float **labelsDistanceMatrix, Graph &labels) override;
};

#endif /* STRUCTURECOUPLING_MI_H */

