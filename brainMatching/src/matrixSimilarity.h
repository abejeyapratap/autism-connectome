/* 
 * File:   matrixSimilarity.h
 * Author: yusuf
 *
 * Created on November 6, 2019, 2:23 PM
 */

#ifndef MATRIXSIMILARITY_H
#define MATRIXSIMILARITY_H

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <algorithm> // std::find
#include <dirent.h>  //for readdir, opendir, closedir
#include "dataset.h"
//#include "graph.h"
#include "edge.h"
//#include "normalizer.h"
#include "parameterSet.h"
#include "initializerParameters.h"
//#include "matrixInitializer.h"


class MatrixSimilarity
{
	public:
		MatrixSimilarity(const MatrixSimilarity &other){dataset= other.dataset;dataset2=other.dataset2;};
		MatrixSimilarity(Dataset *_dataset,Dataset *_dataset2=NULL){dataset= _dataset;dataset2=_dataset2;};
		~MatrixSimilarity(){}

		void correlateGroups(float** rValues, std::vector<int> &group1Orders, std::vector<int> &group2Orders, Edge::Feature edgeType1, Edge::Feature edgeType2);
		void correlateEveryoneWithItself(std::vector<float> &rValues,Edge::Feature edgeType1, Edge::Feature edgeType2);

		void calculateMatrixDistanceOfEveryoneToItself(Edge::Feature edgeType1, Edge::Feature edgeType2,std::vector<float> &similarityScores, std::string distanceMeasure);
		void calculateMatrixDistanceBetweenGroups(float** similarityScores, std::vector<int> &group1IDs, std::vector<int> &group2IDs, std::string distanceMeasure, Edge::Feature edgeType1, Edge::Feature edgeType2);
		void subtractEveryoneFromMembersOfGroup(float*** dissimilarityMatrices, std::vector<int> &groupRowOrders, std::vector<int> &groupColumnOrders,std::vector<std::string> &groupRowIDs, std::vector<std::string> &groupColumnIDs, Edge::Feature edgeType1, Edge::Feature edgeType2);
		void calculateZscoreOfEdgesForEveryoneWRTGroup(float*** dissimilarityMatrices, std::vector<int> &groupRowOrders, std::vector<int> &groupColumnOrders,std::vector<std::string> &groupRowIDs, std::vector<std::string> &groupColumnIDs, Edge::Feature edgeType1, Edge::Feature edgeType2);
	
	private:     
		Dataset *dataset,*dataset2;
};


#endif /* MATRIXSIMILARITY_H */

