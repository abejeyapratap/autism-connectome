/* 
 * File:   normalizer.h
 * Author: yusuf
 *
 * Generated on October 26, 2016, 2:45 PM
 */

#ifndef NORMALIZER_H
#define NORMALIZER_H

#include <map>
#include "graph.h"

//This class keeps record of the data range for features of a set of graphs
//Main purpose of the calculations in here is to normalize the features into a common
//interval, such as [0,1], for all feature types.
//Note that, this class does not modify the edges or he nodes of the graphs.
//It is programmers responsibility to make such adjustments later on based on the 
//data ranges that are calculated by this class.
//In the implementation of mutimodalBrainMatch project, the numbers calculated by this
//class are being utilized in the matrixInitializer class while calculating the 
//terms present in the objective function of the optimization.
class Normalizer
{
	public:
		Normalizer();
		Normalizer(const Normalizer &normalizer);
		
		void normalizeData(std::map<int,Graph> &graphs);
		void print();
		
		void save(std::string filename);
		void load(std::string filename);
		
		inline float getSpatialRange(){return spatialRange;}
		inline float getMinSpatial(){return minSpatial;}
		
		inline float getStrNodeStrengthRange(){return strNodeStrengthRange;}
		inline float getMinStrNodeStrength(){return minStrNodeStrength;}
		
		inline float getFuncNodeStrengthRange(int featureId){return funcNodeStrengthRange[featureId];}
		inline float getMinFuncNodeStrength(int featureId){return minFuncNodeStrength[featureId];}
		inline float getMaxFuncNodeStrength(int featureId){return maxFuncNodeStrength[featureId];}
		
		inline float getSpatialDistanceRange(){return spatialDistanceRange;}
		inline float getMinSpatialDistance(){return minSpatialDistance;}
		inline float getMaxSpatialDistance(){return maxSpatialDistance;}
		
		inline float getStructuralConnectednessRange(){return structuralConnectednessRange;}
		inline float getMinStructuralConnectedness(){return minStructuralConnectedness;}
		inline float getMaxStructuralConnectedness(){return maxStructuralConnectedness;}
		
		inline float getFunctionalCorrelationRange(){return functionalCorrelationRange;}
		inline float getMinFunctionalCorrelation(){return minFunctionalCorrelation;}
		inline float getMaxFunctionalCorrelation(){return maxFunctionalCorrelation;}
	
	private:
		//ranges for node features
		float minSpatial, maxSpatial, spatialRange;
		float minStrNodeStrength, maxStrNodeStrength, strNodeStrengthRange;
		const static int NUM_FUNC_NODE_FEATURES=4;
		float minFuncNodeStrength[NUM_FUNC_NODE_FEATURES], maxFuncNodeStrength[NUM_FUNC_NODE_FEATURES], funcNodeStrengthRange[NUM_FUNC_NODE_FEATURES];
		
		//ranges for edge features
		float minSpatialDistance, maxSpatialDistance, spatialDistanceRange;
		float minStructuralConnectedness, maxStructuralConnectedness, structuralConnectednessRange;
		float minFunctionalCorrelation, maxFunctionalCorrelation, functionalCorrelationRange;
		
		//ranges for brainGraphFeatures
		float minBrainVolume, maxBrainVolume, brainVolumeRange;
		
	//////functions for a set off graphs
		//normalize node features across a set of graphs
		void normalizeNodeFeatures(std::map<int,Graph>& graphs);
		void normalizeSpatial(std::map<int,Graph>& graphs);
		void normalizeStructuralNodeStrength(std::map<int,Graph>& graphs);
		void normalizeFunctionalNodeStrength(std::map<int,Graph>& graphs);
		
		//normalize edge features across a set of graphs
		void normalizeSpatialDistance(std::map<int,Graph>& graphs);
		void normalizeStructuralConnectedness(std::map<int,Graph>& graphs);
		void normalizeFunctionalCorrelation(std::map<int,Graph>& graphs);
		
		//normalize edges with respect to pairwise relations across a set of graphs
		void normalizeStructuralConnectednessPairwise(std::map<int,Graph>& graphs);
		void normalizeSpatialDistancePairwise(std::map<int,Graph>& graphs);
		
	//////functions for a single graph
		//normalize edge features for the edges of a single graph
		void normalizeSpatialDistance(Graph &graph);
		void normalizeStructuralConnectedness(Graph &graph);
		void normalizeFunctionalCorrelation(Graph &graph);
};

#endif /* NORMALIZER_H */

