/* 
 * File:   edge.h
 * Author: yusuf
 *
 * Generated on September 14, 2016, 2:51 PM
 */

#ifndef EDGE_H
#define EDGE_H

#include <vector>
#include <cmath>    //fabs()
#include <cstdlib>  //rand(), exit()
#include <iostream> //std::cout
#include <fstream>  //std::ofstream, std::ifstream
#include <sstream>  //std::istringstream
#include "node.h"

class Edge
{
    public:
        Edge();
        Edge(const Edge &other);
        Edge(int _edgeId, Node &p1, Node &p2);
        Edge(int _edgeId, int _node1Id, int _node2Id, float spatialDistance, float structuralConnectedness, float functionalCorrelation);
        Edge& operator=(const Edge& other);
		
        enum Feature{SPATIAL_DISTANCE=0,STRUCTURAL_CONNECTIVITY,FUNCTIONAL_CONNECTIVITY};

        void print(std::ostream &out=std::cout);
        void saveEdge(std::ofstream &file);
        void loadEdge(std::ifstream &file);

        inline bool isEqual(const Edge & e){return (node1Id == e.node1Id && node2Id == e.node2Id)?true:false;}
        inline bool isEquivalent(const Edge & e){return ((node1Id == e.node1Id && node2Id == e.node2Id)||(node1Id == e.node2Id && node2Id == e.node1Id))?true:false;}
		
        inline int   getEdgeId(){return edgeId;}
        inline int   getNode1Id(){return node1Id;}
        inline int   getNode2Id(){return node2Id;}
        inline float getFeature(int featureId){return features[featureId];}

        inline void setFeature(int featureId, float value){features[featureId]=value;}
        inline void scaleFeature(int featureId, float ratio){features[featureId] *= ratio;}
        inline void logScaleFeature(int featureId){features[featureId] = (features[featureId]<1 ? 0 : (float)std::log(features[featureId])/(float)std::log(2));}

        inline void setNodeIds(int n1, int n2){node1Id=n1;node2Id=n2;}

        static Feature getEdgeType(std::string edgeTypeString);

    private:
        int edgeId;
        int node1Id, node2Id;
        const int numFeatures = 3;
        //FEATURES
        //spatial_distance = Euclidean distance between the locations of nodes
        //structural_connectedness = number of streamlines between two nodes
        //functional_correlation = correlation between two nodes according to the fMRI data
        std::vector<float> features;
};

#endif /* EDGE_H */

