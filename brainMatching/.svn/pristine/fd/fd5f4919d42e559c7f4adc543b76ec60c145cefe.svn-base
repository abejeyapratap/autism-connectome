/* 
 * File:   node.h
 * Author: yusuf
 *
 * Generated on September 14, 2016, 1:52 PM
 */

#ifndef NODE_H
#define NODE_H

#include <vector>
#include <cmath>    //fabs()
#include <cstdlib>  //rand(), exit()
#include <iostream> //std::cout
#include <fstream>  //std::ofstream, std::ifstream
#include <sstream>  //std::istringstream
#include "time.h"
#include "geometry.h" 

class Node
{
    public:
        Node(const Node &other);
        Node(){}
        Node(int _numFeatures, int _numExtraFeatures);
        Node(int _nodeId, int _numFeatures, int _numExtraFeatures);
        Node& operator=(const Node& other);
		
        enum Feature{STR_NODE_DEGREE=0,STR_NODE_STRENGTH,
                    FUNC_NODE_DEGREE_NEG,FUNC_NODE_DEGREE_POS,
                    FUNC_NODE_STRENGTH_NEG,FUNC_NODE_STRENGTH_POS};
        enum {NUM_STR_FEATURES=2, NUM_FUNC_FEATURES=4};
		
        float calculateDistance(Node &other);
        inline void  scaleFeature(int featureId,float ratio){features[featureId] *= ratio;};
        inline void  logScaleFeature(int featureId){features[featureId] = (features[featureId]<1 ? 0 : (float)std::log(features[featureId])/(float)std::log(2));};

        void print(std::string str, std::ostream &out=std::cout);
        void saveNode(std::ofstream &file);
        void loadNode(std::ifstream &file);

        void getExtraFeatures(std::vector<float> &extraFeaturesCopy);

        inline int   getNodeId(){return nodeId;}
        inline float getFeature(int featureId){return features[featureId];}
        inline int   getNumberOfFeatures(){return numFeatures;}
        inline int   getNumberOfExtraFeatures(){return numExtraFeatures;}
        inline int   getIndexOfStructuralFeatures(){return 0;}
        inline int   getIndexOfFunctionalFeatures(){return NUM_STR_FEATURES;}
        inline int   getIndexOfExtraFeatures(){return NUM_STR_FEATURES+NUM_FUNC_FEATURES;}

        inline void setFeature(int featureId, float value){features[featureId]=value;}

        void randomlySetNode(int upperLimitFeature);

    private:
        int nodeId;
		
        int numFeatures, numExtraFeatures;
        //FEATURES
        //structural node strength/degree (0,1) = number of fibers connected to this node/number of regions connected to this node
        //functional node strength/degree (2,3,4,5) = strength/degree of negative/positive and functional correlation between this node and other nodes
        //BCTfeatures (6,7,8,9,10,11) = <degree> <strength> <betweenness centrality> <local efficiency> <participation coefficient> <local assortativity> 
        std::vector<float> features;
};

#endif /* NODE_H */

