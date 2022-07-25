/* 
 * File:   graph.h
 * Author: yusuf
 *
 * Generated on September 14, 2016, 4:15 PM
 */

#ifndef GRAPH_H
#define GRAPH_H

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <ctime>
#include <map>
#include <vector>
#include "edge.h"
#include "node.h"
#include "utility.h"

//#define NIL -1

class Graph
{
    public:
		Graph(bool isDirected,int _numFeatures, int _numExtraFeatures, int numOfNodes, float graphDensity);
        Graph(std::vector<Node> &_nodes, std::vector<Edge> &_edges);
        Graph(const Graph &graph);
        Graph(std::string fileName);
        Graph(std::string streamlineFile, std::string fmriNetworkFile,std::string bctFeaturesFile);
        Graph();
        ~Graph();
        Graph& operator=(const Graph& other);

        //functions for calculating connected components in the graph
        bool doesIncludeIsolatedNode(Node::Feature feature);
        void calculateNumberOfConnectedComponents(std::string connectivityType,std::vector<std::vector<int> > &connectedComponents);

        //shuffle graph
        void shuffleStructuralConnectivityPreservingStructuralNodeStrength(int seedSupplement);
        void shuffleStructuralConnectivityPreservingStructuralNodeDegree(int seedSupplement, int iteration);
        void shuffleFunctionalConnectivityPreservingSignedNodeDegree(int seedSupplement, int iteration);

        //normalize and logScale node/edge features (actions=logScaleNodes,logScaleEdges,normalizeNodes,normalizeEdges)
        void adjustGraph(std::string action);
        //initialize edge features from path calculations
        void calculatePath(std::string pathType,int maxPathLength=0);

        //distance related functions
        float calculateDiameter(Edge::Feature featureType, float **distanceMatrix);
        void  calculateAllPairsShortestPath(Edge::Feature featureType, float **distanceMatrix, int **predecessorMatrix,bool invert=false);
        void  calculateAllPairsShortestPath(Edge::Feature featureType, float **distanceMatrix,bool invert=false);
        static void  calculateAllPairsShortestPath(float **distanceMatrix, int size,bool invert=false);
        static void  calculateWeightedCommunicability(float **connectivityMatrix, std::vector<float> nodeStrengths, int size, int maxPathLength);
        static void  calculateUnweightedCommunicability(float **connectivityMatrix, int size, int maxPathLength);
        void  calculateSearchInformation(float **searchInformationMatrix,int size);
        void  calculatePathTransitivity(float **pathTransitivityMatrix,int size);

        //initialize connectivity matrix from edge features
        void  initializeWeightedConnectivityMatrixFromEdgeFeatures(Edge::Feature featureType, float **distanceMatrix);
        void  initializeBinaryConnectivityMatrixFromEdgeFeatures(Edge::Feature featureType, float **distanceMatrix);
        void  initializeNodeStrengthVectorFromEdgeFeatures(Edge::Feature featureType, std::vector<float> &nodeStrengths);
        void  initializeNodeDegreeVector(Edge::Feature featureType,std::vector<int> &nodeDegrees);
		
        //load/save functions
        void loadGraph(std::string fileName);
        void loadGraph(std::string streamlineFile, std::string fmriNetworkFile, std::string bctFeaturesFile="");
        void saveGraph(std::string fileName);
		
        //print functions
        void printSet(std::vector<Node> &set);
        void printSet(std::vector<Edge> &set);
        void print(std::string level);
        void printGraphAsAVector(std::ostream &out=std::cout);

        //accessors
        inline int getNumNodes() const{return nodes.size();}
        inline int getNumEdges() const{return edges.size();}
        inline std::string getSubjectId(){return subjectId;}
        inline std::vector<Node>::iterator getNodeIteratorBegin(){return nodes.begin();}
        inline std::vector<Node>::iterator getNodeIteratorEnd(){return nodes.end();}
        inline std::vector<Edge>::iterator getEdgeIteratorBegin(){return edges.begin();}
        inline std::vector<Edge>::iterator getEdgeIteratorEnd(){return edges.end();}

        //mutators
        inline void setNumFeatures( int dimFeat ){numFeatures = dimFeat;}
        inline void setSubjectId(std::string _subjectId){subjectId=_subjectId;}
      
    private:
        //constants
        const int BOUNDARY_FEATURE = 250;
        const float FUNC_CORR_RANGE = 0.999;
		
        //fields
        std::vector<Node> nodes;
        std::vector<Edge> edges;
		
        std::string subjectId;
        int   numFeatures; // number of features for each node
        int   numExtraFeatures;//number of extra features, that are calculated with another tool, such as BCT
        std::vector<std::string> extraFeatureNames;//names of the extra features, to keep record of (will use while saving the graph)
        bool  isDirected;       // is graph directed?    
        float diameter;         // diameter of the graph (i.e.,longest distance between any two pair of nodes)

        void initializeNodeFeatures(float **structuralConnectivity,float **functionalConnectivity, int numOfNodes);
        void initializeStructuralNodeFeatures(float **adjacencyMatrix, float threshold=Utility::EPSILON);
        void initializeFunctionalNodeFeatures(float **adjacencyMatrix, float threshold=Utility::EPSILON);
        void initializeExtraNodeFeatures(float **featuresMatrix);
        void initializeEdgeFeatures(float **structuralConnectivity,float **functionalConnectivity);

        //initialize edge features from path calculations
        void  initializeEdgeFeaturesWithDirectEdges(Edge::Feature featureType);
        void  initializeEdgeFeaturesWithStrongestSumPath(Edge::Feature featureType);
        void  initializeEdgeFeaturesWithWeightedShortestPath(Edge::Feature featureType);
        void  initializeEdgeFeaturesWithUnweightedShortestPath(Edge::Feature featureType);
        void  initializeEdgeFeaturesWithWeightedCommunicabilityPath(Edge::Feature featureType, int maxPathLength);
        void  initializeEdgeFeaturesWithUnweightedCommunicabilityPath(Edge::Feature featureType, int maxPathLength);
        void  initializeEdgeFeaturesWithSearchInformation(Edge::Feature featureType);
        void  initializeEdgeFeaturesWithPathTransitivity(Edge::Feature featureType);

        void normalizeEdgeFeatures();
        void normalizeNodeFeatures();
        void logScaleNodeFeatures(std::string action);
        void logScaleEdgeFeatures(std::string action);
        void logScaleStructuralGraph();
        void zScoreStructuralEdgeFeatures();
        void scaleGraph(float ratio); //scale the nodes and edges of the graph wrt the given ratio
        void binarizeEdgesByValue(std::string modality,float threshold); // set edges less than a threshold to zero, and the rest to 1
        void binarizeEdgesByDensity(std::string modality,float density); // set edges less than a threshold to zero, and the rest to 1
        void thresholdEdgesByValue(std::string modality,float threshold); // set edges less than a threshold to zero, and leave the rest as they are
        void thresholdEdgesByDensity(std::string modality,float density); // keep a certain density of nonzero edges that are larger than the rest intact and set the rest to zero
		
        //base code related functions
        bool checkForZeroWeightEdge();
        bool isDuplicateNode(Node &node);//returns true if "Node node" is a duplicate node
        bool checkTriangleInequality();
		
        float getValueOfMinEdge(Edge::Feature featureType); //returns the distance of the shortest edge
		
        void makeCompleteGraphFromExistingNodes();//makes a complete graph by inserting O(n^2) edges between the already existing nodes
        void setRandomGraph(bool _isDirected,int _numFeatures, int _numExtraFeatures, int numOfNodes, float graphDensity);
        void generateRandomNodes(int numOfNodes); //generates random nodes of given size
        void generateRandomEdges(int numOfEdges); //given a graph with existing nodes, generates numOfEdges random edges
        void getSubgraph(int firstNNodes, Graph &targetGraph);//copies the subgraph consisting of first N nodes of this graph to the provided targetGraph
        void getSubgraph(std::vector<int> &nodeList, Graph &grp);//returns a subgraph of input graph grp consisting of nodes whose IDs are listed in nodesList vector	
};

#endif /* GRAPH_H */

