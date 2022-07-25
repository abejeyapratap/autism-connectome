/* 
 * File:   edge.cpp
 * Author: yusuf
 *
 * Generated on September 14, 2016, 3:43 PM
 */

#include "edge.h"
#include "utility.h"

using namespace std;

Edge::Edge()
{
    for(int i=0;i<numFeatures;i++)
        features.push_back(-1);
}

Edge::Edge(const Edge &other)
{
    edgeId = other.edgeId;
    node1Id = other.node1Id;
    node2Id = other.node2Id;
    features = other.features;
}

Edge::Edge(int _edgeId, Node &p1, Node &p2) : Edge()
{
    edgeId = _edgeId;
    node1Id = p1.getNodeId();
    node2Id = p2.getNodeId();

    features[SPATIAL_DISTANCE] = p1.calculateDistance(p2);
    features[STRUCTURAL_CONNECTIVITY] = 1;
    features[FUNCTIONAL_CONNECTIVITY] = 0;
}

Edge::Edge(int _edgeId, int _node1Id, int _node2Id, float spatialDistance, float structuralConnectedness, float functionalCorrelation) : Edge()
{
    edgeId = _edgeId;
    node1Id = _node1Id;
    node2Id = _node2Id;
   
    features[SPATIAL_DISTANCE] = spatialDistance;
    features[STRUCTURAL_CONNECTIVITY] = structuralConnectedness;
    features[FUNCTIONAL_CONNECTIVITY] = functionalCorrelation;
}

Edge& Edge::operator=(const Edge &other)
{
    edgeId = other.edgeId;
    node1Id = other.node1Id;
    node2Id = other.node2Id;
    features = other.features;
    return *this;
}

Edge::Feature Edge::getEdgeType(std::string edgeTypeString)
{
   Feature edgeType;
   if(edgeTypeString.compare("structure")==0 || edgeTypeString.compare("str")==0)
      edgeType=Edge::STRUCTURAL_CONNECTIVITY;
   else if(edgeTypeString.compare("function")==0 || edgeTypeString.compare("func")==0)
      edgeType=Edge::FUNCTIONAL_CONNECTIVITY;
   else
   {
      cerr<<"incorrect edge type in Edge::getEdgeType:"<<edgeTypeString<<"... Exiting..\n";
      exit(1);
   }
   return edgeType;
}

//print the contents of an edge to the screen
void Edge::print(std::ostream &out)
{
    //format: edgeId <tab> node1Id <tab> node2Id <tab> distanceOfEdge <tab> weightOfEdge <tab> functionalCorrelation
    out<< edgeId<< ":\t(" << node1Id << "," << node2Id <<")\t dist:"<<features[SPATIAL_DISTANCE]
        <<"\t w:"<<features[STRUCTURAL_CONNECTIVITY]<<"\t func:"<<features[FUNCTIONAL_CONNECTIVITY];
}

//write the contents of a node into a output stream (i.e., file)
void Edge::saveEdge(ofstream &file)
{
   //format is as follows
   //<edge id> <node1> <node2> <distance> <weight> <functionalCorrelation>
   //<to,from,pos>#<to,from,pos>#...
   
   file<<edgeId<<"\t"<<node1Id<<"\t"<<node2Id<<"\t";
   
   file.precision(4);
   for(int i=0;i<numFeatures;i++)
      file<<features[i]<<"\t";
   file<<endl;
}

//load a node from input stream (i.e., file)
void Edge::loadEdge(ifstream &file)
{
   istringstream is; 
   string line;
  
   line.clear();
   is.clear();
   getline(file, line);
   is.str(line);
   
   is >> edgeId >> node1Id >> node2Id;
   
   float tempVal;
   for(int i=0;i<numFeatures;i++)
   {
      is >> tempVal;
      tempVal = Utility::absoluteValue(tempVal)<Utility::EPSILON?Utility::EPSILON:tempVal;//to avoid divide by zero, we avoid features having value zero
      features[i] = tempVal;
   }
}
