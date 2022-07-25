/* 
 * File:   node.cpp
 * Author: yusuf
 *
 * Generated on September 14, 2016, 2:03 PM
 */

#include "node.h"
#include "utility.h"

using namespace std;

Node::Node(const Node& other)
{
    nodeId = other.nodeId;
    features = other.features;
   
    numFeatures=other.numFeatures;
    numExtraFeatures=other.numExtraFeatures;
}

Node::Node(int _numFeatures, int _numExtraFeatures)
{
    numFeatures = _numFeatures;
    numExtraFeatures = _numExtraFeatures;
}
 
Node::Node(int _nodeId, int _numFeatures, int _numExtraFeatures)
{
    nodeId = _nodeId;
    numFeatures = _numFeatures;
    numExtraFeatures = _numExtraFeatures;
    features.reserve(numFeatures);
    for(int i=0;i<numFeatures;i++)
        features.push_back(-1);
}

Node& Node::operator=(const Node& other)
{
    nodeId = other.nodeId;
    features = other.features;
   
    numFeatures=other.numFeatures;
    numExtraFeatures=other.numExtraFeatures;

    return *this;
}

//keep this function as a place holder
float Node::calculateDistance(Node &other)
{
   return 0;
}

//randomly generates the location/features of a node
//takes the dimension of the node (i.e. 2d, 3d etc) and an upper limit for each dimension to be set to.
void Node::randomlySetNode(int upperLimitFeature)
{
   //randomly set features within the range of [0,UpperLimitFeature]
   Utility::randomlySetVector(features,numFeatures,upperLimitFeature);
    
   //randomly set features among a list of available values
   //vector<float> values{50.0, 125.0, 250.0}; //list of values to be randomly assigned to the features
   //Utility::randomlySetVector(features,values,dimensionFeature);
}

void Node::getExtraFeatures(std::vector<float>& extraFeaturesCopy)
{
    extraFeaturesCopy.clear();
    for(int i=NUM_STR_FEATURES+NUM_FUNC_FEATURES;i<numFeatures;i++)
        extraFeaturesCopy.push_back(features[i]);
}

//print the contents (location and features) of a node to the screen
void Node::print(string str, std::ostream &out)
{
    if(str.compare("id")==0 || str.compare("all")==0)
        out << nodeId << "\t";
   
    int precision = 6;

    if(str.compare("features")==0 || str.compare("all")==0)
        Utility::printVector(features,out,'\t','\n',precision);
}

//write the contents of a node into a output stream (i.e., file)
void Node::saveNode(ofstream &file)
{
    //format is as follows
    //<node id> <location> <features>
    file<<nodeId<<"\t";

    int precision = 4;

    Utility::printVector(features,file,'\t','\n',precision);
}

//load a node from input stream (i.e., file)
void Node::loadNode(ifstream& file)
{
    istringstream is; 
    string line;

    /////read node id, location, and features
    line.clear();
    is.clear();
    getline(file, line);
    is.str(line);

    is >> nodeId;
    float tempVal;

    features.reserve(numFeatures);//preallocate space for vectors so that we do not copy the stuff  again and again as we keep adding new features
    for(int j=0;j<numFeatures;j++)
    {
        is >> tempVal;
        tempVal = Utility::absoluteValue(tempVal)<Utility::EPSILON?Utility::EPSILON:tempVal;//to avoid divide by zero, we avoid features having value zero
        features.push_back(tempVal);
    }
}

