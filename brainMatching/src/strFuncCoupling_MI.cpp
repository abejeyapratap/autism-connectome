/* 
 * File:   strFuncCoupling_MI.cpp
 * Author: yusuf
 *
 * Generated on February 22, 2018, 11:23 AM
 */

#include "strFuncCoupling_MI.h"
#include "node.h"
#include "edge.h"
#include "utility.h"

using namespace std;

//<editor-fold defaultstate="collapsed" desc=" Assignment Cost: initializeAssignmentCostMatrix()">

//calculates "c(p,a)" assignment cost of an object p to a label a for all labels and objects 
//as the Euclidean distance between description of label features and description of object features
//Note: structural and functional connectomes MUST have been edgeNormalized before calling this function
void StrFuncCoupling_MI::initializeAssignmentCostMatrix(float **assignmentCostMatrix, Graph &objects, Graph &labels)
{
   //define a reasonable assignment cost function in here
   float assignmentCost;
   int numNodes = objects.getNumNodes();
   
   Utility::fillMatrix<float>(assignmentCostMatrix,0,numNodes,numNodes);
   
   if(initializerParameters.assignmentCostMode.find("edges")!=std::string::npos)
   {
      float **dist1=Utility::allocate2Dmemory<float>(numNodes,numNodes);
      float **dist2=Utility::allocate2Dmemory<float>(numNodes,numNodes);

      srand((int)clock());
      //STEP 1:adjust structural and functional connectivity matrices
      ////load structural
      objects.initializeWeightedConnectivityMatrixFromEdgeFeatures(Edge::STRUCTURAL_CONNECTIVITY,dist1);
      ////load functional
      labels.initializeWeightedConnectivityMatrixFromEdgeFeatures(Edge::FUNCTIONAL_CONNECTIVITY,dist2);
      //However, we will use either the positive or the negative correlations. Thus, we expect to get a half populated matrix
      //get rid of the negative (or positive) correlations
      if(initializerParameters.functionalConnectivity.compare("negative")==0)
      {//take the absolute value of the negative values only, discard the rest
         for(int i=0;i<numNodes;i++)
            for(int j=0;j<numNodes;j++)
               dist2[i][j] = dist2[i][j]>0 ? 0 : -dist2[i][j];
      }
      else if(initializerParameters.functionalConnectivity.compare("positive")==0)
         Utility::filterOutElementsOfMatrixLessThanThreshold<float>(dist2,0.0,0.0,numNodes,numNodes);
      else
      {
         cerr<<"negative/positive is not defined for functional connectivity. Check the -funcConn flag. Exiting...\n";
         exit(1);
      }

      if(initializerParameters.assignmentCostMode.find("IncludeDiag")!=std::string::npos)
      {
         if(initializerParameters.assignmentCostMode.find("ZeroDiag")!=std::string::npos)
         {
            Utility::fillDiagonalOfTheMatrix<float>(dist1,0.0,numNodes);//fill diagonals with zero
            Utility::fillDiagonalOfTheMatrix<float>(dist2,0.0,numNodes);//fill diagonals with zero
         }
         else if(initializerParameters.assignmentCostMode.find("RandDiag")!=std::string::npos)
         {
            Utility::fillDiagonalOfTheMatrixRandomlyColumnwise<float>(dist1,numNodes,rand());//fill diagonals with random values. Third parameters is a seed value
            Utility::fillDiagonalOfTheMatrixRandomlyColumnwise<float>(dist2,numNodes,rand());//fill diagonals with random values. Third parameters is a seed value
         }
      }

      //Utility::saveMatrixToFile<float>("structuralConnectivity.txt",dist1,numNodes,numNodes,'\t','\n',8);//DEBUG
      //Utility::saveMatrixToFile<float>("functionalConnectivity.txt",dist2,numNodes,numNodes,'\t','\n',8);//DEBUG

      //STEP 2: once we generate the connectivity matrices for each node, we will use their rows as the connectivity feature vector
      //of nodes. Below, we calculate the Euclidean distance between these feature vectors as the similarity (assignment cost)
      //of a structural node to a functional node
      if(initializerParameters.assignmentCostMode.find("IncludeDiag")!=std::string::npos)
      {
         for(int i=0;i<numNodes;i++)
            for(int j=0;j<numNodes;j++)
               assignmentCostMatrix[i][j]+=Utility::calculateL2Distance(dist1[i],dist2[j],numNodes);
      }
      else if(initializerParameters.assignmentCostMode.find("IgnoreDiag")!=std::string::npos)
      {
         for(int i=0;i<numNodes;i++)
            for(int j=0;j<numNodes;j++)
               assignmentCostMatrix[i][j]+=Utility::calculateL2DistanceByIgnoringPairwiseRelations<float>(dist1[i],dist2[j],numNodes,i,j); 
      }

      //We don't the assignment cost matrix This is being handled in doLabeling function 
      
      Utility::free2Dmemory<float>(dist1,numNodes);
      Utility::free2Dmemory<float>(dist2,numNodes);
   }
   
   if(initializerParameters.assignmentCostMode.find("features")!=std::string::npos)
   {
      //since we already guaranteed that we have three nonnegative coefficient for each assignmentCost feature,
      //we can now calculate the assignment cost as a mixture of the three features (if all coefficients are positive)
      for ( vector<Node>::iterator objectsIter = objects.getNodeIteratorBegin(); objectsIter != objects.getNodeIteratorEnd(); ++objectsIter )
      {
         int objId = objectsIter->getNodeId();
         vector<float> features1;
         objectsIter->getExtraFeatures(features1);
         for ( vector<Node>::iterator labelsIter = labels.getNodeIteratorBegin(); labelsIter != labels.getNodeIteratorEnd(); ++labelsIter )
         {
            vector<float> features2;
            labelsIter->getExtraFeatures(features2);

            assignmentCost = Utility::calculateL1Distance<float>(features1,features2);
            assignmentCostMatrix[objId][labelsIter->getNodeId()] += assignmentCost;
         }
      }
   }
   
//   Utility::saveMatrixToFile<float>("assignmentCostMatrix.txt",assignmentCostMatrix,numNodes,numNodes,'\t','\n',8);//DEBUG
//   exit(1);//DEBUG
}
//</editor-fold>

//<editor-fold defaultstate="collapsed" desc=" Separation Cost: initializeObjectsWeightMatrix(), initializeLabelsDistanceMatrix()">

//Separation cost consists of two parts: w_{p,q} * d_{a,b} where (p,q) is an edge between object nodes p,q and 
//and (a,b) is an edge between label nodes a and b, where p,q are mapped to labels a,b respectively via a mapping function f.
// - Value w_{a,b} represents the strength of relation between the two endpoints of a and b of an edge e=(a,b). 
//   w_{a,b} will get a large value if nodes a and b are similar to each other
// - Value d_{a,b} is the distance between the points that p and q get mapped to in the label graph.
//   d_{a,b} will get a large value if nodes x and y are highly dissimilar.
//In a brain graph, we have three measures that we can use for defining pairwise relations between nodes:
//1) spatialDistance between two regions: which is large if the two nodes are spatially separated far away 
//   (in brain graphs, this is Euclidean distance between the locations of the nodes)
//2) structuralConnectedness of two regions, that is the number of streamlines between two regions. 
//   This value is large if the two nodes are similar/strongly correlated
//3) functionalCorrelation between two regions: which is large if the two regions are strongly correlated.
//   (In brain graphs, this is obtained by the fmri data)
//w_{p,q} will be  1/spatialDistance(p,q), structuralConnectedness(p,q), and/or functionalCorrelation(p,q), respectively.
//d_{a,b} will be spatialDistance(p,q), 1/structuralConnectedness(p,q), and/or 1/functionalCorrelation(p,q), respectively. 
//
//Normalization:
//initializeObjectsWeightMatrix, 
//  1) we normalize each spatialDist_{p,q} by making the following update: w_{p,q} = (1/spatialDist_{p,q} - 1/spatialDist_{max}) / (1/spatialDist_{min} - 1/spatialDist_{max})
//  2) we normalize each strConn_{p,q} by making the following update: w_{p,q} = (strConn_{p,q} - strConn_{min}) / (strConn_{max}-strConn_{min})
//  3) we normalize each funcCorr_{p,q} by making the following update: w_{p,q} = (funcCorr_{p,q} - funcCorr_{min}) / (funcCorr_{max}-funcCorr_{min})
void StrFuncCoupling_MI::initializeObjectsWeightMatrix(float **objectsWeightMatrix,Graph &objects)
{
//   float minSpatialDistance = 1.0/(1.0+normalizer.getMaxSpatialDistance());
//   float maxSpatialDistance = 1.0/(1.0+normalizer.getMinSpatialDistance());
//   float spatialDistanceRange = maxSpatialDistance-minSpatialDistance;
//   
//   float minStructuralConnectedness = normalizer.getMinStructuralConnectedness();
//   float structuralConnectednessRange = normalizer.getStructuralConnectednessRange();
//   
//   float minFunctionalCorrelation = normalizer.getMinFunctionalCorrelation();
//   float functionalCorrelationRange = normalizer.getFunctionalCorrelationRange();
   
   //DEBUG: print the ranges of the edge weights
//   cerr<<"----------------initializeObjectsWeightMatrix ranges------------\n";
//   cerr<<"minSpatialDistance:"<<minSpatialDistance<<"\tmaxSpatialDistance:"<<maxSpatialDistance<<"\tspatialDistanceRange:"<<spatialDistanceRange<<endl;
//   cerr<<"minStructuralConnectedness:"<<minStructuralConnectedness<<"\tstructuralConnectednessRange:"<<structuralConnectednessRange<<endl;
//   cerr<<"minFunctionalCorrelation:"<<minFunctionalCorrelation<<"\tfunctionalCorrelationRange:"<<functionalCorrelationRange<<endl;

//   for(std::vector<float>::iterator iter=initializerParameters.parameters.begin();iter!=initializerParameters.parameters.end();iter++)
//   {
//      if(*iter<0)
//      {
//         cerr<<"separationParameters are not set while the separationCostMode is set to be mixed!! Exiting...\n";
//         exit(1);
//      }
//   }

//   if(spatialDistanceRange<Utility::EPSILON || structuralConnectednessRange<Utility::EPSILON || functionalCorrelationRange<Utility::EPSILON)
//   {
//      cerr<<"one of the ranges for the edge features are almost equal to zero. This will cause divide by zero error. Exiting!!!"<<endl;
//      cerr<<"spatialDistanceRange:"<<spatialDistanceRange<<"\tstructuralConnectednessRange:"<<structuralConnectednessRange<<"\tfunctionalCorrelationRange:"<<functionalCorrelationRange<<endl;
//      exit(1);
//   }
//      
   float w_pq;
   float ratio = objects.getNumNodes()*objects.getNumNodes();
   for ( vector<Edge>::iterator edgeIter = objects.getEdgeIteratorBegin(); edgeIter != objects.getEdgeIteratorEnd(); ++edgeIter )
   {
      //to avoid divide by zero in the calculation of spatialDistance, we take this precaution.
//      float spatial = edgeIter->getFeature(Edge::SPATIAL_DISTANCE);
      float structural = edgeIter->getFeature(Edge::STRUCTURAL_CONNECTIVITY);
//      float functional = edgeIter->getFeature(Edge::FUNCTIONAL_CONNECTIVITY);
      
      //use following line if you would like to take weight = 1/distance
      //spatial = spatial  < Utility::EPSILON ? 1 : (1.0/(1.0+spatial)-minSpatialDistance)/spatialDistanceRange;
      
      //use following line if you would like to take weight = 1-distance
//      spatial = 1- spatial;

      //w_pq = initializerParameters.getParameter(ParameterSet::SPATIAL_DISTANCE) * spatial + 
        w_pq = structural/ratio; //initializerParameters.getParameter(ParameterSet::STRUCTURAL_CONNECTIVITY) * structural +
      //       initializerParameters.getParameter(ParameterSet::FUNCTIONAL_CONNECTIVITY) * functional;

      objectsWeightMatrix[edgeIter->getNode1Id()][edgeIter->getNode2Id()] = w_pq;
      objectsWeightMatrix[edgeIter->getNode2Id()][edgeIter->getNode1Id()] = w_pq;
   }
}

//Normalization:
//initializeLabelsDistanceMatrix, 
//  1) we normalize each spatialDist_{a,b} by making the following update: d_{a,b} = (spatialDist_{a,b} - spatialDist_{min}) ) / (spatialDist_{max}-spatialDist_{min})
//  2) we normalize each strConn_{a,b} by making the following update: d_{a,b} = (1/strConn_{a,b} - 1/strConn_{max}) / (1/strConn_{min} - 1/strConn_{max})
//  3) we normalize each funcCorr_{a,b} by making the following update: d_{a,b} = (1/funcCorr_{a,b} - 1/funcCorr_{max}) / (1/funcCorr_{min} - 1/funcCorr_{max})
void StrFuncCoupling_MI::initializeLabelsDistanceMatrix(float** labelsDistanceMatrix, Graph& labels)
{
//   float minSpatialDistance = normalizer.getMinSpatialDistance();
//   float spatialDistanceRange = normalizer.getSpatialDistanceRange();
//         
//   float minStructuralConnectedness = 1.0/(1.0+normalizer.getMaxStructuralConnectedness());
//   float maxStructuralConnectedness = 1.0/(1.0+normalizer.getMinStructuralConnectedness());
//   float structuralConnectednessRange = maxStructuralConnectedness - minStructuralConnectedness;
//   
//   float minFunctionalCorrelation = 1.0/(1.0+normalizer.getMaxFunctionalCorrelation());
//   float maxFunctionalCorrelation = 1.0/(1.0+normalizer.getMinFunctionalCorrelation());
//   float functionalCorrelationRange = maxFunctionalCorrelation - minFunctionalCorrelation;
   
   //DEBUG: print the ranges of the edge weights
//   cerr<<"----------------initializeLabelsDistanceMatrix ranges------------\n";
//   cerr<<"minSpatialDistance:"<<minSpatialDistance<<"\tspatialDistanceRange:"<<spatialDistanceRange<<endl;
//   cerr<<"minStructuralConnectedness:"<<minStructuralConnectedness<<"\tmaxStructuralConnectedness:"<<maxStructuralConnectedness<<"\tstructuralConnectednessRange:"<<structuralConnectednessRange<<endl;
//   cerr<<"minFunctionalCorrelation:"<<minFunctionalCorrelation<<"\tmaxFunctionalCorrelation:"<<maxFunctionalCorrelation<<"\tfunctionalCorrelationRange:"<<functionalCorrelationRange<<endl;

//   for(std::vector<float>::iterator iter=initializerParameters.separationParameters.begin();iter!=initializerParameters.separationParameters.end();iter++)
//   {
//      if(*iter<0)
//      {
//         cerr<<"separationParameters are not set while the separationCostMode is set to be mixed!! Exiting...\n";
//         exit(1);
//      }
//   }
   
//   if(spatialDistanceRange<Utility::EPSILON || structuralConnectednessRange<Utility::EPSILON || functionalCorrelationRange<Utility::EPSILON)
//   {
//      cerr<<"one of the ranges for the edge features are almost equal to zero. This will cause divide by zero error. Exiting!!!"<<endl;
//      cerr<<"spatialDistanceRange:"<<spatialDistanceRange<<"\tstructuralConnectednessRange:"<<structuralConnectednessRange<<"\tfunctionalCorrelationRange:"<<functionalCorrelationRange<<endl;
//      exit(1);
//   }
   
   float d_ab;
   float ratio = labels.getNumNodes()*labels.getNumNodes();
   for ( vector<Edge>::iterator edgeIter = labels.getEdgeIteratorBegin(); edgeIter != labels.getEdgeIteratorEnd(); ++edgeIter )
   {
//      float spatial = edgeIter->getFeature(Edge::SPATIAL_DISTANCE);
//      float structural = edgeIter->getFeature(Edge::STRUCTURAL_CONNECTIVITY);
      float functional = edgeIter->getFeature(Edge::FUNCTIONAL_CONNECTIVITY);

//      spatial = (spatial-minSpatialDistance)/spatialDistanceRange;
      
      //use following two lines if you would like to take weight = 1/distance
      //structural = structural < Utility::EPSILON ? 1 : (1.0/(1.0+structural)-minStructuralConnectedness)/structuralConnectednessRange;
      //functional = functional < Utility::EPSILON ? 1 : (1.0/(1.0+functional)-minFunctionalCorrelation)/functionalCorrelationRange;
      
      //use following two lines if you would like to take weight = 1/distance
//      structural = 1 - structural;
      //since functional takes values in [-1,1] interval, we map it to [0,1] interval and invert the neg/pos correlation
      //due to hisgh distance meaning negative correlation. (thus, -1 becomes 1, 1 becomes 0 in terms of distance)
      functional = (1.0 - functional)/2.0;
      
      //d_ab = initializerParameters.getParameter(ParameterSet::SPATIAL_DISTANCE) * spatial + 
      //        initializerParameters.getParameter(ParameterSet::STRUCTURAL_CONNECTIVITY) * structural + 
      d_ab = functional/ratio; //initializerParameters.getParameter(ParameterSet::FUNCTIONAL_CONNECTIVITY) * functional;

      labelsDistanceMatrix[edgeIter->getNode1Id()][edgeIter->getNode2Id()] = d_ab;
      labelsDistanceMatrix[edgeIter->getNode2Id()][edgeIter->getNode1Id()] = d_ab;
   }
}
//</editor-fold>
