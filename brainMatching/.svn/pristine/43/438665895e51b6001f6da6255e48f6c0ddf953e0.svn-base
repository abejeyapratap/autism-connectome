/* 
 * File:   linearAssignmentSolver.cpp
 * Author: yusuf
 *
 * Generated on July 7, 2017, 11:05 AM
 */

#include "linearAssignmentSolver.h"
#include <dlib/optimization/max_cost_assignment.h>

using namespace std;
using namespace dlib;

//<editor-fold defaultstate="collapsed" desc=" Constructors">
LinearAssignmentSolver::LinearAssignmentSolver( )
{
	objectsWeightMatrix = NULL;
	assignmentCostMatrix = NULL;
	labelsDistanceMatrix = NULL;
	cerr<<"LinearAssignmentSolver default constructor is called: This should never be called\n";
	exit(1);
}

LinearAssignmentSolver::LinearAssignmentSolver( const Graph & objectGraph, const Graph & labelGraph )
{
	objects = objectGraph;
	labels = labelGraph;

	int numLabels = labels.getNumNodes();
	int numObjects = objects.getNumNodes();

	objectsWeightMatrix = Utility::allocate2Dmemory<float>(numObjects,numObjects);
	labelsDistanceMatrix = Utility::allocate2Dmemory<float>(numLabels,numLabels);
	assignmentCostMatrix = Utility::allocate2Dmemory<float>(numObjects,numLabels);
}

LinearAssignmentSolver::~LinearAssignmentSolver()
{
	int numLabels = labels.getNumNodes();
	int numObjects = objects.getNumNodes();

	Utility::free2Dmemory<float>(objectsWeightMatrix,numObjects);
	Utility::free2Dmemory<float>(assignmentCostMatrix,numObjects);
	Utility::free2Dmemory<float>(labelsDistanceMatrix,numLabels);
}
//</editor-fold>

//<editor-fold defaultstate="collapsed" desc=" Public core functions: doLabeling(), saveMathesToArray()">
float LinearAssignmentSolver::doLabeling(MatrixInitializer &initializer)
{
   float result;
   int sizeOfNodeSet = objects.getNumNodes();
   int sizeOfLabelSet = labels.getNumNodes();
   
   double **pairingProbabilities = Utility::allocate2Dmemory<double>(sizeOfNodeSet,sizeOfLabelSet);
      
   //set assignment and separation costs: 
   //@scaleFactor: parameter for adjusting the contribution of separationCost to the cost function
   initializer.initializeAssignmentCostMatrix(assignmentCostMatrix, objects, labels);
   initializer.initializeObjectsWeightMatrix(objectsWeightMatrix, objects);
   initializer.initializeLabelsDistanceMatrix(labelsDistanceMatrix, labels);

   //DEBUG: print content of the cost matrices
   //cout<<"Assignment cost matrix:\n-----------------------------------------\n";
   //Utility::printMatrix(assignmentCostMatrix,sizeOfNodeSet,sizeOfLabelSet);
   //cout<<"Objects weight matrix:\n-----------------------------------------\n";
   //Utility::printMatrix(objectsWeightMatrix,sizeOfNodeSet);
   //cout<<"labels distance matrix:\n-----------------------------------------\n";
   //Utility::printMatrix(labelsDistanceMatrix,sizeOfLabelSet);

   ////////////////////solve metric labeling
   result = solveOptimizationProblem( pairingProbabilities);

   ///////////////////do rounding
   //metricLabelingRounding(pairingProbabilities);  
   
   //////////////////save matching
//   for(int i=0;i<sizeOfNodeSet;i++)
//   {
//      for(int j=0;j<sizeOfLabelSet;j++)
//      {
//         //uncomment following lines to print the matchings
//         //cout<<pairingProbabilities[i][j]<<"\t";
//         if((int)pairingProbabilities[i][j]==1)
//         {
//             //objects.nodes[i].assignedLabel = labels.nodes[j].label;
//             matches.push_back(pair<int,int>(i,j));
//         }
//      }
//      //cout<<endl;
//      //matches.push_back(pair<int,int>(objects.nodes[i].label,objects.nodes[i].assignedLabel));
//   }
   
   //free allocated memory
   Utility::free2Dmemory<double>(pairingProbabilities,sizeOfNodeSet);

   return result;
}

//given a 2D matchesMatrix of size numOfNodes*numOfNodes that keeps a counter of each pairs of nodes that are matched
//this function increases the counter for each matching as a result of this instance of the metricLabelingSolver
void LinearAssignmentSolver::incrementMatchesMatrix(int** matchesMatrix)
{
   for(std::map<int,int>::iterator iter=matches.begin();iter!=matches.end();iter++)
   {
      matchesMatrix[iter->first][iter->second]++;
      matchesMatrix[iter->second][iter->first]++;
   }
}

void LinearAssignmentSolver::saveMatchesToArray(int* matchesVector)
{
   for(std::map<int,int>::iterator iter=matches.begin();iter!=matches.end();iter++)
      matchesVector[iter->first] = iter->second;
}

void LinearAssignmentSolver::saveMatchesToVector(std::vector<int> &matchesVector)
{
   matchesVector.clear();
   matchesVector.reserve(matches.size());
   for(int i=0;i<matches.size();i++)
      matchesVector.push_back(-1);//fill the matches vector with nonsensical value to detect if anything goes wrong
   for(std::map<int,int>::iterator iter=matches.begin();iter!=matches.end();iter++)
      matchesVector[iter->first] = iter->second;
}

//</editor-fold>

//<editor-fold defaultstate="collapsed" desc=" calculateObjectiveValue()">
float LinearAssignmentSolver::calculateObjectiveValue(std::map<int, int> &_matches)
{
   int numNodes=_matches.size();
   float assignmentCost=0;
   for(std::map<int,int>::iterator iter=_matches.begin();iter!=_matches.end();iter++)
      assignmentCost+=assignmentCostMatrix[iter->first][iter->second];

   return assignmentCost;
}

float LinearAssignmentSolver::calculateIdentityMatchesRatio()
{
   int counter=0;
   for(std::map<int,int>::iterator iter=matches.begin();iter!=matches.end();iter++)
      if(iter->first==iter->second)
         counter++;
   return counter/(float)matches.size();
}
//</editor-fold>

float LinearAssignmentSolver::solveOptimizationProblem( double ** pairingProbabilities )
{
   int numObjects = objects.getNumNodes();
   int numLabels = labels.getNumNodes();
   float objectiveCost = -1;
   
   //Note that, assignment cost matrix is assumed to be non-negative
   //we take (maxValue-assignmentCost) since the max_cost_assignment is a maximization problem
   float maxValue=Utility::getMaxValueOfMatrix<float>(assignmentCostMatrix,numObjects,numLabels);
   float** tempAssignmentCostMatrix = Utility::allocate2Dmemory<float>(numObjects,numLabels);
   Utility::subtractMatrixFromScalar<float>(assignmentCostMatrix,tempAssignmentCostMatrix,maxValue,numObjects,numLabels);
   
   //since max_cost_assignment takes integer values as input, we need to scale the assignment costs to take integer values
   float min = Utility::getMinValueOfMatrixGreaterThanThreshold<float>(tempAssignmentCostMatrix,Utility::EPSILON,numObjects,numLabels);
   float max = Utility::getMaxValueOfMatrix<float>(tempAssignmentCostMatrix,numObjects,numLabels);
   
   float scaleValue=1.0;
   //NOTE: linear assignment solver requires that the numbers in assignment cost matrix is integers. Thus, we need to scale the assignment costs to make them distinguishable for the algorithm. 
   //      Otherwise, roundoff errors will make the assignment of otherwise different nodes equal. 
   //      Proper way of doing this is, we should have the difference between any two assignment costs to be greater than one, unless the difference is equal to zero. 
   //      Nested for loop below makes this in a naive way, which takes quite a time to run for all matchings. Experimentally, I observed that the number is being somewhere close to 10^-9 (7.45058e-09, to be precise)
   //      which is probably due to the precision I used for connectomes while saving them to file. Thus, rather than calculating this number each and every time, here I set 10^9 as the scaling factor
   /*
   float minDiff=1;
   for(int i=0;i<numObjects;i++)
       for(int j=0;j<numLabels;j++)
       {
           if(tempAssignmentCostMatrix[i][j]<=Utility::EPSILON)
               continue;
           for(int p=0;p<numObjects;p++)
               for(int q=0;q<numLabels;q++)
               {
                   if(tempAssignmentCostMatrix[p][q]<=Utility::EPSILON || (i==p && j==q))
                       continue;
                   if(Utility::absoluteValue(tempAssignmentCostMatrix[i][j]-tempAssignmentCostMatrix[p][q])<minDiff && Utility::absoluteValue(tempAssignmentCostMatrix[i][j]-tempAssignmentCostMatrix[p][q])>Utility::EPSILON)
                       minDiff=Utility::absoluteValue(tempAssignmentCostMatrix[i][j]-tempAssignmentCostMatrix[p][q]);
               }
       }
    if(minDiff<1)
    {
       scaleValue = 1.0/minDiff;
       //cout<<minDiff<<"  "<<scaleValue<<"  "<<range/(max*scaleValue)<<endl;
    }
   */
    //we want assignmentCost to take values in [1,max(range,maxValue)]
   float range=10e9;//1000000000.0  //see the explanation above for where this magic number is coming from
           
   //first, make sure that the scaling will correct the min value to be larger than 1
   if(min<1)
      scaleValue = 1.0/min;

   //then make sure that the cost will be scaled to get a range large enough
   if(scaleValue*max<range)
      scaleValue = (range/(max*scaleValue));
   
   Utility::multiplyMatrixWithScalar<float>(tempAssignmentCostMatrix,scaleValue,numObjects,numLabels);

   matrix<int> cost(numObjects,numLabels);
   for(int i=0;i<numObjects;i++)
      for(int j=0;j<numLabels;j++)
         cost(i,j)=Utility::round<float>(tempAssignmentCostMatrix[i][j]);
   
   Utility::free2Dmemory<float>(tempAssignmentCostMatrix,numObjects);
   
   std::vector<long> assignment = max_cost_assignment(cost);

   int assignmentSize = assignment.size();
   for(int i=0;i<assignmentSize;i++)
      matches.insert(std::pair<int,int>(i,assignment[i]));
   
//   objectiveCost = assignment_cost(cost, assignment);
   objectiveCost = calculateObjectiveValue(matches);
   
   return objectiveCost;
}

//<editor-fold defaultstate="collapsed" desc=" Print functions: printMatches(), printAssignmentCostMatrix(), printObjectsWeightMatrix(), printLabelDistanceMatrix()">
void LinearAssignmentSolver::printMatches()
{
   for(map <int, int>::iterator iter = matches.begin(); iter!=matches.end();iter++)
      cout<<iter->first<<"\t"<<iter->second<<endl;
}

void LinearAssignmentSolver::printAssignmentCostMatrix()
{
   Utility::printMatrix(assignmentCostMatrix,objects.getNumNodes(),labels.getNumNodes());
}

void LinearAssignmentSolver::printObjectsWeightMatrix()
{
   Utility::printMatrix(objectsWeightMatrix,objects.getNumNodes());
}

void LinearAssignmentSolver::printLabelsDistanceMatrix()
{
   Utility::printMatrix(labelsDistanceMatrix,labels.getNumNodes());
}

//</editor-fold>
