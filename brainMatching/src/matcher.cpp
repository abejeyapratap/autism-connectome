/* 
 * File:   matcher.cpp
 * Author: yusuf
 *
 * Generated on September 15, 2016, 4:44 PM
 */

#include "matcher.h"
#include "dataset.h"
#include "utility.h"
#include "linearAssignmentSolver.h"
#include "matrixInitializer.h"
#include "strFuncCoupling_MI.h"
#include "functionCoupling_MI.h"
#include "structureCoupling_MI.h"
#include <math.h>

using namespace std;

//<editor-fold defaultstate="collapsed" desc=" Match functions: matchEveryoneWithItself(), matchGroups()">

//performs matching over the entire dataset subjectwise (i.e., matches each subject with itself and 
//saves the similarity scores to the provided @similarityScores vector
//The graphs to be matched can be located under the same dataset, or can be located in separate datasets
//@similarityScores is of size numOfSubjects
//@matchesMatrix is of size numOfSubjects X numOfNodes, keeps the account of which nodes is matched to which node
void Matcher::matchEveryoneWithItself(vector<float> &similarityScores, int **matchesMatrix, Edge::Feature edgeType1, Edge::Feature edgeType2)
{
   int numGraphs = dataset->getNumOfSubjects();
   int taskSize = numGraphs;
   int taskCount=0;
   
   MatrixInitializer *initializer = allocateMatrixInitializer(edgeType1, edgeType2);//we use this object for making use of its initializeAssignmentCostMatrix/initializeObjectsWeightMatrix/initializeLabelsDistanceMatrix functions

   // if modalities to be matched are stored in a single dataset, set the pointers to the second dataset to be identical to the first dataset
   // otherwise, set the pointer to second dataset separately
   Dataset *datasetPtr2 = (dataset2==NULL ? dataset : dataset2);
   
   for(int i=0;i<numGraphs;i++)
   {
      LinearAssignmentSolver * ml =  new LinearAssignmentSolver(dataset->getGraph(i), datasetPtr2->getGraph(i));

      similarityScores[i] = ml->doLabeling(*initializer);

      ml->saveMatchesToArray(matchesMatrix[i]);

      delete ml;
      //taskCount++;
      //cerr<<(100.0*taskCount)/(float)taskSize<<"% of the job finished: row "<<iter1->first<<endl;
   }
   delete initializer;
}

//matches each subject from the first group with each subject from the second group
//then saves the the matching nodes for each matching that is performed
//@group1IDs: order number of the subjects that are in the first group. Let group1IDs.size()=m
//@group2IDs: order number of the subjects that are in the second group. Let group1IDs.size()=n
//@similarityScores: a 2D matrix that stores the similarity scores between the matched graphs. similarityScores.size()=m.n
//@matchesMatrix: a 3D matrix that keeps the account of which nodes is matched to which node. 
//                matchesMatrix is of size mxnxnumOfNodes. Each row of the matrix keeps the node ids of the mapping nodes
//@edgeType1/@edgeType2: which type of data is going to be used in matching ('struture','function','strFunc')
void Matcher::matchGroups(float** similarityScores, int*** matchesMatrix, std::vector<int> &group1IDs, std::vector<int> &group2IDs, Edge::Feature edgeType1, Edge::Feature edgeType2)
{

	int numOfGroup1 = group1IDs.size();
	int numOfGroup2 = group2IDs.size();
	int taskSize = numOfGroup1*numOfGroup2;
	int taskCount=0;

	MatrixInitializer *initializer = allocateMatrixInitializer(edgeType1, edgeType2);//we use this object for making use of its initializeAssignmentCostMatrix/initializeObjectsWeightMatrix/initializeLabelsDistanceMatrix functions

	// if modalities to be matched are stored in a single dataset, set the pointers to the second dataset to be identical to the first dataset
	// otherwise, set the pointer to second dataset separately
	Dataset *datasetPtr2 = (dataset2==NULL ? dataset : dataset2);

	for(int i=0;i<numOfGroup1;i++)
	{
		for(int j=0;j<numOfGroup2;j++)
		{//cout<<group1IDs[i]<<"--"<<group1IDs[j]<<endl;
			LinearAssignmentSolver * ml =  new LinearAssignmentSolver(dataset->getGraph(group1IDs[i]), datasetPtr2->getGraph(group2IDs[j]));
			similarityScores[i][j] = ml->doLabeling(*initializer);

			ml->saveMatchesToArray(matchesMatrix[i][j]);

			delete ml;
			taskCount++;
		}
		//cerr<<(100.0*taskCount)/(float)taskSize<<"% of the job finished: row "<<i<<endl;
	}   
	delete initializer;
}

//</editor-fold>

//<editor-fold defaultstate="collapsed" desc=" Utility functions: saveResults(), allocateMatrixInitializer() ">

//this function initializes the MatrixInitializer objects, depending on what type of matching will be done
MatrixInitializer * Matcher::allocateMatrixInitializer(Edge::Feature edgeType1, Edge::Feature edgeType2)
{
   MatrixInitializer *initializer;
   if(edgeType1==Edge::FUNCTIONAL_CONNECTIVITY && edgeType2==Edge::FUNCTIONAL_CONNECTIVITY)
     initializer = new FunctionCoupling_MI(dataset->getInitializerParameters());
   else if(edgeType1==Edge::STRUCTURAL_CONNECTIVITY && edgeType2==Edge::STRUCTURAL_CONNECTIVITY)
     initializer = new StructureCoupling_MI(dataset->getInitializerParameters());
   else if(edgeType1==Edge::STRUCTURAL_CONNECTIVITY && edgeType2==Edge::FUNCTIONAL_CONNECTIVITY)
     initializer = new StrFuncCoupling_MI(dataset->getInitializerParameters());
   else
   {
      cerr<<"incorrect modality entered for matching. Error detected in allocateMatrixInitializer(). Exiting!!\n";
      exit(1);
   }
   return initializer;
}


void Matcher::saveResults(string path, int taskSize, TableRow* results, bool saveMatches)
{
   ofstream outFile;
   outFile.open(path.c_str(), ios::out | ios::app);
   outFile.precision(3);	
   outFile.setf(std::ios::fixed,std::ios::floatfield);
   for(int i=0;i<taskSize;i++)
   {
      if(results[i].obj1==results[i].obj2)
      {
         //if the objects to be compared are the same, we do not keep any results for this case. Simply skip
         continue;
      }
      outFile<<results[i].obj1<<"\t"<<results[i].obj2<<"\t"<<results[i].timeInit<<"\t"<<results[i].timeRun<<"\t"<<results[i].objValue<<"\t";
      if(saveMatches==true)
         Utility::printVector<int>(results[i].matches,outFile);
      else
         outFile<<endl;
   }
   outFile.close();
}

//</editor-fold>

//<editor-fold defaultstate="collapsed" desc=" Match functions from leaveOneOut experiment... Most of the functions in here are deprecated.">

//given a range of rows and columns in the dataset, performs a matching and saves the results 
//to the designated file
//NOTE: you can use matchGroups() function with gorup labels "all" and "all" to perform the same matching as being done in this experiment
void Matcher::matchEveryoneWithEveryone(std::string outputFile, int rowStart, int rowEnd, int columnStart, int columnEnd)
{
   int rowCount = rowEnd-rowStart+1;
   int columnCount = columnEnd-columnStart+1;
   int taskSize = rowCount*columnCount;
   int taskCount=0;
   StructureCoupling_MI initializer(dataset->getInitializerParameters());//we use this object for making use of its initializeAssignmentCostMatrix/initializeObjectsWeightMatrix/initializeLabelsDistanceMatrix functions
   initializer.initializeNormalizer(dataset->getNormalizer());//I should have called normalizeDataset() before calling the match function
   
   long time[3];
   TableRow *results = new (nothrow) TableRow[columnCount];
   
   for(int i=rowStart;i<=rowEnd;i++)
   {
      int item=0;
      for(int j=columnStart;j<=columnEnd;j++)
      {
         results[item].obj1 = i;
         results[item].obj2 = j;
         //cout<<i<<"\t"<<j<<"\t";
         if(i!=j)
         {
            time[0] = Utility::getTime();
            LinearAssignmentSolver * ml =  new LinearAssignmentSolver(dataset->getGraph(i), dataset->getGraph(j));
            time[1]  =Utility::getTime();
            results[item].timeInit = time[1] - time[0];
            
            results[item].objValue = ml->doLabeling(initializer);
            time[2] = Utility::getTime();
            results[item].timeRun = time[2] - time[1];
            
            delete ml;
         }
         else
         {
            results[item].timeInit = 0;
            results[item].objValue = 0;
            results[item].timeRun = 0;
         }
         item++;
         taskCount++;
      }
      //cerr<<(100.0*taskCount)/(float)taskSize<<"% of the job finished: row "<<i<<endl;
      saveResults(outputFile,item,results);
   }
   delete[] results;
}

//performs matching over the entire dataset and saves the results
//to the designated file
//NOTE: you can use matchGroups() function with gorup labels "all" and "all" to perform the same matching as being done in this experiment
void Matcher::matchEveryoneWithEveryone(std::string outputFile, bool saveMatches)
{
   int numGraphs = dataset->getNumOfSubjects();
   int rowCount = numGraphs;
   int columnCount = numGraphs;
   int taskSize = rowCount*columnCount;
   int taskCount=0;
   
   StructureCoupling_MI initializer(dataset->getInitializerParameters());//we use this object for making use of its initializeAssignmentCostMatrix/initializeObjectsWeightMatrix/initializeLabelsDistanceMatrix functions
   //initializer.initializeNormalizer(normalizer);//I should have called normalizeDataset() before calling the match function
   
   long time[3];
   TableRow *results = new (nothrow) TableRow[columnCount];//we save results after matching each graph with all other graphs (i.e., all columns)
   
   for(map<int,Graph>::iterator iter1=dataset->getGraphIteratorBegin();iter1!=dataset->getGraphIteratorEnd();iter1++)
   {
      int item=0;
      for(map<int,Graph>::iterator iter2=dataset->getGraphIteratorBegin();iter2!=dataset->getGraphIteratorEnd();iter2++)
      {
         results[item].obj1 = iter1->first;
         results[item].obj2 = iter2->first;
         
         if(iter1->first!=iter2->first)
         {
            time[0] = Utility::getTime();
            LinearAssignmentSolver * ml = new LinearAssignmentSolver(iter1->second, iter2->second);
            time[1]  =Utility::getTime();
            results[item].timeInit = time[1] - time[0];
            
            results[item].objValue = ml->doLabeling(initializer);
            time[2] = Utility::getTime();
            results[item].timeRun = time[2] - time[1];
            
            ml->saveMatchesToVector(results[item].matches);
            
            delete ml;
         }
         else
         {
            results[item].timeInit = 0;
            results[item].objValue = 0;
            results[item].timeRun = 0;
            results[item].matches.clear();
         }
         item++;
         taskCount++;
      }
      //cerr<<(100.0*taskCount)/(float)taskSize<<"% of the job finished: row "<<iter1->first<<endl;
      saveResults(outputFile,item,results,saveMatches);
   }
   delete[] results;
}

//performs matching over the entire dataset and saves the results
//to the given float matrix (we do not save the timing info in this case)
//matches matrix keeps the counter of which nodes is matched to which node
//resultMatrix is of size numOfObjects X numOfObjects
//matchesMatrix is of size numOfNodes X numOfNodes
//NOTE: you can use matchGroups() function with gorup labels "all" and "all" to perform the same matching as being done in this experiment
void Matcher::matchEveryoneWithEveryone(float **resultMatrix, int ***matchesMatrix)
{
   int numGraphs = dataset->getNumOfSubjects();
   int rowCount = numGraphs;
   int columnCount = numGraphs;
   int taskSize = rowCount*columnCount;
   int taskCount=0;
   int numNodes=dataset->getSizeOfAGraph();
   
   StructureCoupling_MI initializer(dataset->getInitializerParameters());//we use this object for making use of its initializeAssignmentCostMatrix/initializeObjectsWeightMatrix/initializeLabelsDistanceMatrix functions
   initializer.initializeNormalizer(dataset->getNormalizer());//I should have called normalizeDataset() before calling the match function
   
   int count1=0;
   for(map<int,Graph>::iterator iter1=dataset->getGraphIteratorBegin();iter1!=dataset->getGraphIteratorEnd();iter1++)
   {
      int count2=0;

      for(map<int,Graph>::iterator iter2=dataset->getGraphIteratorBegin();iter2!=dataset->getGraphIteratorEnd();iter2++)
      {
         if(iter1->first!=iter2->first)
         {
			 LinearAssignmentSolver * ml =  new LinearAssignmentSolver(iter1->second, iter2->second);

            resultMatrix[iter1->first][iter2->first] = ml->doLabeling(initializer);

            ml->saveMatchesToArray(matchesMatrix[iter1->first][iter2->first]);

            delete ml;
         }
         else
         {
            resultMatrix[iter1->first][iter2->first] = 0;
            Utility::fillVector(matchesMatrix[iter1->first][iter2->first],-1,numNodes);
         }
         count2++;
      }
      count1++;
      //cerr<<(100.0*taskCount)/(float)taskSize<<"% of the job finished: row "<<iter1->first<<endl;
   }
}

void Matcher::matchOneWithEveryone(int testGraphId, float** resultMatrix, int** matchesMatrix)
{
   int numGraphs = dataset->getNumOfSubjects();
   int rowCount = numGraphs;
   int columnCount = numGraphs;
   int taskSize = columnCount-1;
   int taskCount=0;
   
   StructureCoupling_MI initializer(dataset->getInitializerParameters());//we use this object for making use of its initializeAssignmentCostMatrix/initializeObjectsWeightMatrix/initializeLabelsDistanceMatrix functions
   initializer.initializeNormalizer(dataset->getNormalizer());//I should have called normalizeDataset() before calling the match function
   
   for(map<int,Graph>::iterator iter1=dataset->getGraphIteratorBegin();iter1!=dataset->getGraphIteratorEnd();iter1++)
   {
      if(iter1->first!=testGraphId)
      {
		 LinearAssignmentSolver * ml1 =  new LinearAssignmentSolver(iter1->second, dataset->getGraph(testGraphId));
		 LinearAssignmentSolver * ml2 =  new LinearAssignmentSolver(dataset->getGraph(testGraphId), iter1->second);


         resultMatrix[iter1->first][testGraphId] = ml1->doLabeling(initializer);
         resultMatrix[testGraphId][iter1->first] = ml2->doLabeling(initializer);

         if(matchesMatrix!=NULL)
         {
            ml1->incrementMatchesMatrix(matchesMatrix);
            ml2->incrementMatchesMatrix(matchesMatrix);
         }

         delete ml1;
         delete ml2;
      }
      else
      {
         resultMatrix[iter1->first][testGraphId] = 0;
      }
      taskCount++;
      //cerr<<(100.0*taskCount)/(float)taskSize<<"% of the job finished: row "<<iter1->first<<endl;
   }
}

void Matcher::matchTwoGraphs(int g1, int g2,bool printMatches)
{
	long time[8];
	int timeInit, timeRun;
	float objValue;

	StructureCoupling_MI initializer(dataset->getInitializerParameters());

	time[0] = Utility::getTime();
	LinearAssignmentSolver * ml =  new LinearAssignmentSolver(dataset->getGraph(g1), dataset->getGraph(g2));

	time[1]  =Utility::getTime();
	timeInit = time[1] - time[0];

	objValue = ml->doLabeling(initializer);
	time[2] = Utility::getTime();
	timeRun = time[2] - time[1];

	if(printMatches==true)
	   ml->printMatches();
	float identityMatches=ml->calculateIdentityMatchesRatio();

	cout<<"Time Init:"<<timeInit<<" TimeRun:"<<timeRun<<" ObjValue:"<<objValue<<" identityMatches:"<<identityMatches<<endl;

	delete ml; 
}
//</editor-fold>
