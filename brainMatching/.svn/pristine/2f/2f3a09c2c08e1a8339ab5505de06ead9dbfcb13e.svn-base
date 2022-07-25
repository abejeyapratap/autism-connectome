#include "matrixSimilarity.h"
#include "dataset.h"
#include "utility.h"
#include <math.h>

using namespace std;

//<editor-fold defaultstate="collapsed" desc=" calculate distance functions: calculateMatrixDistanceOfEveryoneToItself(), calculateMatrixDistanceBetweenGroups(), subtractEveryoneFromMembersOfGroup(), calculateZscoreOfEdgesForEveryoneWRTGroup(), ">

//calculates the scalar distance (L1, L2, subtract) between a matrix in patients wrt a controls, and saves the result for all patient/control  pairs
//each recorded value is the sum of differences of the values in matrices of a patient and a control
void MatrixSimilarity::calculateMatrixDistanceOfEveryoneToItself(Edge::Feature edgeType1, Edge::Feature edgeType2,std::vector<float>& similarityScores, std::string distanceMeasure)
{
   int numNodes=dataset->getSizeOfAGraph();
   int numSubjects=dataset->getNumOfSubjects();
   float similarityScore;
   similarityScores.clear();
   similarityScores.reserve(numSubjects);
   
   float **connectome1=Utility::allocate2Dmemory<float>(numNodes,numNodes);
   float **connectome2=Utility::allocate2Dmemory<float>(numNodes,numNodes);
   float **differenceMatrix=Utility::allocate2Dmemory<float>(numNodes,numNodes);
   
   // if modalities to be matched are stored in a single dataset, set the pointers to the second dataset to be identical to the first dataset
   // otherwise, set the pointer to second dataset separately
   Dataset *datasetPtr2 = (dataset2==NULL ? dataset : dataset2);
         
   for(int i=0;i<numSubjects;i++)
   {
      dataset->getConnectivityMatrixForGraph(i,edgeType1,connectome1,dataset->getInitializerParameters().functionalConnectivity);
      datasetPtr2->getConnectivityMatrixForGraph(i,edgeType2,connectome2,dataset->getInitializerParameters().functionalConnectivity2);
            
      if(distanceMeasure.compare("subtract")==0)
      {
         Utility::elementwiseSubtractMatrices(connectome1,connectome2,differenceMatrix,numNodes,numNodes);
         similarityScore = Utility::sumElementsOfMatrix(differenceMatrix,numNodes,numNodes);
      }
      else if(distanceMeasure.compare("l1")==0)
         similarityScore = Utility::calculateMatrixDistanceL1(connectome1,connectome2,numNodes,numNodes);
      else if(distanceMeasure.compare("l2")==0)
         similarityScore = Utility::calculateMatrixDistanceL2(connectome1,connectome2,numNodes,numNodes);
      else
      {
         cerr<<"incorrect distance measure:"<<distanceMeasure<<" in Matcher::calculateMatrixDistanceOfEveryoneToItself()...\n";
         exit(1);
      }
      similarityScores.push_back(similarityScore);
   }
   
   Utility::free2Dmemory<float>(connectome1,numNodes);
   Utility::free2Dmemory<float>(connectome2,numNodes);
   Utility::free2Dmemory<float>(differenceMatrix,numNodes);
}


//calculates the scalar distance (L1, L2, subtract) between a matrix in patients wrt a controls, and saves the result for all patient/control  pairs
//each recorded value is the sum of differences of the values in matrices of a patient and a control
void MatrixSimilarity::calculateMatrixDistanceBetweenGroups(float** similarityScores, std::vector<int> &groupRowOrders, std::vector<int> &groupColumnOrders, std::string distanceMeasure, Edge::Feature edgeType1, Edge::Feature edgeType2)
{
   int numOfGroupRow = groupRowOrders.size();
   int numOfGroupColumn = groupColumnOrders.size();
   int numNodes=dataset->getSizeOfAGraph();
   
   float **connectomeRow=Utility::allocate2Dmemory<float>(numNodes,numNodes);
   float **connectomeColumn=Utility::allocate2Dmemory<float>(numNodes,numNodes);
   float **differenceMatrix=Utility::allocate2Dmemory<float>(numNodes,numNodes);
   
   // if modalities to be matched are stored in a single dataset, set the pointers to the second dataset to be identical to the first dataset
   // otherwise, set the pointer to second dataset separately
   Dataset *datasetPtr2 = (dataset2==NULL ? dataset : dataset2);

   //since I am not using dataset wise normalization in this experiment, I do not load the normalizer.
   //initializer.initializeNormalizer(normalizer);//I should have called normalizeDataset() before calling the match function
   for(int i=0;i<numOfGroupRow;i++)
   {
      dataset->getConnectivityMatrixForGraph(groupRowOrders[i],edgeType1,connectomeRow,dataset->getInitializerParameters().functionalConnectivity);
      for(int j=0;j<numOfGroupColumn;j++)
      {
         datasetPtr2->getConnectivityMatrixForGraph(groupColumnOrders[j],edgeType2,connectomeColumn,dataset->getInitializerParameters().functionalConnectivity2);
         if(distanceMeasure.compare("subtract")==0)
         {
            Utility::elementwiseSubtractMatrices(connectomeRow,connectomeColumn,differenceMatrix,numNodes,numNodes);
            similarityScores[i][j] = Utility::sumElementsOfMatrix(differenceMatrix,numNodes,numNodes);
         }
         else if(distanceMeasure.compare("l1")==0)
            similarityScores[i][j] = Utility::calculateMatrixDistanceL1(connectomeRow,connectomeColumn,numNodes,numNodes);
         else if(distanceMeasure.compare("l2")==0)
            similarityScores[i][j] = Utility::calculateMatrixDistanceL2(connectomeRow,connectomeColumn,numNodes,numNodes);
         else
         {
            cerr<<"incorrect distance measure:"<<distanceMeasure<<" in Matcher::calculateMatrixDistanceBetweenGroups()...\n";
            exit(1);
         }
      }
   }
   
   Utility::free2Dmemory<float>(connectomeRow,numNodes);
   Utility::free2Dmemory<float>(connectomeColumn,numNodes);
   Utility::free2Dmemory<float>(differenceMatrix,numNodes);
}

//for each patient, this function subtracts the connectome of the patient from the connectomes of the controls,
//sums the differences, and then normalizes with the number of controls.
//Thus, at the end, we obtain an accumulated difference matrix for each patient wrt the set of controls
//This function differs from calculateMatrixDistanceBetweenGroups() in that, this function produces a matrix per subject, 
//showing the difference of this subject relative to every other subject, 
//whereas the calculateMatrixDistanceBetweenGroups() function returns a single scalar for each subject as a measure of their difference relative to everyone else.
void MatrixSimilarity::subtractEveryoneFromMembersOfGroup(float*** dissimilarityMatrices, std::vector<int> &groupRowOrders, std::vector<int> &groupColumnOrders,std::vector<std::string> &groupRowIDs, std::vector<std::string> &groupColumnIDs, Edge::Feature edgeType1, Edge::Feature edgeType2)
{
   int numOfGroupRow = groupRowIDs.size();
   int numOfGroupColumn = groupColumnIDs.size();
   int numNodes=dataset->getSizeOfAGraph();
   
   float **connectomeRow=Utility::allocate2Dmemory<float>(numNodes,numNodes);
   float **connectomeColumn=Utility::allocate2Dmemory<float>(numNodes,numNodes);
   float **differenceMatrix=Utility::allocate2Dmemory<float>(numNodes,numNodes);
   float **accumulateDifferenceMatrix=Utility::allocate2Dmemory<float>(numNodes,numNodes);
   
   // if modalities to be matched are stored in a single dataset, set the pointers to the second dataset to be identical to the first dataset
   // otherwise, set the pointer to second dataset separately
   Dataset *datasetPtr2 = (dataset2==NULL ? dataset : dataset2);
   
   //since I am not using dataset wise normalization in this experiment, I do not load the normalizer.
   //initializer.initializeNormalizer(normalizer);//I should have called normalizeDataset() before calling the match function
   for(int i=0;i<numOfGroupRow;i++)
   {
      dataset->getConnectivityMatrixForGraph(groupRowOrders[i],edgeType1,connectomeRow,dataset->getInitializerParameters().functionalConnectivity,false);
      Utility::fillMatrix<float>(accumulateDifferenceMatrix,0,numNodes,numNodes);
      float scalar=0;
      for(int j=0;j<numOfGroupColumn;j++)
      {
         if(groupRowIDs[i].compare(groupColumnIDs[j])==0)
            continue;
         datasetPtr2->getConnectivityMatrixForGraph(groupColumnOrders[j],edgeType2,connectomeColumn,dataset->getInitializerParameters().functionalConnectivity2,false);

         Utility::elementwiseSubtractMatrices(connectomeRow,connectomeColumn,differenceMatrix,numNodes,numNodes);
         Utility::elementwiseAddMatrices(differenceMatrix,accumulateDifferenceMatrix,accumulateDifferenceMatrix,numNodes,numNodes);
         scalar+=1;
      }
      Utility::multiplyMatrixWithScalar<float>(accumulateDifferenceMatrix,1.0/scalar,numNodes,numNodes);
      Utility::copyMatrix(accumulateDifferenceMatrix,dissimilarityMatrices[i],numNodes,numNodes);
   }
   
   Utility::free2Dmemory<float>(connectomeRow,numNodes);
   Utility::free2Dmemory<float>(connectomeColumn,numNodes);
   Utility::free2Dmemory<float>(differenceMatrix,numNodes);
   Utility::free2Dmemory<float>(accumulateDifferenceMatrix,numNodes);
}

//for each edge of each patient, this function calculates the zScore of the edge weight wrt the connectomes of the set of controls
//@groupRow : is the control group, with respect to whom the z score is being calculated
//@groupColumn : is the rest of the subjects whose scores are being calculated
void MatrixSimilarity::calculateZscoreOfEdgesForEveryoneWRTGroup(float*** zScoreMatrices, std::vector<int> &groupRowOrders, std::vector<int> &groupColumnOrders,std::vector<std::string> &groupRowIDs, std::vector<std::string> &groupColumnIDs, Edge::Feature edgeType1, Edge::Feature edgeType2)
{
   int numOfGroupRow = groupRowIDs.size();
   int numOfGroupColumn = groupColumnIDs.size();
   int numNodes=dataset->getSizeOfAGraph();
   
   float **connectomeRow=Utility::allocate2Dmemory<float>(numNodes,numNodes);
   float ***connectomesColumn=Utility::allocate3Dmemory<float>(numOfGroupColumn,numNodes,numNodes);
   
   // if modalities to be matched are stored in a single dataset, set the pointers to the second dataset to be identical to the first dataset
   // otherwise, set the pointer to second dataset separately
   Dataset *datasetPtr2 = (dataset2==NULL ? dataset : dataset2);
   
   for(int i=0;i<numOfGroupColumn;i++)
      datasetPtr2->getConnectivityMatrixForGraph(groupColumnOrders[i],edgeType2,connectomesColumn[i],dataset->getInitializerParameters().functionalConnectivity2,false);
   
   //since I am not using dataset wise normalization in this experiment, I do not load the normalizer.
   //initializer.initializeNormalizer(normalizer);//I should have called normalizeDataset() before calling the match function
   for(int subRow=0;subRow<numOfGroupRow;subRow++)
   {
      dataset->getConnectivityMatrixForGraph(groupRowOrders[subRow],edgeType1,connectomeRow,dataset->getInitializerParameters().functionalConnectivity,false);
      
      for(int i=0;i<numNodes;i++)
      {
         for(int j=0;j<numNodes;j++)
         {
            vector<float> distribution;
            for(int subCol=0;subCol<numOfGroupColumn;subCol++)
               distribution.push_back(connectomesColumn[subCol][i][j]);
            
            zScoreMatrices[subRow][i][j] = Utility::calculateZScore(connectomeRow[i][j],distribution);
         }
      }
   }
   
   Utility::free2Dmemory<float>(connectomeRow,numNodes);
   Utility::free3Dmemory<float>(connectomesColumn,numOfGroupColumn,numNodes);
}
//</editor-fold>

//<editor-fold defaultstate="collapsed" desc=" calculate correlation functions: correlateEveryoneWithItself(), correlateGroups() ">

//calculates the Pearson's correlation between the upper triangulars of the two connectomes of the same subject (possibly the correlation between structure and function)
//and saves the 
void MatrixSimilarity::correlateEveryoneWithItself(std::vector<float>& rValues, Edge::Feature edgeType1, Edge::Feature edgeType2)
{
   int numNodes = dataset->getSizeOfAGraph();
   int numSubjects = dataset->getNumOfSubjects();
   rValues.clear();
   rValues.reserve(numSubjects);
   
   float **connectome1=Utility::allocate2Dmemory<float>(numNodes,numNodes);
   float **connectome2=Utility::allocate2Dmemory<float>(numNodes,numNodes);
   
   // if modalities to be matched are stored in a single dataset, set the pointers to the second dataset to be identical to the first dataset
   // otherwise, set the pointer to second dataset separately
   Dataset *datasetPtr2 = (dataset2==NULL ? dataset : dataset2);
   
   for(int i=0;i<numSubjects;i++)
   {
      dataset->getConnectivityMatrixForGraph(i,edgeType1,connectome1,dataset->getInitializerParameters().functionalConnectivity);
      datasetPtr2->getConnectivityMatrixForGraph(i,edgeType2,connectome2,dataset->getInitializerParameters().functionalConnectivity2);
      
      vector<float> conn1 = Utility::copyOffDiagonalUpperTriangularOfMatrixIntoVector<float>(connectome1,numNodes);
      vector<float> conn2 = Utility::copyOffDiagonalUpperTriangularOfMatrixIntoVector<float>(connectome2,numNodes);
      
      float r = Utility::calculateCorrelation(conn1,conn2);
      
      //NOTE: p-value calculation part is missing. Thus, curently, I'm asuming that, the correlation value will be a valid measure of similarity regardless of p-value being significant or not...
      //If you change your mind in the future and would like to only consider significant correlations, then see below!
      // Here are some links: To calculate t-distribution https://jamesmccaffrey.wordpress.com/2016/04/27/implementing-the-students-t-distribution-density-function-in-code/ 
      // To install boost library: https://www.boost.org/doc/libs/1_65_1/libs/math/doc/html/math_toolkit/stat_tut/weg/st_eg/tut_mean_test.html
      // basics of statistics of p-value calculation: https://www.statisticshowto.datasciencecentral.com/probability-and-statistics/correlation-coefficient-formula/
      // and the t-table: https://www.statisticshowto.datasciencecentral.com/tables/ppmc-critical-values/

      rValues.push_back(r);
   }
   Utility::free2Dmemory<float>(connectome1,numNodes);
   Utility::free2Dmemory<float>(connectome2,numNodes);
}

//correlates each subject from the first group with each subject from the second group
//then saves the the matching nodes for each matching that is performed
//@group1IDs: order number of the subjects that are in the first group. Let group1IDs.size()=m
//@group2IDs: order number of the subjects that are in the second group. Let group1IDs.size()=n
//@rValues: a 2D matrix that stores the similarity scores between the matched graphs. similarityScores.size()=m.n
//@edgeType1/@edgeType2: which type of data is going to be used in matching ('struture','function','strFunc')
void MatrixSimilarity::correlateGroups(float** rValues, std::vector<int> &group1Orders, std::vector<int> &group2Orders, Edge::Feature edgeType1, Edge::Feature edgeType2)
{
   int numNodes = dataset->getSizeOfAGraph();
   int numOfGroup1 = group1Orders.size();
   int numOfGroup2 = group2Orders.size();
   int taskSize = numOfGroup1*numOfGroup2;
   int taskCount=0;
   
   float **connectome1=Utility::allocate2Dmemory<float>(numNodes,numNodes);
   float **connectome2=Utility::allocate2Dmemory<float>(numNodes,numNodes);
   
   // if modalities to be matched are stored in a single dataset, set the pointers to the second dataset to be identical to the first dataset
   // otherwise, set the pointer to second dataset separately
   Dataset *datasetPtr2 = (dataset2==NULL ? dataset : dataset2);
   
   for(int i=0;i<numOfGroup1;i++)
   {
      for(int j=0;j<numOfGroup2;j++)
      {
         dataset->getConnectivityMatrixForGraph(group1Orders[i],edgeType1,connectome1,dataset->getInitializerParameters().functionalConnectivity);
         datasetPtr2->getConnectivityMatrixForGraph(group2Orders[j],edgeType2,connectome2,dataset->getInitializerParameters().functionalConnectivity2);

         vector<float> conn1 = Utility::copyOffDiagonalUpperTriangularOfMatrixIntoVector<float>(connectome1,numNodes);
         vector<float> conn2 = Utility::copyOffDiagonalUpperTriangularOfMatrixIntoVector<float>(connectome2,numNodes);

         rValues[i][j] = Utility::calculateCorrelation(conn1,conn2);

         //NOTE: p-value calculation part is missing. Thus, curently, I'm asuming that, the correlation value will be a valid measure of similarity regardless of p-value being significant or not...
         //If you change your mind in the future and would like to only consider significant correlations, then see below!
         // Here are some links: To calculate t-distribution https://jamesmccaffrey.wordpress.com/2016/04/27/implementing-the-students-t-distribution-density-function-in-code/ 
         // To install boost library: https://www.boost.org/doc/libs/1_65_1/libs/math/doc/html/math_toolkit/stat_tut/weg/st_eg/tut_mean_test.html
         // basics of statistics of p-value calculation: https://www.statisticshowto.datasciencecentral.com/probability-and-statistics/correlation-coefficient-formula/
         // and the t-table: https://www.statisticshowto.datasciencecentral.com/tables/ppmc-critical-values/
         
         taskCount++;
      }
      //cerr<<(100.0*taskCount)/(float)taskSize<<"% of the job finished: row "<<i<<endl;
   }
   Utility::free2Dmemory<float>(connectome1,numNodes);
   Utility::free2Dmemory<float>(connectome2,numNodes);
}
//</editor-fold>