/* 
 * File:   test.cpp
 * Author: yusuf
 *
 * Generated on May 1, 2018, 11:58 AM
 */

#include "test.h"

using namespace std;

void Test::shuffleGraph(string inFile, string outFile, int seedSupplement, int iterationBound, bool isFunctional)
{
   Graph graph(inFile);

//   graph.print("nodes");
   
   //graph.shuffleStructuralConnectivityPreservingStructuralNodeStrength(seedSupplement);
   if(isFunctional==true)
      graph.shuffleFunctionalConnectivityPreservingSignedNodeDegree(seedSupplement,iterationBound);
   else
      graph.shuffleStructuralConnectivityPreservingStructuralNodeDegree(seedSupplement,iterationBound);
   
   graph.saveGraph(outFile);
}

void Test::calculateCovarianceOfAGraph(string inFile)
{
   Graph graph(inFile);
   
   int size=graph.getNumNodes();
   
   vector<float> vec1(size,0), vec2(size,0), vec3(size,0);
   int count=0;
   
   for(vector<Node>::iterator iter=graph.getNodeIteratorBegin();iter!=graph.getNodeIteratorEnd();iter++)
   {
      vec1[count] = iter->getFeature(Node::STR_NODE_STRENGTH);
      vec2[count] = iter->getFeature(Node::FUNC_NODE_STRENGTH_POS);
      vec3[count] = iter->getFeature(Node::FUNC_NODE_STRENGTH_NEG);
      count++;
   }
   
   Utility::printVector(vec1);
   Utility::printVector(vec2);
   Utility::printVector(vec3);
   
   cout<<"Cov str-func+: "<<Utility::calculateCovariance(vec1,vec2)<<"\tCov str-func-: "<<Utility::calculateCovariance(vec1,vec3)<<endl;
   cout<<"Corr str-func+: "<<Utility::calculateCorrelation(vec1,vec2)<<"\tCorr str-func-: "<<Utility::calculateCorrelation(vec1,vec3)<<endl;
}

void Test::testRandom(float min, float max, int randCount, int binCount)
{
   vector<int> binVec(binCount,0);
   
   srand(clock());
   
   float binSize = (max-min)/(float)binCount;
   
   for(int i=0;i<randCount;i++)
   {
      float number = Utility::generateRandomNumber(rand(),min,max);
      int bin = (int)((number-min)/binSize);
      binVec[bin]++;
   }
   Utility::printVector(binVec);
}

void Test::testCommunicability(int numNodes,std::string inputMatrix, std::string referenceMatrix, std::string outputMatrix)
{
   int maxPathLength=20;
   
   float **W=Utility::allocate2Dmemory<float>(numNodes,numNodes);
   Utility::loadMatrixFromFile(inputMatrix,W,numNodes,numNodes);
   
   std::vector<float> nodeStrengths(numNodes,0);
   for(int i=0;i<numNodes;i++)
      for(int j=i+1;j<numNodes;j++)
      {
         nodeStrengths[i]+=W[i][j];
         nodeStrengths[j]+=W[i][j];
      }
   cout<<"\n\nStrength Vector\n";
   Utility::printVector(nodeStrengths);
   
   Graph::calculateWeightedCommunicability(W,nodeStrengths,numNodes,maxPathLength);
   Utility::fillDiagonalOfTheMatrix<float>(W,0,numNodes);
   cout<<"\n\nComm Matrix\n";
   Utility::printMatrix(W,numNodes);
   
   float **Aref=Utility::allocate2Dmemory<float>(numNodes,numNodes);
   Utility::loadMatrixFromFile(referenceMatrix,Aref,numNodes,numNodes);
//   cout<<"\n\nA Reference Matrix\n";
//   Utility::printMatrix(Aref,numNodes);
   
   Utility::elementwiseSubtractMatrices(W,Aref,W,numNodes,numNodes);
   Utility::saveMatrixToFile(outputMatrix,W,numNodes,numNodes,',','\n',6);
   cout<<"\n\nDifference Matrix\n";
   Utility::printMatrix(W,numNodes);
   
}

void Test::testCode()
{
//   float **matrix = Utility::allocate2Dmemory<float>(3,3);
//   
//   for(int power=1;power<5;power++)
//   {
//      matrix[0][0]=0;matrix[0][1]=1;matrix[0][2]=1;
//      matrix[1][0]=1;matrix[1][1]=0;matrix[1][2]=0;
//      matrix[2][0]=1;matrix[2][1]=0;matrix[2][2]=0;
//   
//      Utility::matrixPower<float>(matrix,power,3);
//      cout<<"-------------------------"<<endl;
//      cout<<"Power:"<<power<<endl;
//      Utility::printMatrix<float>(matrix,3);
//   }
//   
//   cout<<"-------------------------"<<endl;
//   matrix[0][0]=0;matrix[0][1]=1;matrix[0][2]=1;
//      matrix[1][0]=1;matrix[1][1]=0;matrix[1][2]=0;
//      matrix[2][0]=1;matrix[2][1]=0;matrix[2][2]=0;
//   Utility::matrixExponent<float>(matrix,20,3);
//   Utility::printMatrix<float>(matrix,3);
//
//   Graph::calculateUnweightedCommunicability(matrix,3,4);
//   
//   Utility::free2Dmemory<float>(matrix,3);
   
//   dlib::matrix<float,3,3> mat;
//   dlib::matrix<float> mat2;
//   mat = 0,1,1,1,0,0,1,0,0;
//   
//   mat2=mat*mat*mat;
//   cout << mat << endl;
//   cout << mat2 << endl;
   
//   cout<<Utility::factorial(6)<<endl;
   
   string test="Hello world!", test2="";
   cout<<test<<"  "<<test.find("H")<<endl;
   cout<<test.compare("")<<endl;
   cout<<test2.compare("")<<endl;
   
//   loadFilePathFromFolderSelectivelyPreservingOrder2(std::string objectListFilename, std::string folderPath, std::vector<std::string> &selectedFilenames)
   
}