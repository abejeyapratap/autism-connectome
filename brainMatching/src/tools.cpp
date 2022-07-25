/* 
 * File:   tools.cpp
 * Author: yusuf
 *
 * Generated on May 10, 2017, 1:15 PM
 */

#include "tools.h"
#include "graph.h"
#include <stack>

using namespace std;

//this function applies "threshold consistency" to the structural connectome of a set of subjects
//for which the tracts are calculated using probabilistic tractography.
//@controlsSubjectList: list of subjects according to which we will calculate consistent edges (this should include controls)
//@allSubjectList: list of subjects to whom we will apply the masking of thresholding (this should include controls+patients)
//@streamlineFolder: folder containing the streamline files
//@density: density (in range [0,1]) to calculate the consistent edges
//@outputFolder: folder to save the thresholded connectomes
void Tools::thresholdConsistency(std::string controlsSubjectList,std::string allSubjectList ,std::string streamlineFolder, float density, std::string outputFolder)
{
   if(density>100 || density<=0)
   {
      cerr<<"Density should be either between (0,1] or (0,100] interval!!! Exiting...\n";
      exit(1);
   }
   else if(density>1 && density<=100)
      density/=100.0;
   //<editor-fold defaultstate="collapsed" desc=" calculate thresholding mask from healthy controls ">
   //load control subjects
   vector<string> controlsConnectivityFilenames, controlsSubjectNames;
   
   //load subject names
   Utility::loadFileNamesFromList(controlsSubjectList,controlsSubjectNames);
   int numControlsSubjects = controlsSubjectNames.size();
   
   //load connectivity filenames
   controlsConnectivityFilenames.reserve(numControlsSubjects);
   Utility::loadFilePathFromFolderSelectivelyPreservingOrder(controlsSubjectList,streamlineFolder,controlsConnectivityFilenames);
   
   int numNodes = Utility::countRows(controlsConnectivityFilenames[0]);

   //load connectivity files corresponding to controls into connectivity matrices
   float ***controlsConnectivityMatrices = Utility::allocate3Dmemory<float>(numControlsSubjects,numNodes,numNodes);
   for(int i=0;i<numControlsSubjects;i++)
      Utility::loadMatrixFromFile(controlsConnectivityFilenames[i],controlsConnectivityMatrices[i],numNodes,numNodes);
   
   //allocate space for helper matrices
   float **meanMatrix = Utility::allocate2Dmemory<float>(numNodes,numNodes);
   float **stdDeviationMatrix = Utility::allocate2Dmemory<float>(numNodes,numNodes);
   float **coefficientOfVariationMatrix = Utility::allocate2Dmemory<float>(numNodes,numNodes);
   
   //calculate mean, standard deviation, and coefficient of variation for each node pair across all subjects.
   for(int i=0;i<numNodes;i++)
   {
      for(int j=0;j<numNodes;j++)
      {
         //copy the contents of each connectivity matrix at indexes i,j into a vector of values
         vector<float> values = Utility::copyRowOfMatrixIntoVector(controlsConnectivityMatrices,numControlsSubjects,i,j,0);
         //calculate the mean for the values vector
         meanMatrix[i][j] = Utility::calculateMean(values);
         //calculate the standard deviation for the values vector
         stdDeviationMatrix[i][j] = Utility::calculateStandardDeviation(values,meanMatrix[i][j]);
         //calculate the coefficient of variation(CoV). 
         //Smaller (non-zero) CoV indicates small stdDeviation and large mean, 
         //that is, such a streamline should be kept after the thresholding
         coefficientOfVariationMatrix[i][j] = stdDeviationMatrix[i][j] / (meanMatrix[i][j] + Utility::EPSILON);
      }
   }
   
   //sort the values in the coefficient of variation matrix in ascending order
   vector<Utility::MatrixValueIndex> sorted = Utility::sortValuesOfUpperTriangularMatrixWithIndices<float>(coefficientOfVariationMatrix,numNodes,Utility::MatrixValueIndex::ASCENDING);
   
   //we will make a map that has value 1 for the top density values in the sorted list
   int numToKeep = sorted.size()*density;
   int **maskMatrix = Utility::allocate2Dmemory<int>(numNodes,numNodes);
   Utility::fillMatrix(maskMatrix,0,numNodes,numNodes);
   for(vector<Utility::MatrixValueIndex>::iterator iter=sorted.begin();iter!=sorted.end();iter++)
   {
      if(iter->value > Utility::EPSILON)//we are interested in smallest non-zero entries of the coefficientOfVariationMatrix
      {
         maskMatrix[iter->i][iter->j] = 1;
         numToKeep--;
         
         if(numToKeep<=0)
            break;
      }
   }

   Utility::mirrorUpperTriangleOfTheMatrixToLowerTriangle(maskMatrix,numNodes);
   //</editor-fold>
   
   //<editor-fold defaultstate="collapsed" desc=" apply thresholding mask to all subjects">
   //////now load all subjects and apply masking to the connectomes
   vector<string> allSubjectsConnectivityFilenames, allSubjectNames;
   
   //load subject names
   Utility::loadFileNamesFromList(allSubjectList,allSubjectNames);
   int numAllSubjects = allSubjectNames.size();
   
   //load connectivity filenames
   allSubjectsConnectivityFilenames.reserve(numAllSubjects);
   Utility::loadFilePathFromFolderSelectivelyPreservingOrder(allSubjectList,streamlineFolder,allSubjectsConnectivityFilenames);
   
   //load connectivity files corresponding to controls into connectivity matrices
   float ***allAsubjectsConnectivityMatrices = Utility::allocate3Dmemory<float>(numAllSubjects,numNodes,numNodes);
   for(int i=0;i<numAllSubjects;i++)
      Utility::loadMatrixFromFile(allSubjectsConnectivityFilenames[i],allAsubjectsConnectivityMatrices[i],numNodes,numNodes);
   
   
   //apply the map to each connectivity matrix (this is where actual thresholding takes place)
   for(int i=0;i<numAllSubjects;i++)
      Utility::maskMatrix(allAsubjectsConnectivityMatrices[i],maskMatrix,numNodes,numNodes);
   
   //save the thresholded connectivity matrices into files
   for(int i=0;i<numAllSubjects;i++)
      Utility::saveMatrixToFile(outputFolder+allSubjectNames[i]+".txt",allAsubjectsConnectivityMatrices[i],numNodes,numNodes,' ','\n',8);
   //</editor-fold>
   
   //free the allocated memory
   Utility::free3Dmemory(controlsConnectivityMatrices,numControlsSubjects,numNodes);
   Utility::free3Dmemory(allAsubjectsConnectivityMatrices,numAllSubjects,numNodes);
   Utility::free2Dmemory(meanMatrix,numNodes);
   Utility::free2Dmemory(stdDeviationMatrix,numNodes);
   Utility::free2Dmemory(coefficientOfVariationMatrix,numNodes);
   Utility::free2Dmemory(maskMatrix,numNodes);
}

//this function applies thresholding to the connectomes by removing the edges below a certain threshold that will satisfy 
//a density of @density at the mean connectome among the healthy control set, especially on the probabilistic structural connectomes
//@controlsSubjectList: list of subjects according to which we will calculate consistent edges (this should include controls)
//@allSubjectList: list of subjects to whom we will apply the masking of thresholding (this should include controls+patients)
//@streamlineFolder: folder containing the streamline files
//@density: density (in range [0,1]) to set the edge density in healthy controls
//@outputFolder: folder to save the thresholded connectomes
void Tools::thresholdDensityByGroup(std::string controlsSubjectList,std::string allSubjectList ,std::string streamlineFolder, float density, std::string outputFolder)
{
   if(density>100 || density<=0)
   {
      cerr<<"Density should be either between (0,1] or (0,100] interval!!! Exiting...\n";
      exit(1);
   }
   else if(density>1 && density<=100)
      density/=100.0;
      
   //<editor-fold defaultstate="collapsed" desc=" calculate thresholding mask from healthy controls ">
   //load control subjects
   vector<string> controlsConnectivityFilenames, controlsSubjectNames;

   //load subject names
   Utility::loadFileNamesFromList(controlsSubjectList,controlsSubjectNames);
   int numControlsSubjects = controlsSubjectNames.size();

   //load connectivity filenames
   controlsConnectivityFilenames.reserve(numControlsSubjects);
   Utility::loadFilePathFromFolderSelectivelyPreservingOrder(controlsSubjectList,streamlineFolder,controlsConnectivityFilenames);
   int numNodes = Utility::countRows(controlsConnectivityFilenames[0]);

   //load connectivity files corresponding to controls into connectivity matrices
   float ***controlsConnectivityMatrices = Utility::allocate3Dmemory<float>(numControlsSubjects,numNodes,numNodes);
   for(int i=0;i<numControlsSubjects;i++)
      Utility::loadMatrixFromFile(controlsConnectivityFilenames[i],controlsConnectivityMatrices[i],numNodes,numNodes);

   //allocate space for helper matrices
   float **meanMatrix = Utility::allocate2Dmemory<float>(numNodes,numNodes);
   Utility::fillMatrix<float>(meanMatrix,0,numNodes,numNodes);

   //calculate mean matrix
   for(int i=0;i<numControlsSubjects;i++)
      Utility::elementwiseAddMatrices(meanMatrix,controlsConnectivityMatrices[i],meanMatrix,numNodes,numNodes);

   float scale=1.0/float(numControlsSubjects);
   Utility::multiplyMatrixWithScalar(meanMatrix,scale,numNodes,numNodes);
   
   //////the matrix might include positive and negative values. We would like to threshold the absolute values. Thus we might need two thresholds (positive/negative).
   Utility::absoluteValue(meanMatrix,meanMatrix,numNodes,numNodes);
   
   ////sort elements of the matrix in ascending order, and determine the threshold so as to keep the last portion of the resulting vector based on the desired density
   vector<float> sortedMeanMatrix = Utility::sortValuesOfMatrix(meanMatrix,numNodes,numNodes,Utility::MatrixValueIndex::ASCENDING);
   
   //we will make a map that has value 1 for the top density values in the sorted list
   int numToKeep = numNodes*numNodes*density;
   float threshold=sortedMeanMatrix[sortedMeanMatrix.size() - numToKeep];
   
   //free the allocated memory
   Utility::free3Dmemory(controlsConnectivityMatrices,numControlsSubjects,numNodes);
   Utility::free2Dmemory(meanMatrix,numNodes);
      
   //</editor-fold>
   
   //<editor-fold defaultstate="collapsed" desc=" apply thresholding to all subjects using the calculated threshold, by setting everything less than threshold to zero">
   //////now load all subjects and apply masking to the connectomes
   vector<string> allSubjectsConnectivityFilenames, allSubjectNames;
   
   //load subject names
   Utility::loadFileNamesFromList(allSubjectList,allSubjectNames);
   int numAllSubjects = allSubjectNames.size();
   
   //load connectivity filenames
   allSubjectsConnectivityFilenames.reserve(numAllSubjects);
   Utility::loadFilePathFromFolderSelectivelyPreservingOrder(allSubjectList,streamlineFolder,allSubjectsConnectivityFilenames);
   
   float **connectivityMatrix = Utility::allocate2Dmemory<float>(numNodes,numNodes);
   for(int i=0;i<numAllSubjects;i++)
   {
      //load connectivity files corresponding to controls into connectivity matrices
      Utility::loadMatrixFromFile(allSubjectsConnectivityFilenames[i],connectivityMatrix,numNodes,numNodes);
      //apply the map to each connectivity matrix (this is where actual thresholding takes place)
      Utility::filterOutElementsOfMatrixLessThanThreshold<float>(connectivityMatrix,threshold,0,numNodes,numNodes,true); //true at the end is the parameter to do thresholding with the absolute value
      //save the thresholded connectivity matrices into files
      Utility::saveMatrixToFile(outputFolder+allSubjectNames[i]+".txt",connectivityMatrix,numNodes,numNodes,' ','\n',8);
   }
   //free the allocated memory
   Utility::free2Dmemory(connectivityMatrix,numNodes);
   
   //</editor-fold>
}

//this function applies thresholding to the connectomes by removing the edges below a certain threshold that will satisfy 
//a density of @density per subject, especially on the probabilistic structural connectomes
//@allSubjectList: list of subjects to whom we will apply the masking of thresholding (this should include controls+patients)
//@streamlineFolder: folder containing the streamline files
//@density: density (in range [0,1]) to set the edge density in healthy controls
//@outputFolder: folder to save the thresholded connectomes
void Tools::thresholdDensityPerSubject(std::string allSubjectList ,std::string streamlineFolder, float density, std::string outputFolder)
{
   if(density>100 || density<=0)
   {
      cerr<<"Density should be either between (0,1] or (0,100] interval!!! Exiting...\n";
      exit(1);
   }
   else if(density>1 && density<=100)
      density/=100.0;
      
   vector<string> allSubjectsConnectivityFilenames, allSubjectNames;
   
   //load subject names
   Utility::loadFileNamesFromList(allSubjectList,allSubjectNames);
   int numAllSubjects = allSubjectNames.size();
   
   //load connectivity filenames
   allSubjectsConnectivityFilenames.reserve(numAllSubjects);
   Utility::loadFilePathFromFolderSelectivelyPreservingOrder(allSubjectList,streamlineFolder,allSubjectsConnectivityFilenames);
   
   int numNodes = Utility::countRows(allSubjectsConnectivityFilenames[0]);
   int numToKeep = numNodes*numNodes*density;
   
   float **connectivityMatrix = Utility::allocate2Dmemory<float>(numNodes,numNodes);
   float **connectivityMatrix_absolute = Utility::allocate2Dmemory<float>(numNodes,numNodes);
   for(int i=0;i<numAllSubjects;i++)
   {
      //load connectivity files corresponding to controls into connectivity matrices
      Utility::loadMatrixFromFile(allSubjectsConnectivityFilenames[i],connectivityMatrix,numNodes,numNodes);
      
      //////the matrix might include positive and negative values. We would like to threshold the absolute values. Thus we might need two thresholds (positive/negative).
      Utility::absoluteValue(connectivityMatrix,connectivityMatrix_absolute,numNodes,numNodes);
      
      ////sort elements of the matrix in ascending order, and determine the threshold so as to keep the last portion of the resulting vector based on the desired density
      vector<float> sortedMatrix = Utility::sortValuesOfMatrix(connectivityMatrix_absolute,numNodes,numNodes,Utility::MatrixValueIndex::ASCENDING);
      float threshold=sortedMatrix[sortedMatrix.size() - numToKeep];
      
      //apply the map to each connectivity matrix (this is where actual thresholding takes place)
      Utility::filterOutElementsOfMatrixLessThanThreshold<float>(connectivityMatrix,threshold,0,numNodes,numNodes,true);//true at the end is the parameter to do thresholding with the absolute value
      //save the thresholded connectivity matrices into files
      Utility::saveMatrixToFile(outputFolder+allSubjectNames[i]+".txt",connectivityMatrix,numNodes,numNodes,' ','\n',8);
   }
   //free the allocated memory
   Utility::free2Dmemory(connectivityMatrix,numNodes);
   Utility::free2Dmemory(connectivityMatrix_absolute,numNodes);
}

//this function thresholds the connectomes by setting the edges below a certain threshold to zero and the rest to whatever their value used to be
//@allSubjectList: list of subjects to whom we will apply the masking of threshold (this should include controls+patients)
//@streamlineFolder: folder containing the streamline files
//@outputFolder: folder to save the thresholded connectomes
//@density: density (in range [0,1]) to set the edge density in healthy controls
void Tools::thresholdMatrices(std::string allSubjectList ,std::string streamlineFolder, std::string outputFolder, float threshold)
{
   vector<string> allSubjectsConnectivityFilenames, allSubjectNames;
   
   //load subject names
   Utility::loadFileNamesFromList(allSubjectList,allSubjectNames);
   int numAllSubjects = allSubjectNames.size();
   
   //load connectivity filenames
   allSubjectsConnectivityFilenames.reserve(numAllSubjects);
   Utility::loadFilePathFromFolderSelectivelyPreservingOrder(allSubjectList,streamlineFolder,allSubjectsConnectivityFilenames);
   
   int numNodes = Utility::countRows(allSubjectsConnectivityFilenames[0]);
   
   float **connectivityMatrix = Utility::allocate2Dmemory<float>(numNodes,numNodes);
   for(int i=0;i<numAllSubjects;i++)
   {
      //load connectivity files corresponding to controls into connectivity matrices
      Utility::loadMatrixFromFile(allSubjectsConnectivityFilenames[i],connectivityMatrix,numNodes,numNodes);
      //set values less than threshold to zero
      Utility::filterOutElementsOfMatrixLessThanThreshold<float>(connectivityMatrix,threshold,0,numNodes,numNodes,true);//true at the end is the parameter to do thresholding with the absolute value
      //save the thresholded connectivity matrices into files
      Utility::saveMatrixToFile(outputFolder+allSubjectNames[i]+".txt",connectivityMatrix,numNodes,numNodes,' ','\n',8);
   }
   //free the allocated memory
   Utility::free2Dmemory(connectivityMatrix,numNodes);
}


//this function binarizes the connectomes by setting the edges below a certain threshold to zero and the rest to 1
//@allSubjectList: list of subjects to whom we will apply the masking of binarization (this should include controls+patients)
//@streamlineFolder: folder containing the streamline files
//@outputFolder: folder to save the binarized connectomes
//@density: density (in range [0,1]) to set the edge density in healthy controls
void Tools::binarizeMatrices(std::string allSubjectList ,std::string streamlineFolder, std::string outputFolder, float threshold)
{
   vector<string> allSubjectsConnectivityFilenames, allSubjectNames;
   
   //load subject names
   Utility::loadFileNamesFromList(allSubjectList,allSubjectNames);
   int numAllSubjects = allSubjectNames.size();
   
   //load connectivity filenames
   allSubjectsConnectivityFilenames.reserve(numAllSubjects);
   Utility::loadFilePathFromFolderSelectivelyPreservingOrder(allSubjectList,streamlineFolder,allSubjectsConnectivityFilenames);
   
   int numNodes = Utility::countRows(allSubjectsConnectivityFilenames[0]);
   
   float **connectivityMatrix = Utility::allocate2Dmemory<float>(numNodes,numNodes);
   for(int i=0;i<numAllSubjects;i++)
   {
      //load connectivity files corresponding to controls into connectivity matrices
      Utility::loadMatrixFromFile(allSubjectsConnectivityFilenames[i],connectivityMatrix,numNodes,numNodes);
      
      //set values less than threshold to zero
      Utility::filterOutElementsOfMatrixLessThanThreshold<float>(connectivityMatrix,threshold,0,numNodes,numNodes,true);//true at the end is the parameter to do thresholding with the absolute value
      //set the rest of the values to 1 or -1
      Utility::filterOutElementsOfMatrixGreaterThanThreshold<float>(connectivityMatrix,0,1,numNodes,numNodes,false);//false at the end is the parameter to do thresholding without the absolute value
      Utility::filterOutElementsOfMatrixLessThanThreshold<float>(connectivityMatrix,0,-1,numNodes,numNodes,false);//false at the end is the parameter to do thresholding without the absolute value
      //save the thresholded connectivity matrices into files
      Utility::saveMatrixToFile(outputFolder+allSubjectNames[i]+".txt",connectivityMatrix,numNodes,numNodes,' ','\n',8);
   }
   //free the allocated memory
   Utility::free2Dmemory(connectivityMatrix,numNodes);
}

//given a folder @matricesFolder containing matrices, this function calculates the mean matrix of all and saves in @outFile
//if the @row and @column of the matrices are not provided, function will read the files (assuming that the file only contains the matrix)
// and get the row/column numbers from there
//if @skipLine is not provided, function assumes that the files contain only matrices. If @skipline is provided, then the function
// loads the matrix after skipping this many lines from the beginning of the file.
//sample calls:
//calculateMeanMatrix("folder/","out.txt");
//calculateMeanMatrix("folder/","out.txt",row,col);
//calculateMeanMatrix("folder/","out.txt",row,col,skipLine);
void Tools::calculateMeanMatrix(string matricesFolder, string outFile, bool append, int _row, int _col, int skipline)
{
   vector<string> filenames;
   Utility::loadFilePathFromFolder(matricesFolder,filenames);
   
   int row = (_row == 0 ? Utility::countRowsAfterSkipline(filenames[0],skipline) : _row);
   int col = (_col == 0 ? Utility::countColumnsAfterSkipline(filenames[0],skipline) : _col);
   int numMatrices = filenames.size();
   
   float **meanMatrix = Utility::allocate2Dmemory<float>(row,col);
   float **tempMatrix = Utility::allocate2Dmemory<float>(row,col);
   Utility::fillMatrix<float>(meanMatrix,0,row,col);
   
   for(int i=0;i<numMatrices;i++)
   {
      Utility::loadMatrixFromFileAfterSkipLine<float>(filenames[i],tempMatrix,row,col,skipline);
      Utility::elementwiseAddMatrices<float>(meanMatrix,tempMatrix,meanMatrix,row,col);
   }
   
   float ratio = 1.0/(float)numMatrices;
   Utility::multiplyMatrixWithScalar(meanMatrix,ratio,row,col);
   
   //Utility::fillDiagonalOfTheMatrix<float>(meanMatrix,0.0,row); //uncomment if you would like to set the diagonal entries to 0
   
   if(append==false)
      Utility::saveMatrixToFile(outFile,meanMatrix,row,col,'\t','\n',8);//,'\t','\n',8
   else
      Utility::appendMatrixToFile(outFile,meanMatrix,row,col,'\t','\n',8);//,'\t','\n',8
   
   Utility::free2Dmemory(meanMatrix,row);
   Utility::free2Dmemory(tempMatrix,row);
}

//given a folder @vectorsFolder containing files that has vectors in them that we would like to take average of, 
//this function calculates the mean of all designated vectors and saves it in @outFile
//if the @size of the vector is not provided, function will read the files (assuming that the file only contains the vectors)
// and get the size from there
//if @skipLine is not provided, function assumes that the files contain only vectors. If @skipline is provided, then the function
// loads the vector after skipping this many lines from the beginning of the file.
//sample calls:
//calculateMeanVector("folder/","out.txt");
//calculateMeanVector("folder/","out.txt",size);
//calculateMeanVector("folder/","out.txt",size,skipLine);
void Tools::calculateMeanVector(string vectorsFolder, string outFile, bool append, int _size, int skipline)
{
   vector<string> filenames;
   Utility::loadFilePathFromFolder(vectorsFolder,filenames);

   int size = (_size == 0 ? Utility::countRowsAfterSkipline(filenames[0],skipline) : _size);
   int numVectors = filenames.size();
   
   vector<float> meanVector(size,0), tempVector;
   
   for(int i=0;i<numVectors;i++)
   {
      Utility::loadVectorFromFileAfterSkipLine<float>(filenames[i],tempVector,size,skipline);
      Utility::elementwiseAddVectors<float>(meanVector,tempVector,meanVector);
   }
   
   float ratio = 1.0/(float)numVectors;
   Utility::multiplyVectorWithScalar(meanVector,ratio);
   
   if(append==false)
      Utility::saveVectorToFile(outFile,meanVector);//,'\t','\n',8
   else
      Utility::appendVectorToFile(outFile,meanVector);//,'\t','\n',8
}


//given two files @inFile1 and @inFile2 which contains matrices (assuming that the files contain only matrix),
//this function loads them, subtracts the second from the first elementwise
//and saves the resulting difference matrix into @outFile
void Tools::subtractMatrices(string inFile1, string inFile2, string outFile)
{
   int row = Utility::countRows(inFile1);
   int column = Utility::countColumns(inFile1);
   float **matrix1=Utility::allocate2Dmemory<float>(row,column);
   float **matrix2=Utility::allocate2Dmemory<float>(row,column);
   float **matrix3=Utility::allocate2Dmemory<float>(row,column);
   
   Utility::loadMatrixFromFile<float>(inFile1,matrix1,row,column);
   Utility::loadMatrixFromFile<float>(inFile2,matrix2,row,column);
   
   Utility::elementwiseSubtractMatrices(matrix1,matrix2,matrix3,row,column);
//   Utility::multiplyMatrixWithScalar<float>(matrix3,1000.0,row,column);
   Utility::saveMatrixToFile(outFile,matrix3,row,column,'\t','\n',8);
   
   Utility::printVector<float>(Utility::copyDiagonalOfMatrixIntoVector<float>(matrix1,row),std::cout,'\t','\n',8);
   Utility::printVector<float>(Utility::copyDiagonalOfMatrixIntoVector<float>(matrix2,row),std::cout,'\t','\n',8);
   Utility::printVector<float>(Utility::copyDiagonalOfMatrixIntoVector<float>(matrix3,row),std::cout,'\t','\n',8);
}

//given a folder containing matrices in separate files, this function calculates the
//density of each matrix (i.e., ratio of #nonzero cells to #all cells) and prints statistics of each to the screen
void Tools::calculateDensityOfMatrices(std::string subjectList,std::string matrixFolder)
{
   vector<string> filenames, subjectIDs;
   
   if(subjectList.compare("")==0)
      Utility::loadFilePathFromFolder(matrixFolder,filenames);
   else
      Utility::loadFilePathFromFolderSelectivelyPreservingOrder(subjectList,matrixFolder,filenames);
    
   int numFiles = filenames.size();
   int numNodes = Utility::countColumns(filenames[0]);
   
   float **matrix = Utility::allocate2Dmemory<float>(numNodes,numNodes);
   
   //Check if the matrix contains negative values as well. If so, we will calculate density of positive and negative entries separately
   bool isSigned=false;
   float negativeThreshold=-0.0001;
   Utility::loadMatrixFromFile<float>(filenames[0],matrix,numNodes,numNodes);
   for(int i=0;i<numNodes && isSigned==false;i++)
      for(int j=0;j<numNodes && isSigned==false;j++)
         if(matrix[i][j]<negativeThreshold)
            isSigned=true;
   
   cout.precision(4);
   
   if(isSigned==false)
   {//then, we are dealing with a non-negative matrix, thus calculate the ratio of numbers greater than zero to total number of cells
      vector<float> densities;
      for(std::vector<string>::iterator filenameIter=filenames.begin();filenameIter!=filenames.end();filenameIter++)
      {
         Utility::loadMatrixFromFile<float>(*filenameIter,matrix,numNodes,numNodes);
         float density = Utility::calculateMatrixDensity(matrix,numNodes,"positive");
         densities.push_back(density);
         cout<<Utility::pruneFilenameFromPath(*filenameIter)<<":\t"<<density*100<<"%"<<endl;
      }
      float min,max,mean,median,stdDev;
      max = Utility::getMaxValueOfVector(densities);
      min = Utility::getMinValueOfVector(densities);
      median = Utility::getMedianValueOfVector(densities);
      mean = Utility::calculateMean(densities);
      stdDev = Utility::calculateStandardDeviation(densities,mean);

      cout<<"min : "<<min*100.0<<"%"<<endl;
      cout<<"median : "<<median*100.0<<"%"<<endl;
      cout<<"max : "<<max*100.0<<"%"<<endl;
      cout<<"mean : "<<mean*100.0<<"%"<<endl;
      cout<<"stdDev : "<<stdDev*100.0<<"%"<<endl;
   }
   else
   {//then, we are dealing with a matrix which has positive and negative values. We will calculate densities separately for each
      vector<float> densities[3];
      float density[3];
      float min[3],max[3],mean[3],median[3],stdDev[3];
      std:string types[3]={"positive","negative","absolute"};
      for(std::vector<string>::iterator filenameIter=filenames.begin();filenameIter!=filenames.end();filenameIter++)
      {
         Utility::loadMatrixFromFile<float>(*filenameIter,matrix,numNodes,numNodes);
         for(int i=0;i<3;i++)
         {
            density[i] = Utility::calculateMatrixDensity(matrix,numNodes,types[i]);
            densities[i].push_back(density[i]);
         }

         cout<<Utility::pruneFilenameFromPath(*filenameIter)<<":\t"<<types[0]<<": "<<density[0]*100<<"%\t"<<types[1]<<": "<<density[1]*100<<"%\t"<<types[2]<<": "<<density[2]*100<<endl;
      }
      
      for(int i=0;i<3;i++)
      {
         max[i] = Utility::getMaxValueOfVector(densities[i]);
         min[i] = Utility::getMinValueOfVector(densities[i]);
         median[i] = Utility::getMedianValueOfVector(densities[i]);
         mean[i] = Utility::calculateMean(densities[i]);
         stdDev[i] = Utility::calculateStandardDeviation(densities[i],mean[i]);
         
         cout<<"====="<<types[i]<<"====="<<endl;
         cout<<"min : "<<min[i]*100.0<<"%"<<endl;
         cout<<"median : "<<median[i]*100.0<<"%"<<endl;
         cout<<"max : "<<max[i]*100.0<<"%"<<endl;
         cout<<"mean : "<<mean[i]*100.0<<"%"<<endl;
         cout<<"stdDev : "<<stdDev[i]*100.0<<"%"<<endl;
      }
   }
}

//given a folder of connectomes, and a file containing hemisphere assignment per ROI in conenctomes, this function calculates strength of total,inter,and intra hemispheric connections per connectome
//@hemisphereFile: path to the file that contains hemispehere assignments per ROI as 0/1 per line. Should contain -1 if an ROI should not be counted in either hemisphere(as in cerebellum). If this file is not provided, connectome is assumed to split in half for left/right
void Tools::calculateStrengthOfMatrices(std::string subjectList,std::string matrixFolder,std::string hemisphereFile)
{
   vector<string> filenames, subjectIDs;
   vector<int> hemispheres;
   
   if(subjectList.compare("")==0)
      Utility::loadFilePathFromFolder(matrixFolder,filenames);
   else
      Utility::loadFilePathFromFolderSelectivelyPreservingOrder(subjectList,matrixFolder,filenames);
    
   int numFiles = filenames.size();
   int numNodes = Utility::countColumns(filenames[0]);
   
   if(hemisphereFile.compare("")==0)
   {
       for(int i=0;i<numNodes/2;i++)
           hemispheres.push_back(0);
       for(int i=numNodes/2;i<numNodes;i++)
           hemispheres.push_back(1);
   }
   else
      Utility::loadVectorFromFile(hemisphereFile,hemispheres);
   
   float **matrix = Utility::allocate2Dmemory<float>(numNodes,numNodes);
   
   //Check if the matrix contains negative values as well. If so, we will calculate density of positive and negative entries separately
   std::vector<string> connectivityTypes;
   connectivityTypes.push_back("positive");//we want to check for positive connectivity in either case (i.e., structural or functional connectome))
   
   bool isSigned=false;
   float negativeThreshold=-0.0001;
   Utility::loadMatrixFromFile<float>(filenames[0],matrix,numNodes,numNodes);
   for(int i=0;i<numNodes && isSigned==false;i++)
      for(int j=0;j<numNodes && isSigned==false;j++)
         if(matrix[i][j]<negativeThreshold)
         {
            isSigned=true;
            connectivityTypes.push_back("negative");//we want to check for negative connectivity if there is any negative value in the connectome
         }
   cout.setf(std::ios::fixed,std::ios::floatfield);
   
   
   cout<<"File name\t"<<"Total Connectivity\t"<<"InterHemispheric Connectivity\t"<<"IntraHemispheric Connectivity\n";
   for(std::vector<std::string>::iterator connectivityTypeIter=connectivityTypes.begin();connectivityTypeIter!=connectivityTypes.end();connectivityTypeIter++)
   {
      vector<float> interConnectivities(numFiles,0),intraConnectivities(numFiles,0),totalConnectivities(numFiles,0);
      
      cout<<"===== "<< *connectivityTypeIter <<" connectivities =====\n";
      for(int filenameCounter=0;filenameCounter<numFiles;filenameCounter++)
      {
         Utility::loadMatrixFromFile<float>(filenames[filenameCounter],matrix,numNodes,numNodes);
         totalConnectivities[filenameCounter]=Utility::calculateMatrixConnectivityStrength(matrix,numNodes,hemispheres,"total",*connectivityTypeIter);
         interConnectivities[filenameCounter]=Utility::calculateMatrixConnectivityStrength(matrix,numNodes,hemispheres,"interhemispheric",*connectivityTypeIter);
         intraConnectivities[filenameCounter]=Utility::calculateMatrixConnectivityStrength(matrix,numNodes,hemispheres,"intrahemispheric",*connectivityTypeIter);
         
         cout<<Utility::pruneFilenameFromPath(filenames[filenameCounter])<<":\t"<<totalConnectivities[filenameCounter]<<"\t"<<interConnectivities[filenameCounter]<<"\t"<<intraConnectivities[filenameCounter]<<endl;
      }
      std::vector<float> min(3),max(3),mean(3),median(3),stdDev(3);
      max[0] = Utility::getMaxValueOfVector(totalConnectivities); max[1] = Utility::getMaxValueOfVector(interConnectivities); max[2] = Utility::getMaxValueOfVector(intraConnectivities);
      min[0] = Utility::getMinValueOfVector(totalConnectivities); min[1] = Utility::getMinValueOfVector(interConnectivities); min[2] = Utility::getMinValueOfVector(intraConnectivities);
      median[0] = Utility::getMedianValueOfVector(totalConnectivities); median[1] = Utility::getMedianValueOfVector(interConnectivities); median[2] = Utility::getMedianValueOfVector(intraConnectivities);
      mean[0] = Utility::calculateMean(totalConnectivities); mean[1] = Utility::calculateMean(interConnectivities); mean[2] = Utility::calculateMean(intraConnectivities);
      stdDev[0] = Utility::calculateStandardDeviation(totalConnectivities,mean[0]); stdDev[1] = Utility::calculateStandardDeviation(interConnectivities,mean[1]); stdDev[2] = Utility::calculateStandardDeviation(intraConnectivities,mean[2]);

      cout<<"min : "<<min[0]<<"\t"<<min[1]<<"\t"<<min[2]<<"\t"<<endl;
      cout<<"median : "<<median[0]<<"\t"<<median[1]<<"\t"<<median[2]<<"\t"<<endl;
      cout<<"max : "<<max[0]<<"\t"<<max[1]<<"\t"<<max[2]<<"\t"<<endl;
      cout<<"mean : "<<mean[0]<<"\t"<<mean[1]<<"\t"<<mean[2]<<"\t"<<endl;
      cout<<"stdDev : "<<stdDev[0]<<"\t"<<stdDev[1]<<"\t"<<stdDev[2]<<"\t"<<endl;
   }
}

//given the results of pairwise matchings between structural and functional graphs in the @testMatrix
//and the mapping of brain regions to the functional networks as @trueMapping where the number of networks is @newSize
//this function reduces the oldMatchingMatrix to a matrix of newSize X newSize and saves it into @newMatchingMatrix
void Tools::reduceMatricesIntoFunctionalNetworks(std::string resultsFolder, std::string functionalMappingsFile, std::string outputFolder, bool normalize)
{
   vector<std::string> filePaths,subjectIDs;
   vector<int> functionMaps;
   
   Utility::loadFilePathFromFolder(resultsFolder,filePaths);
   Utility::loadVectorFromFile(functionalMappingsFile,functionMaps);
   
   int numNodes = functionMaps.size();
   int numSystems = Utility::getMaxValueOfVector(functionMaps);
   vector<int> systemSizes(numSystems,0);
   for(vector<int>::iterator mapIter=functionMaps.begin();mapIter!=functionMaps.end();mapIter++)
      systemSizes[*mapIter-1]++;
   
   float **oldMatrix = Utility::allocate2Dmemory<float>(numNodes,numNodes);
   float **newMatrix = Utility::allocate2Dmemory<float>(numSystems,numSystems);
   
   for(vector<std::string>::iterator fileIter=filePaths.begin();fileIter!=filePaths.end();fileIter++)
   {
      Utility::loadMatrixFromFile(*fileIter,oldMatrix,numNodes,numNodes);
      Utility::fillMatrix<float>(newMatrix,0,numSystems,numSystems);
      
      for(int i=0;i<numNodes;i++)
      {
         for(int j=0;j<numNodes;j++)
         {
            int sourceNetwork = functionMaps[i]-1;
            int targetNetwork = functionMaps[j]-1;
            newMatrix[sourceNetwork][targetNetwork]+=oldMatrix[i][j];
         }
      }
      
      if(normalize==true)
      {
         for(int i=0;i<numSystems;i++)
            for(int j=0;j<numSystems;j++)
               newMatrix[i][j] /= (systemSizes[i]+systemSizes[j]);
      }
      
      string filePath = outputFolder+Utility::pruneFilenameFromPath(*fileIter)+".txt";
      Utility::saveMatrixToFile(filePath,newMatrix,numSystems,numSystems);
   }
   
   Utility::free2Dmemory(oldMatrix,numNodes);
   Utility::free2Dmemory(newMatrix,numSystems);
//   Utility::fillMatrix<float>(newMatchingMatrix,0,newSize,newSize);
//   int oldSize=trueMapping.size();
//   for(int i=0;i<oldSize;i++)
//   {
//      for(int j=0;j<oldSize;j++)
//      {
//         int sourceNetwork = trueMapping[i]-1;
//         int targetNetwork = trueMapping[j]-1;
//         newMatchingMatrix[sourceNetwork][targetNetwork]+=oldMatchingMatrix[i][j];
//      }
//   }
}

//given a folder of either graph files or connectivity matrices, 
//this function calculates the number of connected components for each file and prints the count to the stdout.
void Tools::checkGraphConnectivity(std::string connectivityType,std::string graphsFolder,bool printComponents,std::string type,std::string samplesFile)
{
   vector<string> filenames;
   if(connectivityType.compare("structure")!=0 && connectivityType.compare("positiveFunction")!=0 && connectivityType.compare("negativeFunction")!=0)
   {
      cerr<<"wrong feature type entered for Tools::checkGraphConnectivity():"<<connectivityType<<endl;
      cerr<<"Exiting!..."<<endl;
      exit(1);
   }
   if(samplesFile.compare("")==0)
      Utility::loadFilePathFromFolder(graphsFolder,filenames);
   else
      Utility::loadFilePathFromFolderSelectively(samplesFile,graphsFolder,filenames);
   
//   Utility::printVector(filenames,std::cout,'\n','\n');//for DEBUG
   
   bool noDisconnectedGraph=true;//flag for printing "no disconnected graph detected" message, just in case
   
   if(type.compare("graph")==0)
   {
      for(vector<string>::iterator iter = filenames.begin(); iter!=filenames.end();iter++)
      {
         Graph graph(*iter);
         vector<vector<int> > connectedComponents;
         graph.calculateNumberOfConnectedComponents(connectivityType,connectedComponents);
         if(connectedComponents.size()>1)
         {
            noDisconnectedGraph=false;
            cout<<Utility::pruneFilenameFromPath(*iter)<<":"<<connectedComponents.size()-1<<" ( ";
         
            if(printComponents==true)
            {
               //We don't want to print the major connected component. Thus, we get the sizes of all the components, find the largest one, and do not print its components in the section below
               int majorComponentSize=connectedComponents[0].size();
               int majorComponentId=0;
               for(int i=1;i<connectedComponents.size();i++)
               {
                  if(connectedComponents[i].size()>majorComponentSize)
                  {
                     majorComponentSize=connectedComponents[i].size();
                     majorComponentId=i;
                  }
               }
               for(int i=0;i<connectedComponents.size();i++)
               {
                  if(i==majorComponentId)
                     continue;
                  
                  cout<<"<  ";
                  Utility::printVector(connectedComponents[i],std::cout,' ',' ');
                  cout<<">";
               }
            }
            cout<<")\n";
         }
      }
      if(noDisconnectedGraph==true)
         cout<<"All graphs are connected in the given dataset...\n";
   }
   else if(type.compare("matrix")==0)
   {
      int numNodes=Utility::countRows(filenames[0]);
      float **connectivityMatrix=Utility::allocate2Dmemory<float>(numNodes,numNodes);

      for(vector<string>::iterator iter = filenames.begin(); iter!=filenames.end();iter++)
      {
         Utility::loadMatrixFromFile(*iter,connectivityMatrix,numNodes,numNodes);
         Utility::symmetrizeMatrixByMean(connectivityMatrix,numNodes);
         Utility::fillDiagonalOfTheMatrix<float>(connectivityMatrix,0.0,numNodes);
         if(connectivityType.compare("positiveFunction")==0)
         {
            Utility::filterOutElementsOfMatrixLessThanThreshold<float>(connectivityMatrix,0.0,0.0,numNodes,numNodes);
         }
         else if(connectivityType.compare("negativeFunction")==0)
         {
            Utility::filterOutElementsOfMatrixGreaterThanThreshold<float>(connectivityMatrix,0.0,0.0,numNodes,numNodes);
            Utility::multiplyMatrixWithScalar<float>(connectivityMatrix,-1.0,numNodes,numNodes);
         }

         vector<bool> isVisited(numNodes,false);
         vector<vector<int> > connectedComponents;
         int counter=0;
         stack<int> S;

         while(true)
         {
            //pick an unvisited node as the starting point of a DFS (this will constitute a connected component)
            vector<int> component;
            for(int i=0;i<numNodes;i++)
            {
               if(isVisited[i]==false)
               {
                  S.push(i);
                  counter++;
                  break;
               }
            }
            if(S.empty()==true)
               break;

            while(S.empty()==false)
            {
               int node=S.top();
               S.pop();

               if(isVisited[node]==false)
               {
                  isVisited[node]=true;
                  component.push_back(node);
                  for(int i=0;i<numNodes;i++)
                     if(connectivityMatrix[node][i]>0 && isVisited[i]==false)
                        S.push(i);
               }           
            }
            connectedComponents.push_back(component);
         }
         if(connectedComponents.size()>1)
         {
            noDisconnectedGraph=false;
            cout<<Utility::pruneFilenameFromPath(*iter)<<":"<<connectedComponents.size()-1<<" ( ";
         
            if(printComponents==true)
            {
               //We don't want to print the major connected component. Thus, we get the sizes of all the components, find the largest one, and do not print its components in the section below
               int majorComponentSize=connectedComponents[0].size();
               int majorComponentId=0;
               for(int i=1;i<connectedComponents.size();i++)
               {
                  if(connectedComponents[i].size()>majorComponentSize)
                  {
                     majorComponentSize=connectedComponents[i].size();
                     majorComponentId=i;
                  }
               }
               for(int i=0;i<connectedComponents.size();i++)
               {
                  if(i==majorComponentId)
                     continue;
                  
                  cout<<"<  ";
                  Utility::printVector(connectedComponents[i],std::cout,' ',' ');
                  cout<<">";
               }
            }
            cout<<")\n";
         }
      }
      Utility::free2Dmemory<float>(connectivityMatrix,numNodes);
      if(noDisconnectedGraph==true)
         cout<<"All graphs are connected in the given dataset...\n";
   }
   else
   {
      cerr<<"Wrong 'type' entered for Tools::checkGraphConnectivity() function:"<<type<<endl;
      exit(1);
   }
}

void Tools::checkMissingNode(std::string graphsFolder,int numNodes,std::string type,std::string samplesFile)
{
   vector<string> filenames;
   
   if(samplesFile.compare("")==0)
      Utility::loadFilePathFromFolder(graphsFolder,filenames);
   else
      Utility::loadFilePathFromFolderSelectively(samplesFile,graphsFolder,filenames);
   
//   Utility::printVector(filenames,std::cout,'\n','\n');//for DEBUG
   
   bool noMissingNode=true;//flag for printing "no disconnected graph detected" message, just in case
   
   if(type.compare("graph")==0)
   {
      for(vector<string>::iterator iter = filenames.begin(); iter!=filenames.end();iter++)
      {
         Graph graph(*iter);
         if(graph.getNumNodes()!=numNodes)
         {
             noMissingNode=false;
             cout<<Utility::pruneFilenameFromPath(*iter)<<":"<<graph.getNumNodes()<<endl;
         }         
      }
      if(noMissingNode==true)
         cout<<"All graphs have exactly "<<numNodes<<" number of nodes...\n";
   }
   else if(type.compare("matrix")==0)
   {
      for(vector<string>::iterator iter = filenames.begin(); iter!=filenames.end();iter++)
      {
         int tempCount=Utility::countRows(*iter);
         
         if(tempCount!=numNodes)
         {
             noMissingNode=false;
             cout<<Utility::pruneFilenameFromPath(*iter)<<":"<<tempCount<<endl;
         }
      }
      if(noMissingNode==true)
         cout<<"All graphs have exactly "<<numNodes<<" number of nodes...\n";
   }
   else
   {
      cerr<<"Wrong 'type' entered for Tools::checkMissingNode() function:"<<type<<endl;
      exit(1);
   }
}

void Tools::symmetrizeMatrices(std::string matricesFolder, std::string outputFolder)
{
   vector<string> filenames;
   Utility::loadFilePathFromFolder(matricesFolder,filenames);
   
   int numNodes=Utility::countRows(filenames[0]);
   float **connectivityMatrix=Utility::allocate2Dmemory<float>(numNodes,numNodes);

   for(vector<string>::iterator iter = filenames.begin(); iter!=filenames.end();iter++)
   {
      string filename=Utility::pruneFilenameFromPath(*iter);
      cout<<"-------"<<filename<<"-------"<<endl;

      Utility::loadMatrixFromFile(*iter,connectivityMatrix,numNodes,numNodes);
      Utility::symmetrizeMatrixByMean(connectivityMatrix,numNodes);
      Utility::fillDiagonalOfTheMatrix<float>(connectivityMatrix,0.0,numNodes);
      Utility::saveMatrixToFile(outputFolder+filename+".txt",connectivityMatrix,numNodes,numNodes);
   }
   Utility::free2Dmemory<float>(connectivityMatrix,numNodes);
}

//this function takes a set of matrices, where it is assumed that the matrices consists of values >=0
//the goal is to scale up the nonzero values of each matrix such that the smallest nonzero value will become 1
//generated matrices are written into @outputFolder with same filenames from @matricesFolder
void Tools::scaleUpMatricesMakingMinNonZeroEntryUnity(std::string matricesFolder, std::string outputFolder)
{
   vector<string> filenames;
   Utility::loadFilePathFromFolder(matricesFolder,filenames);
   
   int numNodes=Utility::countRows(filenames[0]);
   float **connectivityMatrix=Utility::allocate2Dmemory<float>(numNodes,numNodes);

   for(vector<string>::iterator iter = filenames.begin(); iter!=filenames.end();iter++)
   {
      string filename=Utility::pruneFilenameFromPath(*iter);
      cout<<"-------"<<filename<<"-------"<<endl;

      Utility::loadMatrixFromFile(*iter,connectivityMatrix,numNodes,numNodes);
      float scalar = 1.0/Utility::getNonZeroMinValueOfMatrix(connectivityMatrix,numNodes,numNodes);
      Utility::multiplyMatrixWithScalar(connectivityMatrix,scalar,numNodes,numNodes);
      Utility::saveMatrixToFile(outputFolder+filename+".txt",connectivityMatrix,numNodes,numNodes);
   }
   Utility::free2Dmemory<float>(connectivityMatrix,numNodes);
}

void Tools::saveGraphsAsVectors(Dataset &dataset, std::string outputFile)
{
   ofstream file;
   file.open(outputFile.c_str());
   file.setf(std::ios::fixed,std::ios::floatfield);
   
   for(std::map<int,Graph>::iterator iter=dataset.getGraphIteratorBegin(); iter!=dataset.getGraphIteratorEnd();iter++)
      iter->second.printGraphAsAVector(file);
   
   file.close();
}


//This is a generic function, that is bound to change according to the need
//main idea is, given a folder of matrices, do modifications to each matrix and save them in a target folder
void Tools::modifyMatrices(std::string matricesFolder, std::string outputFolder)
{
   vector<string> filenames;
   Utility::loadFilePathFromFolder(matricesFolder,filenames);
   
   int numNodes=Utility::countRows(filenames[0]);
   float **connectivityMatrix=Utility::allocate2Dmemory<float>(numNodes,numNodes);

   for(vector<string>::iterator iter = filenames.begin(); iter!=filenames.end();iter++)
   {
      string filename=Utility::pruneFilenameFromPath(*iter);
      cout<<"-------"<<filename<<"-------"<<endl;

      Utility::loadMatrixFromFile(*iter,connectivityMatrix,numNodes,numNodes);
      Utility::fillDiagonalOfTheMatrix<float>(connectivityMatrix,0.0,numNodes);
      Utility::saveMatrixToFile(outputFolder+filename+".txt",connectivityMatrix,numNodes,numNodes);
   }
   Utility::free2Dmemory<float>(connectivityMatrix,numNodes);
}