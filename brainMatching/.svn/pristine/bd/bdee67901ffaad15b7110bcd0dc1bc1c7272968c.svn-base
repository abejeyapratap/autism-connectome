/* 
 * File:   groupwiseExperiment.cpp
 * Author: yusuf
 *
 * Generated on February 22, 2018, 3:58 PM
 */


#include <complex>
#include <math.h>
#include "groupwiseExperiment.h"
#include "matcher.h"
#include "matrixSimilarity.h"
#include "dataset.h"
#include "utility.h"

using namespace std;

//given a set of patient and controls subjects and a modality (structure/function)
//this function runs a matching between the connectomes of the same modality among one group (such as all or patients) 
//compared to another group of subjects (such as healthy controls or all subjects)
//===inputs===
//@modality: <structure | function>
//@numPermutation: number of times to repeat the experiment (we would like to do that 
//                1) if the experiment is normal matching, to control for the randomness in the assignment cost 
//                2)if the experiment is the permutation test, to randomly shuffle the graphs and repeat the experiment
//@controlsLabel: text that identifies the group of control subjects 
//@patientsLabel: text that identifies the group of patient subjects (note that, these two groups might be genders, different time points, etc.)
//===products of the function===
//@similarityScores: matching score between the two graphs for each subject
//@matchesMatrix: for each patient, to which node j in the control does the node i in patient get mapped to
void GroupwiseExperiment::groupwiseMatchingExperiment(int numPermutation,std::string controlsLabel,std::string patientsLabel)
{
    dataset->preprocessGraphs();
    Matcher matcher(dataset);

    vector<int> patientsOrder,controlsOrder;
    vector<string> patientIDs,controlIDs;

    identifyGroups(samplesFile,controlsLabel,patientsLabel,controlsOrder,patientsOrder,controlIDs,patientIDs);

    int numOfPatients = patientsOrder.size();
    int numOfControls = controlsOrder.size();
    int numOfNodes= dataset->getSizeOfAGraph();

    //allocate memory to record the results of objective values and the corresponding mappings for per node per subject
    //for a matchDatasetSubjectwise call to the dataset with a given set of parameters
    float** similarityScores=Utility::allocate2Dmemory<float>(numOfPatients,numOfControls);
    //to be used for calculating pairwise matching accuracy among subjects (which we will eventually use for calculating average matching accuracy per subject
    float** matchingAccuracy=Utility::allocate2Dmemory<float>(numOfPatients,numOfControls);

    //allocate memory to keep the account of which nodes matches which node for different beta values
    //we initialize the contents of the matchesMatrix to -1 just to make sure ever node is mapped 
    //to a non-zero value at the end of mapping
    int ***matchesMatrix = Utility::allocate3Dmemory<int>(numOfPatients,numOfControls,numOfNodes);
    Utility::fillMatrix<int>(matchesMatrix,-1,numOfPatients,numOfControls,numOfNodes);

    //we repeat the experiment, since the diagonals are being randomly set in the matrix initializer class
    for(int perm=0;perm<numPermutation;perm++)
    {
        //run the matcher for the entire dataset for this value of the beta, and store the matching scores in the resultMatrix
        matcher.matchGroups(similarityScores,matchesMatrix,patientsOrder,controlsOrder,edgeType1,edgeType2);

        //calculate average matching score of each patient relative to the control set at overall graph level
        vector<float> avgMatchingAccuracy(numOfPatients);//initialize an array of all zeros as the average score
        Utility::fillMatrix<float>(matchingAccuracy,0.0,numOfPatients,numOfControls);
        for(int i=0;i<numOfPatients;i++)
        {
            int count=0;
            for(int j=0;j<numOfControls;j++)
            {
                if(patientIDs[i].compare(controlIDs[j])==0)
                    continue;
                for(int k=0;k<numOfNodes;k++)
                {
                    if(matchesMatrix[i][j][k]==k)
                        matchingAccuracy[i][j]+=1;
                }
                matchingAccuracy[i][j]=(matchingAccuracy[i][j]*100.0)/float(numOfNodes);
                avgMatchingAccuracy[i] += matchingAccuracy[i][j];
                count++;
            }
            if(count!=0)
                avgMatchingAccuracy[i] /= float(count);
        }

        saveMatchingResults(numOfPatients,numOfControls,numOfNodes,matchesMatrix,NULL,similarityScores,&avgMatchingAccuracy,matchingAccuracy,outputPath+"_"+Utility::stringify(perm)+".res");
    }

    //free the memory that we have allocated for the resultMatrix
    Utility::free2Dmemory<float>(similarityScores,numOfPatients);
    Utility::free2Dmemory<float>(matchingAccuracy,numOfPatients);
    Utility::free3Dmemory<int>(matchesMatrix,numOfPatients,numOfControls);
}

//given a set of patient and controls subjects and a modality (structure/function)
//this function calculates the correlation of each subject in one group relative to each subject in the other group
//at the end, the correlation matrix and also the average correlation of every patient relative to every control are saved into file
//@correlationMatrix: correlation of patients to controls (or everybody to everybody if both groups are labelled as "all") are stored in this matrix
void GroupwiseExperiment::groupwiseCorrelationExperiment(std::string controlsLabel,std::string patientsLabel)
{
    dataset->preprocessGraphs();
    MatrixSimilarity matrixSimilarity(dataset);

    vector<int> patientsOrder,controlsOrder;
    vector<string> patientIDs,controlIDs;
    identifyGroups(samplesFile,controlsLabel,patientsLabel,controlsOrder,patientsOrder,controlIDs,patientIDs);

    int numOfPatients = patientsOrder.size();
    int numOfControls = controlsOrder.size();
    int numOfNodes= dataset->getSizeOfAGraph();

    //allocate memory to record the results of objective values and the corresponding mappings for per node per subject
    //for a matchDatasetSubjectwise call to the dataset with a given set of parameters
    float** correlationMatrix=Utility::allocate2Dmemory<float>(numOfPatients,numOfControls);

    //calculate pairwise correlation between patients and controls, and store the correlation scores in the correlationMatrix
    matrixSimilarity.correlateGroups(correlationMatrix,patientsOrder,controlsOrder,edgeType1,edgeType2);

    //calculate average score of each patient relative to the control set
    vector<float> avgCorrelation(numOfPatients);//initialize an array of all zeros as the average score
    for(int i=0;i<numOfPatients;i++)
    {
        int count=0;
        for(int j=0;j<numOfControls;j++)
        {
            //if the patient and the control are the same subject, then do not take this self correlation into account in calculation of the average
            if(patientIDs[i].compare(controlIDs[j])!=0)
            {
                avgCorrelation[i] += correlationMatrix[i][j];
                count++;
            }
        }
        if(count!=0)
            avgCorrelation[i] /= float(count);
    }

    saveSimilarityResults(numOfPatients,numOfControls,numOfNodes,&avgCorrelation,correlationMatrix,outputPath+".res",true);

    //free the memory that we have allocated for the resultMatrix
    Utility::free2Dmemory<float>(correlationMatrix,numOfPatients);
}

//given a set of patient and controls subjects and a modality (structure/function)
//this function calculates the distance between the connectomes of the same modality among patients compared to control subjects
//using subtraction, l1, or l2 distance between matrices and returns a single number for each subject pair.
//===inputs===
//@modality: <structure | function>
//@distanceMeasure: <subtract | l1 | l2> the method to use for calculating the distance between two matrix
//@controlsLabel: text that identifies the group of control subjects 
//@patientsLabel: text that identifies the group of patient subjects (note that, these two groups might be genders, different time points, etc.)
//===products of the function===
//@similarityScores: matching score between the two graphs for each subject
void GroupwiseExperiment::groupwiseMatrixDistanceExperiment(std::string distanceMeasure,std::string controlsLabel,std::string patientsLabel)
{
    dataset->preprocessGraphs();
    MatrixSimilarity matrixSimilarity(dataset);

    vector<int> patientsOrder,controlsOrder;
    vector<string> patientIDs,controlIDs;
    identifyGroups(samplesFile,controlsLabel,patientsLabel,controlsOrder,patientsOrder,controlIDs,patientIDs);

    int numOfPatients = patientsOrder.size();
    int numOfControls = controlsOrder.size();
    int numOfNodes= dataset->getSizeOfAGraph();

    //allocate memory to record the results of objective values and the corresponding mappings for per node per subject
    //for a matchDatasetSubjectwise call to the dataset with a given set of parameters
    float** similarityScores=Utility::allocate2Dmemory<float>(numOfPatients,numOfControls);

    //run the matcher for the entire dataset and store the similarity scores in the similarityScores matrix
    matrixSimilarity.calculateMatrixDistanceBetweenGroups(similarityScores,patientsOrder,controlsOrder,distanceMeasure,edgeType1,edgeType2);

    //calculate average score of each patient relative to the control set
    vector<float> avgSimilarity(numOfPatients);//initialize an array of all zeros as the average score
    for(int i=0;i<numOfPatients;i++)
    {
        int count=0;
        for(int j=0;j<numOfControls;j++)
        {
            //if the patient and the control are the same subject, then do not take this self correlation into account in calculation of the average
            if(patientIDs[i].compare(controlIDs[j])!=0)
            {
                avgSimilarity[i] += similarityScores[i][j];
                count++;
            }
        }
        if(count!=0)
            avgSimilarity[i] /= float(count);
    }

    saveSimilarityResults(numOfPatients,numOfControls,numOfNodes,&avgSimilarity,similarityScores,outputPath+".res",true);

    //free the memory that we have allocated for the resultMatrix
    Utility::free2Dmemory<float>(similarityScores,numOfPatients);
}

//given a set of patient and controls subjects and a modality (structure/function)
//this function subtracts the connectome of each subject from the rest of the group and aerages the edgewise differences
//at the end, this function produces a distance connectome relative to a population for each subject. Thus, it does not produce a scalar but a matrix for each subject.
//@dissimilarityMatrices: for each patient, average of subtraction of each patient from all controls
void GroupwiseExperiment::groupwiseSubjectSpecificMatrixSubtractionExperiment(std::string controlsLabel,std::string patientsLabel)
{
    dataset->preprocessGraphs();
    MatrixSimilarity matrixSimilarity(dataset);

    vector<int> patientsOrder,controlsOrder;
    vector<string> patientIDs,controlIDs;
    identifyGroups(samplesFile,controlsLabel,patientsLabel,controlsOrder,patientsOrder,controlIDs,patientIDs);

    int numOfPatients = patientsOrder.size();
    int numOfNodes= dataset->getSizeOfAGraph();

    //allocate memory to record the results of objective values and the corresponding mappings for per node per subject
    //for a matchDatasetSubjectwise call to the dataset with a given set of parameters
    float*** dissimilarityMatrices=Utility::allocate3Dmemory<float>(numOfPatients,numOfNodes,numOfNodes);

    //run the matcher for the entire dataset for this value of the beta, and store the matching scores in the resultMatrix
    matrixSimilarity.subtractEveryoneFromMembersOfGroup(dissimilarityMatrices,patientsOrder,controlsOrder,patientIDs,controlIDs,edgeType1,edgeType2);

    for(int i=0;i<numOfPatients;i++)
        Utility::saveMatrixToFile(outputPath+patientIDs[i]+".diff",dissimilarityMatrices[i],numOfNodes,numOfNodes,'\t','\n',4);

    //free the memory that we have allocated for the resultMatrix
    Utility::free3Dmemory<float>(dissimilarityMatrices,numOfPatients,numOfNodes);
}

//given a set of patient and controls subjects and a modality (structure/function)
//for each edge of each patient, this function calculates the zScore of the edge weight wrt the connectomes of the set of controls
//and produces the following data
//@zScoreMatrices: for each patient, zScore of each edge relative to control population
void GroupwiseExperiment::groupwiseSubjectSpecificMatrixZScoreExperiment(std::string controlsLabel,std::string patientsLabel)
{
    dataset->preprocessGraphs();
    MatrixSimilarity matrixSimilarity(dataset);

    vector<int> patientsOrder,controlsOrder;
    vector<string> patientIDs,controlIDs;
    identifyGroups(samplesFile,controlsLabel,patientsLabel,controlsOrder,patientsOrder,controlIDs,patientIDs);

    int numOfPatients = patientsOrder.size();
    int numOfNodes= dataset->getSizeOfAGraph();

    //allocate memory to record the results of objective values and the corresponding mappings for per node per subject
    //for a matchDatasetSubjectwise call to the dataset with a given set of parameters
    float*** zScoreMatrices=Utility::allocate3Dmemory<float>(numOfPatients,numOfNodes,numOfNodes);

    //run the matcher for the entire dataset for this value of the beta, and store the matching scores in the resultMatrix
    matrixSimilarity.calculateZscoreOfEdgesForEveryoneWRTGroup(zScoreMatrices,patientsOrder,controlsOrder,patientIDs,controlIDs,edgeType1,edgeType2);

    for(int i=0;i<numOfPatients;i++)
        Utility::saveMatrixToFile(outputPath+patientIDs[i]+".zSc",zScoreMatrices[i],numOfNodes,numOfNodes,'\t','\n',4);

    //free the memory that we have allocated for the resultMatrix
    Utility::free3Dmemory<float>(zScoreMatrices,numOfPatients,numOfNodes);
}

//This function saves the beta/numodes/numSubjects and the matches matrix of an experiment into either file or to stdout
//@matchesMatrix: each row of the matrix contains the mapping info for a subject which consists of NUM_OF_NODES (number of regions) numbers
//                where ith column having value j means ith structural region is mapped to jth functional region for this subject
void GroupwiseExperiment::saveMatchingResults(int numOfGroup1, int numOfGroup2, int numOfNodes, int ***matchesMatrix, std::vector<float> *averageSimilarityScores, float **similarityScores, std::vector<float> *averageMatchingAccuracies, float **matchingAccuracies, string outfile, bool printInformationLines)
{
    stringstream outputStream;

    if(printInformationLines==true)
        outputStream<<"#Command line\n";
    outputStream<<commandLine;
    if(printInformationLines==true)
        outputStream<<"#numOfNodes, numOfPatients, numOfControls: "<<endl;
    outputStream<<numOfNodes<<"\t"<<numOfGroup1<<"\t"<<numOfGroup2<<endl;

    if(printAverageSimilarityScores==true && averageSimilarityScores!=NULL) ////// TODO: check HERE!!!!
    {
        if(printInformationLines==true)
            outputStream<<"#Average similarity score of each patient relative to the controls:"<<endl;
        Utility::printVector<float>(*averageSimilarityScores,outputStream,'\t','\n',4);
    }


    if(printSimilarityScores==true && similarityScores!=NULL)
    {
        if(printInformationLines==true)
            outputStream<<"#Pairwise similarity scores:"<<endl;
        Utility::printMatrix(similarityScores,numOfGroup1,numOfGroup2,outputStream);
    }

    if(printAverageMatchingAccuracies==true && averageMatchingAccuracies!=NULL) ////// TODO: check HERE!!!!
    {
        if(printInformationLines==true)
            outputStream<<"#Average matching accuracy score of each patient relative to the controls:"<<endl;
        Utility::printVector<float>(*averageMatchingAccuracies,outputStream,'\t','\n',4);
    }


    if(printMatchingAccuracies==true && matchingAccuracies!=NULL)
    {
        if(printInformationLines==true)
            outputStream<<"#Pairwise matching accuracies:"<<endl;
        Utility::printMatrix(matchingAccuracies,numOfGroup1,numOfGroup2,outputStream);
    }

    if(printMatches==true && matchesMatrix!=NULL)
    {
        if(printInformationLines==true)
            outputStream<<"#Matching pairs for each subject:"<<endl;
        for(int i=0;i<numOfGroup1;i++)
        {
            for(int j=0;j<numOfGroup2;j++)
            {
                outputStream<<i<<"\t"<<j<<"\t";
                Utility::printVector<int>(matchesMatrix[i][j],numOfNodes,outputStream);
            }
        }
    }
    outputStream<<endl;

    if(outfile.compare("")!=0)
        Utility::printOstreamToFile(outfile,outputStream);//appendOstreamToFile(outfile,outputStream);
    else
        cout<<outputStream.rdbuf();
}

//This function saves the beta/numodes/numSubjects and the matches matrix of an experiment into either file or to stdout
//@matchesMatrix: each row of the matrix contains the mapping info for a subject which consists of NUM_OF_NODES (number of regions) numbers
//                where ith column having value j means ith structural region is mapped to jth functional region for this subject
void GroupwiseExperiment::saveSimilarityResults(int numOfGroup1, int numOfGroup2, int numOfNodes, std::vector<float> *averageSimilarities, float **similarityScores, string outfile, bool printInformationLines)
{
    stringstream outputStream;

    if(printInformationLines==true)
       outputStream<<"#Command line\n";
    outputStream<<commandLine;
    if(printInformationLines==true)
       outputStream<<"#numOfNodes, numOfPatients, numOfControls: "<<endl;
    outputStream<<numOfNodes<<"\t"<<numOfGroup1<<"\t"<<numOfGroup2<<endl;

    if(printAverageSimilarityScores==true && averageSimilarities!=NULL) ////// TODO: check HERE!!!!
    {
       if(printInformationLines==true)
          outputStream<<"#Average similarity score of each patient relative to the controls:"<<endl;
       Utility::printVector<float>(*averageSimilarities,outputStream,'\t','\n',4);
    }

    if(printSimilarityScores==true && similarityScores!=NULL)
    {
       if(printInformationLines==true)
          outputStream<<"#Pairwise similarity scores:"<<endl;
       Utility::printMatrix(similarityScores,numOfGroup1,numOfGroup2,outputStream);
    }

    outputStream<<endl;

    if(outfile.compare("")!=0)
       Utility::printOstreamToFile(outfile,outputStream);//appendOstreamToFile(outfile,outputStream);
    else
       cout<<outputStream.rdbuf();
}

//given the results of a set of groupwise matching experiments (that is, same experiment is repeated several times with permutation), 
//where each file contains the similarity score of each subject relative to everyone else, and the matching nodes between pairs of subjects,
//this function calculates 1) an average similarity score matrix 2)matching matrixc file for each subject containing which node matched with which node
void GroupwiseExperiment::calculateAverageResults(std::string resultsFolder, std::string controlsLabel, std::string patientsLabel)
{
    vector<string> filenames;

    vector<int> patientsOrder,controlsOrder;
    vector<string> patientIDs,controlIDs;
    identifyGroups(samplesFile,controlsLabel,patientsLabel,controlsOrder,patientsOrder,controlIDs,patientIDs);

    Utility::loadFilePathFromFolder(resultsFolder,filenames);

    int skipline1=2;

    //First get the numNodes, numPatients, and numControls
    int numNodes, numPatients, numControls;

    std::ifstream inFile; 
    std::istringstream is;
    std::string line;
    inFile.open(filenames[0].c_str());
    if(inFile.fail()){std::cout << "File not found!!!!\n";exit(1);}
    for(int i=0;i<skipline1;i++)
       getline(inFile, line);//ignore skiplines
    getline(inFile, line);is.clear();is.str(line);
    is>>numNodes>>numPatients>>numControls;
    inFile.close();

    int skipline2=4,skipline3=skipline2+numPatients+1;

    float **similarityScores = Utility::allocate2Dmemory<float>(numPatients,numControls);
    float **similarityTemp = Utility::allocate2Dmemory<float>(numPatients,numControls);
    Utility::fillMatrix<float>(similarityScores,0.0,numPatients,numControls);

    float ***matchingMatrices = Utility::allocate3Dmemory<float>(numPatients,numNodes,numNodes);
    int **matchingTemp = Utility::allocate2Dmemory<int>(numPatients*numControls,numNodes+2);
    Utility::fillMatrix<float>(matchingMatrices,0,numPatients,numNodes,numNodes);

    for(std::vector<string>::iterator fileIter=filenames.begin();fileIter!=filenames.end();fileIter++)
    {
       Utility::loadMatrixFromFileAfterSkipLine(*fileIter,similarityTemp,numPatients,numControls,skipline2);
       Utility::elementwiseAddMatrices(similarityScores,similarityTemp,similarityScores,numPatients,numControls);

       Utility::loadMatrixFromFileAfterSkipLine(*fileIter,matchingTemp,numPatients*numControls,numNodes+2,skipline3);
    
       for(int i=0;i<numPatients;i++)
          for(int j=0;j<numControls;j++)
          {
             //if we are comparing a subject with itself, do not account for it in the matching matrix
             if(matchingTemp[i*numControls+j][0]==matchingTemp[i*numControls+j][1])
                continue;
             for(int k=0;k<numNodes;k++)
                matchingMatrices[i][k][matchingTemp[i*numControls+j][k+2]]+=1;
          }
    }
    Utility::multiplyMatrixWithScalar<float>(similarityScores,1.0/(float)filenames.size(),numPatients,numControls);
    Utility::saveMatrixToFile(outputPath+"similarityScores.txt",similarityScores,numPatients,numControls); 

    for(int i=0;i<numPatients;i++)
    {
       Utility::multiplyMatrixWithScalar<float>(matchingMatrices[i],1.0/(float)(filenames.size()*numControls),numNodes,numNodes);
       Utility::saveMatrixToFile(outputPath+patientIDs[i]+"_matching.txt",matchingMatrices[i],numNodes,numNodes,'\t','\n',4); 
    }

    Utility::free2Dmemory<float>(similarityScores,numPatients);
    Utility::free2Dmemory<float>(similarityTemp,numPatients);
    Utility::free3Dmemory<float>(matchingMatrices,numPatients,numNodes);
    Utility::free2Dmemory<int>(matchingTemp,numPatients*numControls);
}
