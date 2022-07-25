/* 
 * File:   couplingExperiment.cpp
 * Author: yusuf
 *
 * Generated on May 1, 2018, 2:02 PM
 */

#include "subjectwiseExperiment.h"
#include "dataset.h"
#include "matcher.h"
#include "matrixSimilarity.h"

using namespace std;

//<editor-fold defaultstate="collapsed" desc=" Experiments for matching structural and functional graphs: structureFunctionMappingExperiment()">

//given a set of subjects and a single parameter set
//this function runs a matching between the structural and functional graphs for each subjects
//and collects the following data
//@similarityScores: matching score between the two graphs for each subject
//@matchingAccuraciesPerNode: how many times each structural node mapped to its corresponding counterpart in the functional graph
//@matchingAccuraciesPerSubject: how many structural node mapped to its corresponding functional node for each subject
//@matchesAccMatrix: how many times each pair of structure-function node get mapped to each other
//@matchesMatrix: for each subject, to which functional node does structural node i get mapped to
void SubjectwiseExperiment::subjectwiseMatchingExperiment(int numPermutation, bool saveAvgConnectome, int seedSupplement, string shuffleType,int shuffleDataset)
{
    int numSubjects = dataset->getNumOfSubjects();
    int numNodes= dataset->getSizeOfAGraph();
   
    //allocate memory to record the results of objective values and the corresponding mappings for per node per subject
    //for a matchDatasetSubjectwise call to the dataset with a given set of parameters
    vector<float> similarityScores(numSubjects,0);

    //allocate memory to keep the account of which nodes matches which node for different beta values
    //we initialize the contents of the matchesMatrix to -1 just to make sure ever node is mapped 
    //to a non-zero value at the end of mapping
    int **matchesMatrix = Utility::allocate2Dmemory<int>(numSubjects,numNodes);

    if(numPermutation<1)
    {
        cerr<<"Error caught in subjectwiseMatchingExperiment():"<<endl;
        cerr<<"numPermutation should be >=1. However it was set to <"<<numPermutation<<">... Exiting!!!\n";
        exit(1);
    }

    for(int i=0;i<numPermutation;i++)
    {
        Dataset *dataset_, *dataset2_=NULL;
        //if we are going to repeat the experiment for whatever reason and we do shuffling at each iteration, then generate a copy of the matcher and set the matcher pointer to point to it
        //we do it like that instead of loading data from hard disk  in each iteration in order to save time (memory copy is faster than loading from disk)
        if(numPermutation>1 && shuffleType.compare("")!=0)
        {
            dataset_=new Dataset(*dataset);//I need to reload the dataset from the original globalMatcher object each and every time since shuffling changes the edge weights etc.
            if(dataset2!=NULL)
                dataset2_=new Dataset(*dataset2);//I need to reload the dataset from the original globalMatcher object each and every time since shuffling changes the edge weights etc.
        }
        else
        {
            dataset_=dataset;
            if(dataset2!=NULL)
                dataset2_=dataset2;
        }

        if(shuffleType.compare("")!=0)
        {
            if(shuffleDataset==1)
                dataset_->shuffleGraphs(i+seedSupplement,shuffleType);
            else if(shuffleDataset==2)
                dataset2_->shuffleGraphs(i+seedSupplement,shuffleType);
        }

        //Save average connectomes across subjects at each iteration if there is shuffling, and save them once at the initial iteration if no shuffling is done
        if(saveAvgConnectome==true)
        {//save the graph each time if there is shuffling, and once if no shuffling is done
            if(shuffleType.compare("")!=0 || (shuffleType.compare("")==0 && i==0))
            {
                dataset_->saveAverageConnectomeAcrossSubjects("structure",outputPath+"structureAvgConn"+Utility::stringify(i));
                dataset_->saveAverageConnectomeAcrossSubjects("function",outputPath+"functionPosAvgConn"+Utility::stringify(i));
            }
        }

        //preprocess graphs at each iteration if there is shuffling, and process it once at the initial iteration if no shuffling is done
        if(shuffleType.compare("")!=0 || (shuffleType.compare("")==0 && i==0))
        {
            dataset_->preprocessGraphs();
            if(dataset2_!=NULL)
                dataset2_->preprocessGraphs();
        }


        Matcher matcher(dataset_,dataset2_);
        Utility::fillMatrix<int>(matchesMatrix,-1,numSubjects,numNodes);
        matcher.matchEveryoneWithItself(similarityScores,matchesMatrix,edgeType1,edgeType2);

        saveSubjectwiseMatchingResults(numNodes,numSubjects,matchesMatrix,similarityScores,outputPath+"_"+Utility::stringify(i)+".res");

        //if there was more than one iteration, and the shuffling was done, this implies that we used local copies of the connectomes, which we need to free here
        if(numPermutation>1 && shuffleType.compare("")!=0)
        {
            delete dataset_;
            if(dataset2!=NULL)
                delete dataset2_;
        }
    }

    //free the memory that we have allocated for the resultMatrix
    Utility::free2Dmemory<int>(matchesMatrix,numSubjects);
}

//given a set of subjects and a single parameter set
//this function runs a matching between the structural and functional graphs for each subjects
//and collects the following data
//@similarityScores: matching score between the two graphs for each subject
//@matchingAccuraciesPerNode: how many times each structural node mapped to its corresponding counterpart in the functional graph
//@matchingAccuraciesPerSubject: how many structural node mapped to its corresponding functional node for each subject
//@matchesAccMatrix: how many times each pair of structure-function node get mapped to each other
//@matchesMatrix: for each subject, to which functional node does structural node i get mapped to
void SubjectwiseExperiment::subjectwiseDistanceExperiment(std::string distanceMeasure)
{
    dataset->preprocessGraphs();
    MatrixSimilarity matrixSimilarity(dataset);

    int numSubjects = dataset->getNumOfSubjects();
    int numNodes= dataset->getSizeOfAGraph();

    //allocate memory to record the results of objective values and the corresponding mappings for per node per subject
    //for a matchDatasetSubjectwise call to the dataset with a given set of parameters
    vector<float> similarityScores;//(numSubjects,0);

    matrixSimilarity.calculateMatrixDistanceOfEveryoneToItself(edgeType1,edgeType2,similarityScores,distanceMeasure);

    saveSubjectwiseMatchingResults(numNodes,numSubjects,NULL,similarityScores,outputPath+".res");
}

//This function saves the beta/numodes/numSubjects and the matches matrix of an experiment into either file or to stdout
//@matchesMatrix: each row of the matrix contains the mapping info for a subject which consists of NUM_OF_NODES (number of regions) numbers
//                where ith column having value j means ith structural region is mapped to jth functional region for this subject
void SubjectwiseExperiment::saveSubjectwiseMatchingResults(int numNodes, int numSubjects, int **matchesMatrix, vector<float>& similarityScores, string outfile)
{
    stringstream outputStream;

    outputStream<<"#Command line\n";
    outputStream<<commandLine;
    outputStream<<"#numNodes, numSubjects: "<<endl;
    outputStream<<numNodes<<"\t"<<numSubjects<<endl;

    if(printSimilarityScores==true && similarityScores.size()!=0)
    {
        outputStream<<"#Similarity scores:"<<endl;
        Utility::printVector<float>(similarityScores,outputStream,'\t',' ',6);//for DEBUG
        outputStream<<endl;
    }

    if(printMatches==true && matchesMatrix!=NULL)
    {
        outputStream<<"#Matching pairs for each subject:"<<endl;
        Utility::printMatrix<int>(matchesMatrix,numSubjects,numNodes,outputStream);
    }
    outputStream<<endl;

    if(outfile.compare("")!=0)
        Utility::printOstreamToFile(outfile,outputStream);//appendOstreamToFile(outfile,outputStream);
    else
        cout<<outputStream.rdbuf();
}
//</editor-fold>

//<editor-fold defaultstate="collapsed" desc=" str-func correlation: saveStructureFunctionCorrelationofSingleSubject()">
//this function saves the upper triangular of the structural and functional connectomes into a single file in two lines in vector form
void SubjectwiseExperiment::saveStructureFunctionCorrelationOfSingleSubject(int subjectId)
{
    int numSubjects = dataset->getNumOfSubjects();
    int numNodes= dataset->getSizeOfAGraph();//number of nodes in a graph

    float **structuralConnectome = Utility::allocate2Dmemory<float>(numNodes,numNodes);
    float **functionalConnectome = Utility::allocate2Dmemory<float>(numNodes,numNodes);

    dataset->preprocessGraphs();

    dataset->getConnectivityMatrixForGraph(subjectId,Edge::STRUCTURAL_CONNECTIVITY,structuralConnectome,"notImportant");

    dataset->getConnectivityMatrixForGraph(subjectId,Edge::FUNCTIONAL_CONNECTIVITY,functionalConnectome,"full");

    vector<float> vec1,vec2;
    vec1 = Utility::copyOffDiagonalUpperTriangularOfMatrixIntoVector(structuralConnectome,numNodes);
    vec2 = Utility::copyOffDiagonalUpperTriangularOfMatrixIntoVector(functionalConnectome,numNodes);

    outputPath += "_" + initializerParameters.functionalConnectivity + "_" + initializerParameters.pathType;

    Utility::saveVectorToFile(outputPath+"_vec.txt",vec1,'\t','\n',6);
    Utility::appendVectorToFile(outputPath+"_vec.txt",vec2,'\t','\n',6);

    Utility::free2Dmemory(structuralConnectome,numNodes);
    Utility::free2Dmemory(functionalConnectome,numNodes);
}

//this function saves the average of structural and functional connectivity matrices for the entire dataset into two separate files
//also saves the upper triangular of the average structural and functional connectomes into a single file in two lines in vector form
//For the structural connectivity, it saves path as direct,shortest,strongest, or communicability.
void SubjectwiseExperiment::saveAverageStructureFunctionCorrelationOfAllSubjects()
{
    int numSubjects = dataset->getNumOfSubjects();
    int numNodes= dataset->getSizeOfAGraph();//number of nodes in a graph

    float **connectome = Utility::allocate2Dmemory<float>(numNodes,numNodes);
    float **structuralConnectome = Utility::allocate2Dmemory<float>(numNodes,numNodes);
    float **functionalConnectome = Utility::allocate2Dmemory<float>(numNodes,numNodes);

    Utility::fillMatrix<float>(structuralConnectome,0.0,numNodes,numNodes);
    Utility::fillMatrix<float>(functionalConnectome,0.0,numNodes,numNodes);

    dataset->preprocessGraphs();

    for(int i=0;i<numSubjects;i++)
    {
        dataset->getConnectivityMatrixForGraph(i,Edge::STRUCTURAL_CONNECTIVITY,connectome,"notImportant");
        Utility::elementwiseAddMatrices(connectome,structuralConnectome,structuralConnectome,numNodes,numNodes);

        dataset->getConnectivityMatrixForGraph(i,Edge::FUNCTIONAL_CONNECTIVITY,connectome,"full");
        Utility::elementwiseAddMatrices(connectome,functionalConnectome,functionalConnectome,numNodes,numNodes);
    }
    Utility::multiplyMatrixWithScalar<float>(structuralConnectome,1.0/(float)numSubjects,numNodes,numNodes);
    Utility::multiplyMatrixWithScalar<float>(functionalConnectome,1.0/(float)numSubjects,numNodes,numNodes);

    vector<float> vec1,vec2;
    vec1 = Utility::copyOffDiagonalUpperTriangularOfMatrixIntoVector(structuralConnectome,numNodes);
    vec2 = Utility::copyOffDiagonalUpperTriangularOfMatrixIntoVector(functionalConnectome,numNodes);

    outputPath += "_" + initializerParameters.functionalConnectivity + "_" + initializerParameters.pathType;

    Utility::saveVectorToFile(outputPath+"_vec.txt",vec1,'\t','\n',6);
    Utility::appendVectorToFile(outputPath+"_vec.txt",vec2,'\t','\n',6);

    Utility::free2Dmemory(connectome,numNodes);
    Utility::free2Dmemory(structuralConnectome,numNodes);
    Utility::free2Dmemory(functionalConnectome,numNodes);
}

//</editor-fold>
