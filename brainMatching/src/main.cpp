/* 
 * File:   main.cpp
 * Author: yusuf
 *
 * Generated on September 14, 2016, 1:46 PM
 */

//use http://tclap.sourceforge.net/ for accepting command line parameters at some point

#include <iostream> //cin,cout
#include <cstdlib> //atoi,stof, etc.
#include <string> //std::stof Note: stof requires C++11 as the compiler standard
#include "time.h"
#include "graph.h"
#include "dataset.h"
#include "matcher.h"
#include "experiment.h"
#include "subjectwiseExperiment.h"
#include "groupwiseExperiment.h"
#include "utility.h"
#include "tools.h"
#include "test.h"

using namespace std;

void commandLineHelp(string helpSubcategory);

int main(int argc, char** argv) 
{
   srand(time(NULL));
   
   //<editor-fold defaultstate="collapsed" desc=" Read parameters">
   InitializerParameters initializerParameters;
   //first load the parameters if they exist, if not, they will be assigned with ""
   initializerParameters.pathType = string(Utility::getCmdOption(argv,argv+argc,"-pathType"));//-pathType: direct,uShortest,wShortest,wCommunicability,uCommunicability,searchInformation,pathTransitivity
   initializerParameters.functionalConnectivity=string(Utility::getCmdOption(argv,argv+argc,"-funcConn"));//positive/negative
   initializerParameters.assignmentCostMode=string(Utility::getCmdOption(argv,argv+argc,"-assCost"));//features,edges<IgnoreDiag|IncludeDiag<ZeroDiag|RandDiag> >,
   initializerParameters.adjustGraphParameters=string(Utility::getCmdOption(argv,argv+argc,"-preprocessGraphs"));//traffic,logScale<All,Edges,Nodes><Str,Fun>,normalize<All,Edges,Nodes>
   string modality = string(Utility::getCmdOption(argv,argv+argc,"-modality")); // <str_str | str_func | func_func>
   string samplesFile = string(Utility::getCmdOption(argv,argv+argc,"-samples"));
   string outputPath = string(Utility::getCmdOption(argv,argv+argc,"-outputPath"));
   string resultsPath = string(Utility::getCmdOption(argv,argv+argc,"-results"));
   string truthFile = string(Utility::getCmdOption(argv,argv+argc,"-truth"));
   string normalizerFile = string(Utility::getCmdOption(argv,argv+argc,"-normalizer"));
   bool printMatches = Utility::cmdOptionExists( argv, argv+argc, "-printMatches" );
   bool printSimilarityScores = Utility::cmdOptionExists( argv, argv+argc, "-printSimilarity" );
   bool printAverageSimilarityScores = Utility::cmdOptionExists( argv, argv+argc, "-printAverageSimilarity" );
   bool printMatchingAccuracies = Utility::cmdOptionExists( argv, argv+argc, "-printMatchingAccuracies" );
   bool printAverageMatchingAccuracies = Utility::cmdOptionExists( argv, argv+argc, "-printAverageMatchingAccuracies" );
   bool verbose = Utility::cmdOptionExists( argv, argv+argc, "-verbose" );
      
   if(modality.compare("func_func")==0)
      initializerParameters.functionalConnectivity2=string(Utility::getCmdOptionWithOrder(argv, argv+argc,"-funcConn",2));//same type of functional connectomes will be compared groupwise
   
   string commandLine = Utility::getCommandLine(argc,argv);//to be used for keeping record of which command line is run to produce certain results
   
   if(Utility::cmdOptionExists( argv, argv+argc, "-pathLength" ))
      initializerParameters.maxPathLength=std::atoi(Utility::getCmdOption(argv, argv+argc,"-pathLength"));
   else
      initializerParameters.maxPathLength=-1;//this is used for the communicability calculations. If not set, default value is -1, indicating that it should not be used
      
   //</editor-fold>
   
   //<editor-fold defaultstate="collapsed" desc=" Load dataset(s)">
   Dataset dataset[2];
   int numDataset;
   if( Utility::cmdOptionExists( argv, argv+argc, "-data" ) )
   {
      // experiments can load one or two datasets, and these datasets could be loaded from previously generated .grp files or from connectome matrix files
      // after <-data> flag, enter the type of the datasets to be loaded in <graph | matrix> and then enter the number of datasets to be loaded
      // use cases: -data graph 2 <graphsPath1> <graphsPath2>
      //            -data graph 1 <graphsPath>
      //            -data matrix 2 -dti <dtiPath1> <dtiPath2> -fmri <fmriPath1> <fmriPath2> -feature <featuresPath1> <featuresPath2>
      //            -data matrix 1 -dti <dtiPath> -fmri <fmriPath> -feature <featuresPath>
      string datasetType = string(Utility::getCmdOptionWithOrder(argv, argv+argc,"-data",1)); 
      numDataset= std::atoi(Utility::getCmdOptionWithOrder(argv, argv+argc,"-data",2));   

      for(int i=0;i<numDataset;i++)
      {
         //initialize parameters
         dataset[i].setInitializerParameters(initializerParameters); 

         //load graphs or matrices
         if(datasetType.compare("graph")==0)
            dataset[i].loadDataset(string(Utility::getCmdOptionWithOrder(argv, argv+argc,"-data",i+3)),samplesFile,verbose);
         else if(datasetType.compare("matrix")==0)
            dataset[i].loadDataset(string(Utility::getCmdOptionWithOrder(argv, argv+argc,"-dti",i+1)),string(Utility::getCmdOptionWithOrder(argv, argv+argc,"-fmri",i+1)),string(Utility::getCmdOptionWithOrder(argv, argv+argc,"-feature",i+1)),samplesFile,verbose);

         //preprocess graphs inside experiment functions
      }

      if(numDataset==2 && dataset[0].getNumOfSubjects()!=dataset[1].getNumOfSubjects())
      {
         cerr<<"Two datasets are loaded, albeit containing different number of subjects. Check your data!!! Exiting !...\n";
         exit(1);
      }
   }

   //</editor-fold>
   
   //now, go through the decision tree of program parameters
   if( Utility::cmdOptionExists( argv, argv+argc, "-help" ) )
      commandLineHelp(string(Utility::getCmdOption(argv,argv+argc,"-help"))); //print help for a certain help subcategory
   else if( Utility::cmdOptionExists( argv, argv+argc, "-experiment" ) )
   //<editor-fold defaultstate="collapsed" desc=" Experiments: -experiment <groupwise, subjectwise, subjectwisePerm, sampleGraphs, listOfSamples>">
	{  
      string experimentType = string(Utility::getCmdOptionWithOrder(argv, argv+argc,"-experiment",1));//usage: -experiment groupwise match
      string experimentSubtype = string(Utility::getCmdOptionWithOrder(argv, argv+argc,"-experiment",2));//usage: -experiment subjectwise match
      bool saveAvgConnectome = Utility::cmdOptionExists(argv,argv+argc,"-saveAvgConnectome");
      
      int numPermutation = 1;
      if(Utility::cmdOptionExists(argv,argv+argc,"-permutation"))
         numPermutation = std::atoi(Utility::getCmdOption(argv,argv+argc,"-permutation"));
      
      if( experimentType.compare("groupwise" )==0 )
      {
         string controlsLabel = string(Utility::getCmdOptionWithOrder(argv, argv+argc,"-groups",1)); //usage: -groups c p
         string patientsLabel = string(Utility::getCmdOptionWithOrder(argv, argv+argc,"-groups",2)); //usage: -groups c all
         //if we are comparing function across subjects, which function of the first is being matched with which function of the second
        
         GroupwiseExperiment exp(initializerParameters,commandLine);
         exp.setFiles(samplesFile,outputPath,normalizerFile);
         exp.setDataset(&dataset[0]); //groupwise experiments uses a single dataset where subjects are being compared with respect to each other
         exp.setModality(modality); //Sets edgeType1 and edgeType2 //modality should look like: str_func, str_str, func_func or structure_structure, structure_function, function_function
         exp.setPrint(printMatches,printSimilarityScores,printAverageSimilarityScores,printMatchingAccuracies,printAverageMatchingAccuracies);//shall I print which object nodes match to which label nodes
         
         if( experimentSubtype.compare("match" )==0 )
         {
            //this function runs a matching between the connectomes of the same modality among one group (such as all or patients) 
            //compared to another group of subjects (such as healthy controls or all subjects)
            exp.groupwiseMatchingExperiment(numPermutation,controlsLabel,patientsLabel);
            //sample: ./brainMatch -experiment groupwise match -groups c p -data graph 1 ~/data/TBI/graphs/Desikan86/thresh15/ -samples ~/data/TBI/resources/commonList_allTimePoints_1_Desikan86.txt -printMatches -printSimilarity -modality str_str -outputPath results -pathType direct -assCost features_edgesIgnoreDiag -preprocessGraphs logScaleEdgesStr_normalizeEdges
            //        ./brainMatch -experiment groupwise match -groups c p -data matrix 1 ~/data/TBI/connectomes/Desikan86/thresh15/ -samples ~/data/TBI/resources/commonList_allTimePoints_1_Desikan86.txt -printMatches -printSimilarity -modality str_str -outputPath results -pathType direct -assCost features_edgesIgnoreDiag -preprocessGraphs logScaleEdgesStr_normalizeEdges
            //        ./brainMatch -experiment groupwise match -groups c all -permutation 3 -data graph 1 ~/data/TBI/graphs/Desikan86/thresh15/ -samples ~/data/TBI/resources/commonList_allTimePoints_1_Desikan86.txt -printMatches -printSimilarity -modality str_str -outputPath results -pathType direct
         }
         else if( experimentSubtype.compare("correlation" )==0 )
         {
            //this function runs a matching between the connectomes of the same modality among one group (such as all or patients) 
            //compared to another group of subjects (such as healthy controls or all subjects)
            exp.groupwiseCorrelationExperiment(controlsLabel,patientsLabel);
            //sample: ./brainMatch -experiment groupwise correlation -groups c p -data graph 1 ~/data/TBI/graphs/Desikan86/thresh15/ -samples ~/data/TBI/resources/commonList_allTimePoints_1_Desikan86.txt -printMatches -printSimilarity -modality str_str -outputPath results -pathType direct -assCost features_edgesIgnoreDiag -preprocessGraphs logScaleEdgesStr_normalizeEdges
            //        ./brainMatch -experiment groupwise correlation -groups c p -data matrix 1 ~/data/TBI/connectomes/Desikan86/thresh15/ -samples ~/data/TBI/resources/commonList_allTimePoints_1_Desikan86.txt -printMatches -printSimilarity -modality str_str -outputPath results -pathType direct -assCost features_edgesIgnoreDiag -preprocessGraphs logScaleEdgesStr_normalizeEdges
            //        ./brainMatch -experiment groupwise correlation -groups c all -permutation 3 -data graph 1 ~/data/TBI/graphs/Desikan86/thresh15/ -samples ~/data/TBI/resources/commonList_allTimePoints_1_Desikan86.txt -printMatches -printSimilarity -modality str_str -outputPath results -pathType direct
         }
         else if( experimentSubtype.compare("distance" )==0 )
         {
            //this function calculates the distance between the connectomes of the same modality among patients compared to control subjects
            //using subtraction, l1, or l2 distance between matrices and returns a single number for each subject pair.
            exp.groupwiseMatrixDistanceExperiment(string(Utility::getCmdOption(argv,argv+argc,"-distanceMeasure")),controlsLabel,patientsLabel);
            //sample: ./brainMatch -experiment groupwise distance -groups c all -data graph 1 ~/data/TBI/graphs/Desikan86/thresh15/ -samples ~/data/TBI/resources/commonList_allTimePoints_1_Desikan86.txt -printMatches -printSimilarity -modality structure -outputPath results -pathType direct
         }
         else if( experimentSubtype.compare("subSpecZScore" )==0 )
         {
            //For each edge of each patient, this function calculates the zScore of the edge weight wrt the connectomes of the set of controls
            exp.groupwiseSubjectSpecificMatrixZScoreExperiment(controlsLabel,patientsLabel);
            //sample: ./brainMatch -experiment groupwise subSpecZScore -groups c all -data graph 1 ~/data/TBI/graphs/Desikan86/thresh15/ -samples ~/data/TBI/resources/commonList_allTimePoints_1_Desikan86.txt -printMatches -printSimilarity -modality structure -outputPath results -pathType direct
         }
         else if( experimentSubtype.compare("subSpecDifference" )==0 )
         {
            //this function subtracts the connectome of each subject from the rest of the group and aerages the edgewise differences
            //at the end, this function produces a distance connectome relative to a population for each subject. Thus, it does not produce a scalar but a matrix for each subject.
            exp.groupwiseSubjectSpecificMatrixSubtractionExperiment(controlsLabel,patientsLabel);
            //sample: ./brainMatch -experiment groupwise subSpecDifference -groups c all -data graph 1 ~/data/TBI/graphs/Desikan86/thresh15/ -samples ~/data/TBI/resources/commonList_allTimePoints_1_Desikan86.txt -printMatches -printSimilarity -modality structure -outputPath results -pathType direct
         }
         else if( experimentSubtype.compare("average" )==0 )
         {
            //given the results of a set of groupwise matching experiments (that is, same experiment is repeated several times with permutation), 
            //where each file contains the similarity score of each subject relative to everyone else, and the matching nodes between pairs of subjects,
            //this function calculates 1) an average similarity score matrix 2)matching matrix file for each subject containing which node matched with which node
            exp.calculateAverageResults(string(Utility::getCmdOption(argv,argv+argc,"-resultsFolder")),controlsLabel,patientsLabel);
            //sample: ./brainMatch -experiment groupwise average -groups Control all -resultsFolder functionFull/ -samples ~/data/SmithDOS/resources/dosAll.txt -outputPath ./
         }
      }
      else if( experimentType.compare("subjectwise" )==0 )
      {
         SubjectwiseExperiment exp(initializerParameters,commandLine);
         exp.setFiles(samplesFile,outputPath,normalizerFile);
         if(numDataset==1)
            exp.setDataset(&dataset[0]);
         else if(numDataset==2)
            exp.setDataset(&dataset[0],&dataset[1]);
         
         exp.setModality(modality); //Sets edgeType1 and edgeType2
         exp.setPrint(printMatches,printSimilarityScores,printAverageSimilarityScores,printMatchingAccuracies,printAverageMatchingAccuracies);//shall I print which object nodes match to which label nodes
         
         if(experimentSubtype.compare("match" )==0 )
         {
            if(Utility::cmdOptionExists(argv,argv+argc,"-shuffle")==false)//that is, if we are not making permutation test, but actual matching
            {
               exp.subjectwiseMatchingExperiment(numPermutation,saveAvgConnectome);
               //sample: =====str_func=====
               //        ./brainMatch -experiment subjectwise match -modality str_func -funcConn positive -data graph 1 ~/data/PNC/graphs/structureFunctionMappingProject/notNormalized/ -samples ~/data/PNC/src/samples/samples_4.txt -printMatches -printSimilarity -outputPath LinAss -pathType direct -assCost edgesIgnoreDiag -preprocessGraphs logScaleEdgesStructure_normalizeEdges
               //        ./brainMatch -experiment subjectwise match -modality str_func -funcConn positive -data graph 1 ~/data/PNC/graphs/structureFunctionMappingProject/notNormalized/ -samples ~/data/PNC/src/samples/allSamples.txt -printMatches -printSimilarity -outputPath QAPD -pathType wShortest
               //        ./brainMatch -experiment subjectwise match -modality str_func -funcConn positive -data graph 1 ~/data/PNC/graphs/structureFunctionMappingProject/notNormalized/ -samples ~/data/PNC/src/samples/allSamples.txt -printMatches -printSimilarity -outputPath QAPD -pathType wCommunicability -pathLength 3
               //        ./brainMatch -experiment subjectwise match -modality str_func -funcConn positive -data graph 1 ~/data/PNC/graphs/structureFunctionMappingProject/notNormalized/ -samples ~/data/PNC/src/samples/allSamples.txt -printMatches -printSimilarity -outputPath QAPD -pathType wShortest -saveAvgConnectome
               //        ./brainMatch -experiment subjectwise match -modality str_func -funcConn positive -permutation 333 -data graph 1 ~/data/PNC/graphs/structureFunctionMappingProject/notNormalized/ -samples ~/data/PNC/src/samples/allSamples.txt -printMatches -printSimilarity -outputPath QAPD -pathType wShortest
               //sample: =====func_func=====
               //        ./brainMatch -experiment subjectwise match -modality func_func  -funcConn positive negative -data graph 1 ~/data/PNC/graphs/structureFunctionMappingProject/notNormalized/ -samples ~/data/PNC/src/samples/samples_4.txt -printMatches -printSimilarity -outputPath LinAss -solver LinAss
               
            }
            else//that is, we are calculating permutation testing, which involves shuffling the graphs
            {
               string shuffleType=string(Utility::getCmdOptionWithOrder(argv, argv+argc,"-shuffle",1)); //<structure | function>
               int shuffleDataset=std::atoi(Utility::getCmdOptionWithOrder(argv, argv+argc,"-shuffle",2)); //<1 | 2> , which dataset shall we shuffle?
               if(shuffleDataset!=1 && shuffleDataset!=2)
               {
                  cerr<<"Incorrect number entered for the shuffle dataset ID! Exiting...\n";
                  exit(1);
               }
               int seed=rand();
               if(Utility::cmdOptionExists(argv,argv+argc,"-seed"))
                  seed = std::atoi(Utility::getCmdOption(argv,argv+argc,"-seed")); //<some int val>
               
               exp.subjectwiseMatchingExperiment(numPermutation,saveAvgConnectome,seed,shuffleType,shuffleDataset);
               //sample: =====str_func=====
               //        ./brainMatch -experiment subjectwise match -modality str_func -permutation 333 -seed 10 -shuffle structure 1 -data graph 1 ~/data/PNC/graphs/structureFunctionMappingProject/notNormalized/ -samples ~/data/PNC/src/samples/allSamples.txt -printMatches -printSimilarity -funcConn positive -outputPath results/QAPD -pathType wShortest -saveAvgConnectome
               //        ./brainMatch -experiment subjectwise match -modality str_func -permutation 333 -seed 10 -shuffle function 1 -data graph 1 ~/data/PNC/graphs/structureFunctionMappingProject/notNormalized/ -samples ~/data/PNC/src/samples/allSamples.txt -printMatches -printSimilarity -funcConn positive -outputPath results/QAPD -pathType wShortest -saveAvgConnectome
               //sample: =====func_func=====
               //        ./brainMatch -experiment subjectwise match -modality func_func -permutation 333 -seed 10 -shuffle function 1 -data graph 1 ~/data/PNC/graphs/structureFunctionMappingProject/notNormalized/ -samples ~/data/PNC/src/samples/allSamples.txt -printMatches -printSimilarity -funcConn positive negative -outputPath results/QAPD -solver LinAss
            }
         }
         else if(experimentSubtype.compare("correlation")==0)
         {
            //pathType == uCommunicability, wCommunicability, uShortest, wShortest, direct
            //corrType == positive, negative
            //strFuncCorrelation == average, single
            string allSingle = string(Utility::getCmdOptionWithOrder(argv, argv+argc,"-experiment",3));
            if(allSingle.compare("average")==0)
               exp.saveAverageStructureFunctionCorrelationOfAllSubjects();
            else if(allSingle.compare("single")==0)
               exp.saveStructureFunctionCorrelationOfSingleSubject(std::atoi(Utility::getCmdOption(argv, argv+argc,"-subjectId")));

            //sample: ./brainMatch -experiment subjectwise correlation average -data graph 1 ~/data/PNC/graphs/structureFunctionMappingProject/notNormalized/ -samples ~/data/PNC/src/samples/allSamples.txt -pathType wShortest -funcConn positive -outputPath correlation
            //        ./brainMatch -experiment subjectwise correlation single -subjectId 1 -data graph 1 ~/data/PNC/graphs/structureFunctionMappingProject/notNormalized/ -samples ~/data/PNC/src/samples/allSamples.txt -pathType wShortest -funcConn positive -outputPath correlation
         }
         else if(experimentSubtype.compare("distance" )==0 )
         {
            exp.subjectwiseDistanceExperiment(string(Utility::getCmdOption(argv,argv+argc,"-distanceMeasure")));
            //sample: ./brainMatch -experiment subjectwise distance -data graph 1 ~/data/TBI/graphs/Desikan86/deterministic/ -samples ~/data/TBI/resources/qa_commonList_Deterministic_Desikan86.txt -printSimilarity -funcConn positive -outputPath ./outFile -pathType direct -preprocessGraphs logScaleEDgesStructure_normalizeEdges -distanceMeasure l1
         }
      }
      else if( experimentType.compare("info")==0 )
      {
         Experiment exp(initializerParameters,commandLine);
         exp.setFiles(samplesFile,outputPath,normalizerFile);
         
         if( string(Utility::getCmdOption(argv,argv+argc,"info")).compare("listOfSamples" )==0 )
         {
            exp.listSubjectIDsOfSamplesUsedInExperiment(string(Utility::getCmdOption(argv,argv+argc,"-data")));
            //sample: ./brainMatch -experiment info listOfSamples -data ~/data/PNC/graphs/structureFunctionMappingProject/notNormalized/ -samples ~/data/PNC/src/samples/allSamples.txt
         }
      }
      else
         cerr<<"Incorrect experiment parameter: "<<experimentType<<"  Exiting...\n";
   }
   //</editor-fold>
   else if( Utility::cmdOptionExists( argv, argv+argc, "-utility" ) )
   //<editor-fold defaultstate="collapsed" desc=" Utility functions: -utility <-combine, -thresholdConsistency, -matrixDensity, -meanMatrix, -normalizeDataset, -pick3Columns, -filterResults, -mergeResults, -load/saveGraph(), -shuffle, -diffMatrix">
	{
      string utilityType=string(Utility::getCmdOption( argv, argv+argc, "-utility" ));

      if(utilityType.compare("generateGraph")==0)
      {
         //given a dataset where the structural and functional connectomes are kept in separate files,
         //this function combines them into a single graph file that can be loaded by the matcher program.
         bool printWarning=true;//print warning if data is missing
         Dataset dataset;
         dataset.loadDataset(string(Utility::getCmdOption(argv,argv+argc,"-dti")),string(Utility::getCmdOption(argv,argv+argc,"-fmri")),
                 string(Utility::getCmdOption(argv,argv+argc,"-feature")),string(Utility::getCmdOption(argv,argv+argc,"-samples")),printWarning);
         dataset.saveDataset(string(Utility::getCmdOption(argv,argv+argc,"-outputPath")));
         //sample: ./brainMatch -utility generateGraph -samples ~/data/PNC/src/samples/allSamples.txt -dti ~/data/PNC/dti/streamline/ -fmri ~/data/PNC/fmri/network/ -feature ~/data/PNC/features/ -outputPath ~/data/PNC/graphs/notNormalized/
      }
      else if(utilityType.compare("saveConnectome")==0)
      {
         //given a dataset of graphs, this function produces the structural or functional connectomes of single subject, average of subjects,
         //or all of the subjects as a separate matrix.

         //pathType == pathTransitivity, searchInformation, uCommunicability, wCommunicability, uShortest, wShortest, direct
         //funcConn == positive, negative, full
         string connectomeType = string(Utility::getCmdOptionWithOrder(argv, argv+argc, "saveConnectome",1));//structure or function
         string allSingle = string(Utility::getCmdOptionWithOrder(argv, argv+argc, "saveConnectome",2));// all, single, or average
         if(allSingle.compare("single")==0)
         {
            dataset[0].preprocessGraphs(std::atoi(Utility::getCmdOption(argv, argv+argc,"-subjectId")));
            dataset[0].saveConnectomeOfSingleSubject(std::atoi(Utility::getCmdOption(argv, argv+argc,"-subjectId")),connectomeType,outputPath);
         }
         else if(allSingle.compare("average")==0)
         {
            dataset[0].preprocessGraphs();
            dataset[0].saveAverageConnectomeAcrossSubjects(connectomeType,outputPath);
         }
         else if(allSingle.compare("all")==0)
         {
            dataset[0].preprocessGraphs();
            dataset[0].saveConnectomesOfAllSubjects(connectomeType,outputPath);
         }
            
         //sample: ./brainMatch -utility saveConnectome structure all -data graph 1 ~/dataPath/ -samples ~/dataSamples.txt -pathType wShortest -preprocessGraphs logScaleEdgesStr_normalizeEdges -outputPath connectome
         //        ./brainMatch -utility saveConnectome structure single -subjectId 1  -data matrix 1 -dti ~/dtiPath/ -samples ~/dataSamples.txt -pathType wCommunicability -pathLength 2 -preprocessGraphs logScaleEdgesStr_normalizeEdges -outputPath connectome         
         //        ./brainMatch -utility saveConnectome function single -funcConn full -subjectId 1  -data matrix 1 -fmri ~/fmriPath/ -samples ~/dataSamples.txt -outputPath connectome         
      }
      else if(utilityType.compare("saveVectorizedConnectome")==0)
      {//this function is to be used for printing the node features and edge weights of structural graphs after preprocessing the graphs in a dataset, for the purpose of feeding it into SVM
         Tools::saveGraphsAsVectors(dataset[0],outputPath);
         //sample: ./brainMatch -utility saveVectorizedConnectome -data graph 1 ~/data/TBI/graphs/Desikan86/deterministic/ -pathType wCommunicability -pathLength 2 -preprocessGraphs logScaleEdgesStr_normalizeEdges -outputPath vec.txt 
      }
      else if(utilityType.compare("thresholdConsistency")==0)
      {
         //given a dataset where the weighted structural connectomes are kept in individual files, 
         //this function thresholds each connectome using consistency thresholding and saves thresholded connectomes to the target folder
         //consistency thresholding: across the dataset, the edges that are consistently present are kept in each subject's connectome
         Tools::thresholdConsistency(string(Utility::getCmdOption(argv,argv+argc,"-controls")),string(Utility::getCmdOption(argv,argv+argc,"-allSubjects")),
                 string(Utility::getCmdOption(argv,argv+argc,"-connectomes")),std::stof(Utility::getCmdOption(argv,argv+argc,"-density")),string(Utility::getCmdOption(argv,argv+argc,"-outputPath")));
         
         //sample: ./brainMatch -utility thresholdConsistency -controls ~/data/PNC/src/samples/allSamples.txt -allSubjects ~/data/PNC/src/samples/allSamples.txt -connectomes ~/data/PNC/dti/streamline/ -density 0.15 -outputPath /home/yusuf/data/PNC/dti/streamlineThresholded/
      }
      else if(utilityType.compare("thresholdDensity")==0)
      {
         if((string(Utility::getCmdOption(argv,argv+argc,"thresholdDensity"))).compare("group")==0)
         {
            //given a dataset where the weighted structural connectomes are kept in individual files, 
            //this function thresholds each connectome by establishing the mean conenctome of the healthy controls to be at a certain density
            //and saves thresholded connectomes to the target folder
            //density thresholding: across the healthy controls, the mean connectome is set to have a certain density (i.e., the threshold is determined to keep the density of the mean connectome at a certain level)
            Tools::thresholdDensityByGroup(string(Utility::getCmdOption(argv,argv+argc,"-controls")),string(Utility::getCmdOption(argv,argv+argc,"-allSubjects")),
                    string(Utility::getCmdOption(argv,argv+argc,"-connectomes")),std::stof(Utility::getCmdOption(argv,argv+argc,"-density")),string(Utility::getCmdOption(argv,argv+argc,"-outputPath")));

            //sample: ./brainMatch -utility thresholdDensity group -controls ~/data/PNC/src/samples/allSamples.txt -allSubjects ~/data/PNC/src/samples/allSamples.txt -connectomes ~/data/PNC/dti/streamline/ -density 0.15 -outputPath /home/yusuf/data/PNC/dti/streamlineThresholded/
         }
         else if((string(Utility::getCmdOption(argv,argv+argc,"thresholdDensity"))).compare("perSubject")==0)
         {
            //given a dataset where the weighted structural connectomes are kept in individual files, 
            //this function thresholds each connectome by establishing the mean conenctome of the healthy controls to be at a certain density
            //and saves thresholded connectomes to the target folder
            //density thresholding: across the healthy controls, the mean connectome is set to have a certain density (i.e., the threshold is determined to keep the density of the mean connectome at a certain level)
            Tools::thresholdDensityPerSubject(string(Utility::getCmdOption(argv,argv+argc,"-allSubjects")),string(Utility::getCmdOption(argv,argv+argc,"-connectomes")),
                    std::stof(Utility::getCmdOption(argv,argv+argc,"-density")),string(Utility::getCmdOption(argv,argv+argc,"-outputPath")));

            //sample: ./brainMatch -utility thresholdDensity perSubject -allSubjects ~/data/PNC/src/samples/allSamples.txt -connectomes ~/data/PNC/dti/streamline/ -density 0.15 -outputPath /home/yusuf/data/PNC/dti/streamlineThresholded/
         }
      }
      else if(utilityType.compare("thresholdMatrices")==0)
      {
         //given a dataset of matrices, this function produces the binarized version of the matrices taking the provided threshold as the reference point.
         //if the provided matrix contains positive and negative values, the threshold will be treated as +-threshold, and the numbers above +threshold will be set to 1
         //the values lewss than -threshold will be set to -1, and the rest will be set to 0.
         Tools::thresholdMatrices(string(Utility::getCmdOption(argv,argv+argc,"-allSubjects")),string(Utility::getCmdOption(argv,argv+argc,"-connectomes")),
                 string(Utility::getCmdOption(argv,argv+argc,"-outputPath")),std::stof(Utility::getCmdOption(argv,argv+argc,"-threshold")));
         
         //sample: ./brainMatch -utility thresholdMatrices -allSubjects ~/data/PNC/src/samples/allSamples.txt -connectomes ~/data/PNC/dti/streamline/ -threshold 0.15 -outputPath /home/yusuf/data/PNC/dti/streamlineThresholded/
      }
      else if(utilityType.compare("binarizeMatrices")==0)
      {
         //given a dataset of matrices, this function produces the binarized version of the matrices taking the provided threshold as the reference point.
         //if the provided matrix contains positive and negative values, the threshold will be treated as +-threshold, and the numbers above +threshold will be set to 1
         //the values lewss than -threshold will be set to -1, and the rest will be set to 0.
         Tools::binarizeMatrices(string(Utility::getCmdOption(argv,argv+argc,"-allSubjects")),string(Utility::getCmdOption(argv,argv+argc,"-connectomes")),
                 string(Utility::getCmdOption(argv,argv+argc,"-outputPath")),std::stof(Utility::getCmdOption(argv,argv+argc,"-threshold")));
         
         //sample: ./brainMatch -utility binarizeMatrices -allSubjects ~/data/PNC/src/samples/allSamples.txt -connectomes ~/data/PNC/dti/streamline/ -threshold 0.15 -outputPath /home/yusuf/data/PNC/dti/streamlineThresholded/
      }
      else if(utilityType.compare("matrixDensity")==0)
      {
         //given a dataset where the weighted structural connectomes are kept in individual files, 
         //this function calculates the density of the connectome in terms of number of nonzero weighted edges
         Tools::calculateDensityOfMatrices(string(Utility::getCmdOption(argv,argv+argc,"-samples")),string(Utility::getCmdOption(argv,argv+argc,"-folderPath")));
         
         //sample: ./brainMatch -utility matrixDensity -folderPath ~/data/PNC/dti/streamline/ 
         //        ./brainMatch -utility matrixDensity -samples ~/data/PNC/src/samples/samples_10.txt -folderPath ~/data/PNC/dti/streamlineThresholded/
      }
      else if(utilityType.compare("matrixConnectivityStrength")==0)
      {
         //given a dataset where the weighted structural connectomes are kept in individual files, 
         //this function calculates the density of the connectome in terms of number of nonzero weighted edges
         Tools::calculateStrengthOfMatrices(string(Utility::getCmdOption(argv,argv+argc,"-samples")),string(Utility::getCmdOption(argv,argv+argc,"-folderPath")),string(Utility::getCmdOption(argv,argv+argc,"-hemispheres")));
         
         //sample: ./brainMatch -utility matrixConnectivityStrength -folderPath ~/data/PNC/dti/streamline/ -hemispheres ~/data/PNC/src/samples/hemispheres.txt
         //        ./brainMatch -utility matrixConnectivityStrength -samples ~/data/PNC/src/samples/samples_10.txt -folderPath ~/data/PNC/dti/streamlineThresholded/ -hemispheres ~/data/PNC/src/samples/hemispheres.txt
      }
      else if(utilityType.compare("graphConnectivity")==0)
      {
         //given a dataset where the weighted structural connectomes are kept in individual files, 
         //this function calculates the number of connected components in a graph
         //@type: <graph> or <matrix>
         Tools::checkGraphConnectivity(string(Utility::getCmdOption(argv,argv+argc,"-connectivity")),string(Utility::getCmdOption(argv,argv+argc,"-connectomePath")),Utility::cmdOptionExists(argv,argv+argc,"-printComponents"),string(Utility::getCmdOption(argv,argv+argc,"-type")),string(Utility::getCmdOption(argv,argv+argc,"-samples")));
         
        
         //sample: ./brainMatch -utility graphConnectivity -connectomePath ~/data/PNC/dti/streamline/ -connectivity structure -type matrix -printComponents
         //        ./brainMatch -utility graphConnectivity -samples ~/data/PNC/src/samples/samples_10.txt -connectomePath ~/data/PNC/dti/streamlineThresholded/ -connectivity structure
      }
      else if(utilityType.compare("missingNode")==0)
      {
         //given a dataset where the weighted structural connectomes are kept in individual files, 
         //this function calculates the number of connected components in a graph
         //@type: <graph> or <matrix>
         Tools::checkMissingNode(string(Utility::getCmdOption(argv,argv+argc,"-connectomePath")),std::atoi(Utility::getCmdOption(argv,argv+argc,"-numNodes")),string(Utility::getCmdOption(argv,argv+argc,"-type")),string(Utility::getCmdOption(argv,argv+argc,"-samples")));
         
        
         //sample: ./brainMatch -utility missingNode -connectomePath ~/data/PNC/dti/streamline/ -type matrix -numNodes 118
         //        ./brainMatch -utility missingNode -samples ~/data/PNC/src/samples/samples_10.txt -connectomePath ~/data/PNC/dti/streamline/ -type matrix -numNodes 118
      }
      else if(utilityType.compare("reduceConnectome")==0)
      {
         //given a folder where the connectomes are kept in individual files, along with a vector that maps the regions of this connectome into functional groups
         //this function reduces the connectomes into functional systems and saves them in a separate file
         Tools::reduceMatricesIntoFunctionalNetworks(string(Utility::getCmdOption(argv,argv+argc,"-connectomePath")),string(Utility::getCmdOption(argv,argv+argc,"-mappingFile")),string(Utility::getCmdOption(argv,argv+argc,"-outputFolder")),Utility::cmdOptionExists(argv, argv+argc, "-normalize"));
         
         //sample: ./brainMatch -utility reduceConnectome -connectomePath ~/repo/projects/tbiStructureAndFunction/results/2_groupDifference_systems_thresh15_Desikan86/connectomes/direct/ -mappingFile ~/data/atlases/Desikan86/desikan86_function.csv -outputFolder ./connectome/ -normalize
      }
      else if(utilityType.compare("meanMatrix")==0)
      {
         //given a folder containing matrices of same size, this function calculates the mean of the matrices and saves it to a file
         int skipLine = 0;//number of lines to skip before starting to read the matrix
         int row = 0, col=0;//number of rows and columns to read after skipline (practically this is the size of the matrix to read, if set to 0 will read whatever is in the file)
         bool append=false;
         Tools::calculateMeanMatrix(string(Utility::getCmdOption(argv,argv+argc,"-matricesPath")),string(Utility::getCmdOption(argv,argv+argc,"-outputPath")),append,row,col,skipLine);
         //sample: ./brainMatch -utility meanMatrix -matricesPath ~/data/PNC/dti/streamline/ -outputPath out.txt
      }
      else if(utilityType.compare("symmetrizeMatrices")==0)
      {
         //given a folder containing matrices of same size, this function symmetrizes the matrices by taking 
         //the average of upper/lower triangular and saves it as a separate file
         Tools::symmetrizeMatrices(string(Utility::getCmdOption(argv,argv+argc,"-matricesPath")),string(Utility::getCmdOption(argv,argv+argc,"-outputPath")));
         //sample: ./brainMatch -utility symmetrizeMatrices -matricesPath ~/data/PNC/dti/streamline/ -outputPath ./out/
      }
      else if(utilityType.compare("scaleUpMatrices")==0)
      {
         //given a folder containing matrices of same size with nonnegative values, this function scales up values of each matrix
         //to make the smallest nonzero value equal to 1
         Tools::scaleUpMatricesMakingMinNonZeroEntryUnity(string(Utility::getCmdOption(argv,argv+argc,"-matricesPath")),string(Utility::getCmdOption(argv,argv+argc,"-outputPath")));
         //sample: ./brainMatch -utility scaleUpMatrices -matricesPath ~/data/PNC/dti/streamline/ -outputPath ./out/
      }
      else if(utilityType.compare("modifyMatrices")==0)
      {
         //a generic function to make changes to matrices contained under a certain folder. Modify the code in this function to do desired changes to matrices...
         Tools::modifyMatrices(string(Utility::getCmdOption(argv,argv+argc,"-matricesPath")),string(Utility::getCmdOption(argv,argv+argc,"-outputPath")));
         //sample: ./brainMatch -utility modifyMatrices -matricesPath ~/data/PNC/dti/streamline/ -outputPath ./out/
      }
      else if(utilityType.compare("subtractMatrices")==0)
      {
         Tools::subtractMatrices(string(Utility::getCmdOptionWithOrder(argv,argv+argc,"-inFile",1)),string(Utility::getCmdOptionWithOrder(argv,argv+argc,"-inFile",2)),string(Utility::getCmdOption(argv,argv+argc,"-outputPath")));
         //sample: ./brainMatch -utility diffMatrices -inFile matrix1.txt matrix2.txt -outputPath out.txt
      }
      else if(utilityType.compare("normalizeDataset")==0)
      {
         //given a dataset, calculates the ranges for each feature type that is to be used in normalization
         //later on, and saves these ranges into file
         dataset[0].normalizeDataset();  
         dataset[0].saveNormalizer(outputPath);
         
         //sample: ./brainMatch -utility normalizeDataset -graphs ~/data/PNC/graphs/ -samples ~/data/PNC/src/samples70.txt -outputPath normalizedRanges.txt
      }
      else if(utilityType.compare("mergeFiles")==0)
      {
         //given a list of files that contain results of the same experiment, this function merges them into a single file
         Utility::mergeFiles(string(Utility::getCmdOption(argv, argv+argc,"-folderPath")),string(Utility::getCmdOption(argv, argv+argc,"-outputPath")),string(Utility::getCmdOption(argv, argv+argc,"-list")));
         
         //sample input: ./brainMatch -utility mergeFiles -folderPath ./test3/ -list ./test3/test3List.txt -outputPath ./test3/test3Results.txt
      }
      else
      {
         cerr<<"Incorrect utility parameter: "<<utilityType<<"  Exiting...\n";
         exit(1);
      }
   }
   //</editor-fold>
   else if( Utility::cmdOptionExists( argv, argv+argc, "-test" ) )
   //<editor-fold defaultstate="collapsed" desc=" Test code: -test <loadGraph, load/saveGraph(), shuffle, random, communicability, sandbox, matcher>">
	{
      string testType=string(Utility::getCmdOption( argv, argv+argc, "-test" ));
      if( testType.compare("loadGraph" )==0 )
      {
         //this is to be used for testing the load/save function of the graph class
         Graph graph;
         if( Utility::cmdOptionExists( argv, argv+argc, "-graph" ) )
            graph.loadGraph(string(Utility::getCmdOption(argv,argv+argc,"-graph")));
         if( Utility::cmdOptionExists( argv, argv+argc, "-save" ) )
            graph.saveGraph(string(Utility::getCmdOption(argv,argv+argc,"-save")));
         
         //sample: ./brainMatch -test loadGraph -input ~/repo/results/PNC/2_referencePoint/exp2/QAPD.txt -output max.txt -columns 7 0 15
         //        ./brainMatch -test loadGraph -input ~/repo/results/PNC/2_referencePoint/exp2/QAPD.txt -output avg.txt -columns 7 0 14
      }
      else if( testType.compare("shuffle")==0 )
      {
         Test::shuffleGraph(string(Utility::getCmdOption(argv,argv+argc,"-inFile")),string(Utility::getCmdOption(argv,argv+argc,"-outputPath")),std::atoi(Utility::getCmdOption(argv,argv+argc,"-seed")),std::atoi(Utility::getCmdOption(argv,argv+argc,"-iteration")), Utility::cmdOptionExists( argv, argv+argc, "-function" ));
         //sample: ./brainMatch -test shuffle -inFile in.grp -outputPath out.grp -seed 10 -iteration 10
      }
      else if( testType.compare("random")==0 )
      {
         Test::testRandom(std::stof(Utility::getCmdOption(argv, argv+argc,"-min")),std::stof(Utility::getCmdOption(argv, argv+argc,"-max")),std::atoi(Utility::getCmdOption(argv, argv+argc,"-count")),std::atoi(Utility::getCmdOption(argv, argv+argc,"-bin")));
         //sample: ./brainMatch -test random -min 0 -max 10 -count 150 -bin 35
      }
      else if( testType.compare("communicability")==0 )
      {
         Test::testCommunicability(std::atoi(Utility::getCmdOption(argv, argv+argc,"-numNodes")),string(Utility::getCmdOption(argv,argv+argc,"-aMatrix")),string(Utility::getCmdOption(argv,argv+argc,"-aRefMatrix")),string(Utility::getCmdOption(argv,argv+argc,"-diff")));
         //sample: ./brainMatch -test communicability -numNodes 250 -aMatrix matrix1.txt -aRefMatrix matrix2.txt -diffMatrix diff.txt
      }
      else if( testType.compare("sandbox")==0 )
      {
         Test::testCode();
         //sample: ./brainMatch -test sandbox
      }
      //<editor-fold defaultstate="collapsed" desc=" Tests for matcher: -match <datasetFull, datasetPartial, twoGraphs>">
      else if( Utility::cmdOptionExists( argv, argv+argc, "matcher" ) ) 
      {
         long time[4];

         time[0] = Utility::getTime();
         Dataset dataset(initializerParameters,string(Utility::getCmdOption(argv,argv+argc,"-graphs")),samplesFile);
         time[1] = Utility::getTime();
         dataset.preprocessGraphs();
         time[2] = Utility::getTime();
         Matcher matcher(&dataset);

         string matcherType=string(Utility::getCmdOption(argv, argv+argc,"matcher"));
         
         if(matcherType.compare("datasetPartial")==0)
         {
            int rowStart = std::atoi(Utility::getCmdOption(argv, argv+argc,"-rowStart"));
            int rowEnd = std::atoi(Utility::getCmdOption(argv, argv+argc,"-rowEnd"));
            int columnStart = std::atoi(Utility::getCmdOption(argv, argv+argc,"-columnStart"));
            int columnEnd = std::atoi(Utility::getCmdOption(argv, argv+argc,"-columnEnd"));

            matcher.matchEveryoneWithEveryone(outputPath,rowStart,rowEnd,columnStart,columnEnd);
            //sample: ./brainMatch -test matcher datasetPartial -graphs ~/data/TBI/graphs/ -samples QA_pass_subject_list_s1 -solver MLQP  -adjustMode none -assignmentCostMode volume -separationCostMode functional -outputPath results.txt -rowStart 0 -rowEnd 4 -columnStart 0 -columnEnd 4
         }
         else if(matcherType.compare("datasetFull")==0)
         {
            matcher.matchEveryoneWithEveryone(outputPath,printMatches);
            //sample: ./brainMatch -test matcher datasetFull -graphs ~/data/TBI/graphs/ -samples QA_pass_subject_list_s1 -solver MLPD  -adjustMode none -assignmentCostMode volume -separationCostMode functional -out results.txt
         }
         else if(matcherType.compare("graphPair")==0) 
         {
            int graph1 = std::atoi(Utility::getCmdOption(argv, argv+argc,"-graph1"));
            int graph2 = std::atoi(Utility::getCmdOption(argv, argv+argc,"-graph2"));

            matcher.matchTwoGraphs(graph1,graph2,printMatches);
            time[3] = Utility::getTime();

            cerr<<"loading time:"<<time[1]-time[0]<<"\t normalize time:"<<time[2]-time[1]<<"\t matchTime:"<<time[3]-time[2]<<endl;
            //sample: ./brainMatch -test matcher graphPair -graphs ~/data/PNC/graphs/structureFunctionMappingProject/notNormalized/ -samples ~/data/PNC/src/samples/allSamples.txt -graph1 0 -graph2 0 -pathType direct -printMatches -printSimilarity
         }
      }
   //</editor-fold>
   }
   //</editor-fold>
   else
   {
      cerr<<"Command line couldn't be parsed properly. Check your input command. To see possible commands, use -help <experiments|utility>"<<endl;
      exit(1);
   }
}

void commandLineHelp(std::string subcategory)
{
   cout<<"\n______Command Line Help______"<<endl;
   cout<<"\t DisplayHelp: ./brainMatch -help <experiments|utility>"<<endl;
   
   if(subcategory.compare("")==0 || subcategory.compare("experiments")==0)
   {
      cout<<"\n\t[[While Loading data]] : <DATA> =  -data <graph 1 <~/graphFolder/> > | <matrix 1 -dti </dtiPath/> -fmri </fmriPath/> -features </featuresPath/> >"<<endl;
      cout<<"\n\t[[While entering funcConn]] : if -modality is str_func, then -funcConn can be positive or negative, whereas if -modality is func_func, then -funcConn can be positive, negative, or full."<<endl;
      cout<<"\n\t[[While doing permutation test]] : -shuffle parameter needs to take <shuffleType> and <shuffleDatasetId> as the first and second parameters, where the first and the second needs to be <structure | function> and <1|2> respectively."<<endl;
      
      cout<<"\n______Experiments______"<<endl;
      cout<<"\n\t______Groupwise______"<<endl;
      cout<<"\t\t __match: \n\t\t\t./brainMatch -experiment groupwise match -groups c p <-permutation 10 <-seed 10 -shuffle <structure | function> <1 | 2> > > -modality <str_str | str_func | func_func> -funcConn <positive | negative | full> <DATA> -samples <~/samples.txt> -printMatches -printSimilarity -outputPath <~/outputFolder/> -pathType direct"<<endl;
      cout<<"\t\t __distance: \n\t\t\t./brainMatch -experiment groupwise distance -distanceMeasure <subtract | l1 | l2> -groups c all -modality <str_str | str_func | func_func> -funcConn <positive | negative | full> <DATA> -samples <~/samples.txt> -printSimilarity -outputPath <~/outputFolder/> -pathType direct"<<endl;
      cout<<"\t\t __subSpecZScore: \n\t\t\t./brainMatch -experiment groupwise subSpecZScore -groups c all -modality <str_str | str_func | func_func> -funcConn <positive | negative | full> <DATA> -samples <~/samples.txt> -outputPath <~/outputFolder/> -pathType direct"<<endl;
      cout<<"\t\t __subSpecDifference: \n\t\t\t./brainMatch -experiment groupwise subSpecDifference -groups c all -modality <str_str | str_func | func_func> -funcConn <positive | negative | full> <DATA> -samples <~/samples.txt> -outputPath <~/outputFolder/> -pathType direct"<<endl;
      cout<<"\t\t __averageResults: \n\t\t\t./brainMatch -experiment groupwise average -groups c all -resultsFolder <~/inResultsFolder/> -samples <~/samples.txt> -outputPath <~/outputFolder/>"<<endl;

      cout<<"\n\t______Subjectwise______"<<endl;
      cout<<"\t\t __match: \n\t\t\t./brainMatch -experiment subjectwise match <-permutation 333 <-seed 10 -shuffle <structure | function> <1 | 2> > > -modality <str_str | str_func | func_func> -funcConn <positive | negative | full> <DATA> -samples <~/samples.txt> -printMatches -printSimilarity -outputPath <~/outputFolder/> -pathType wCommunicability <-pathLength 3> <-saveAvgConnectome>"<<endl;
      cout<<"\t\t __distance: \n\t\t\t./brainMatch -experiment subjectwise distance -modality <str_str | str_func | func_func> -funcConn <positive | negative> <DATA> -samples <~/samples.txt> -printMatches -printSimilarity -outputPath <~/outputFolder/> -pathType wShortest <-saveAvgConnectome>"<<endl;
      cout<<"\t\t __correlation: \n\t\t\t./brainMatch -experiment subjectwise correlation <single -subjectId 1 | average> <DATA> -samples <~/samples.txt> -pathType wShortest -funcConn <positive | negative> -outputPath <~/outputFolder/>"<<endl;

      
      
      cout<<"\n\t______Info______"<<endl;
      cout<<"\t\t __SampleGraphs: \n\t\t\t./brainMatch -experiment info sampleGraphs -graphs <~/graphFolder/> -samples <~/samples.txt> -outputPath <~/outputFolder/> -pathType wShortest"<<endl;
      cout<<"\t\t __ListOfSamples: \n\t\t\t./brainMatch -experiment info listOfSamples -graphs <~/graphFolder/> -samples <~/samples.txt>"<<endl;
   }
   if(subcategory.compare("")==0 || subcategory.compare("utility")==0)
   {
      cout<<"\n______Utility______"<<endl;
      cout<<"\t __GenerateGraph: \n\t\t./brainMatch -utility generateGraph -samples <~/samples.txt> < -dti <~/dtiConnectomeFolder/> | -fmri <~/fmriConnectomeFolder/> | -feature <~/nodeFeaturesFolder/> > -outputPath <~/outputFolder/>"<<endl;
      cout<<"\t __SaveConnectome: \n\t\t./brainMatch -utility saveConnectome <structure | function> <all | single -subjectId 1 | average> <DATA> -samples <~/samples.txt> -outputPath <~/outputFolder/> <-pathType wCommunicability -pathLength 2 | -funcConn positive|negative|full > <-preprocessGraphs logScaleEdgesStructure_traffic_normalizeAll>"<<endl;
      cout<<"\t __SaveVectorizedConnectome: \n\t\t./brainMatch -utility saveVectorizedConnectome <DATA> -samples <~/samples.txt> -outputPath <~/outputFolder/> <-pathType wCommunicability -pathLength 2> <-preprocessGraphs logScaleEdgesStructure_traffic_normalizeEdges>"<<endl;
      cout<<"\t __ThresholdConsistency: \n\t\t./brainMatch -utility thresholdConsistency -controls <~/controls.txt> -allSubjects <~/allSubjects.txt> -connectomes <~/connectomesFolder/> -density 0.15 -outputPath <~/outputFolder/>"<<endl;
      cout<<"\t __ThresholdDensity: \n\t\t./brainMatch -utility thresholdDensity <group | perSubject> < -controls <~/controls.txt> > -allSubjects <~/allSubjects.txt> -connectomes <~/connectomesFolder/> -density 0.15 -outputPath <~/outputFolder/>"<<endl;
      cout<<"\t __ThresholdMatrices: \n\t\t./brainMatch -utility thresholdMatrices -allSubjects <~/allSubjects.txt> -connectomes <~/connectomesFolder/> -threshold 5 -outputPath <~/outputFolder/>"<<endl;
      cout<<"\t __BinarizeMatrices: \n\t\t./brainMatch -utility binarizeMatrices -allSubjects <~/allSubjects.txt> -connectomes <~/connectomesFolder/> -threshold 5 -outputPath <~/outputFolder/>"<<endl;
      cout<<"\t __MatrixDensity: \n\t\t./brainMatch -utility matrixDensity -folderPath <~/matricesFolder/> <-samples <~/samples.txt> >"<<endl;
      cout<<"\t __GraphConnectivity: \n\t\t./brainMatch -utility graphConnectivity -connectomePath <~/connectomeFolder/> <-samples <~/samples.txt> > -connectivity <structure|positiveFunction|negativeFunction> -type <matrix|graph> <-printComponents>"<<endl;
      cout<<"\t __ReduceConnectome: \n\t\t./brainMatch -utility reduceConnectome -connectomePath <~/connectomeFolder/> -mappingFile <~/networkMappingFile.csv> -outputPath <~/outputFolder/> <-normalize>"<<endl;
      cout<<"\t __MeanMatrix: \n\t\t./brainMatch -utility meanMatrix -matricesPath <~/matricesFolder/> -outputPath <~/outputFile.txt>"<<endl;
      cout<<"\t __SymmetrizeMatrices: \n\t\t./brainMatch -utility symmetrizeMatrices -matricesPath <~/matricesFolder/> -outputPath <~/outputFolder>"<<endl;
      cout<<"\t __ScaleUpMatrices: \n\t\t./brainMatch -utility scaleUpMatrices -matricesPath <~/matricesFolder/> -outputPath <~/outFile.txt>"<<endl;
      cout<<"\t __SubtractMatrices: \n\t\t./brainMatch -utility subtractMatrices -inFile <~/file1.txt> <~/file1.txt> -outputPath <~/outFile.txt>"<<endl;
      cout<<"\t __MergeFiles: \n\t\t./brainMatch -utility mergeFiles -folderPath <~/inputFilesPath/> <-list <~/listOfFilenamesToMerge.txt> > -outputPath <~/mergedFile.txt>"<<endl;
   }
}