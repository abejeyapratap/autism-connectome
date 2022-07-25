/* 
 * File:   experiment.cpp
 * Author: yusuf
 *
 * Generated on September 21, 2016, 5:35 PM
 */

#include <complex>
#include <math.h>
#include "experiment.h"
#include "dataset.h"
#include "graph.h"
#include "matcher.h"
#include "utility.h"
#include "tools.h"
//#include <dlib/matrix.h>

using namespace std;

void Experiment::setFiles(string _samplesFile, string _outputPath, string _normalizerFile)
{
   samplesFile = _samplesFile;
   outputPath = _outputPath;
   normalizerFile = _normalizerFile;
}

void Experiment::setDataset(Dataset *_dataset,Dataset *_dataset2)
{
   dataset = _dataset;
   dataset2 = _dataset2;
}

void Experiment::setModality(std::string modality)
{
   edgeType1 = Edge::getEdgeType(Utility::splitString(modality,"_",0));//modality should look like: str_func, str_str, func_func
   edgeType2 = Edge::getEdgeType(Utility::splitString(modality,"_",1));//it can also be in the form of: structure_structure, structure_function, function_function
}

void Experiment::setPrint(bool _printMatches,bool _printSimilarityScores,bool _printAverageSimilarityScores, bool _printMatchingAccuracies, bool _printAverageMatchingAccuracies)
{
   printMatches = _printMatches;
   printSimilarityScores = _printSimilarityScores;
   printAverageSimilarityScores = _printAverageSimilarityScores;
   printMatchingAccuracies = _printMatchingAccuracies;
   printAverageMatchingAccuracies = _printAverageMatchingAccuracies;
}

//For an experiment, this function prints the names of the subjects that has a graph file associated with them
//Note: there might be cases where the subject is included in the sample set while 
//      graph in a certain parcellation does not exist for this subject.
void Experiment::listSubjectIDsOfSamplesUsedInExperiment(std::string dataFolder)
{
   vector<string> filenames,subjectNames;
   Utility::loadFileNamesFromList(samplesFile,subjectNames);
   Utility::loadFilePathFromFolderSelectivelyPreservingOrder(samplesFile,dataFolder,filenames);

   for(vector<string>::iterator filenameIter=filenames.begin();filenameIter!=filenames.end();filenameIter++)
   {
      if((*filenameIter).compare("")!=0)
      {
         //prune the subject name from the path
         // /home/yusuf/data/SmithDOS/graphs/Desikan86/thresh15/DOS015.grp
         cout<<Utility::pruneFilenameFromPath(*filenameIter)<<endl;
      }
   }
}

void Experiment::identifyGroups(std::string samplesFilename,string controlText, string patientText, std::vector<int> &controlsOrder, std::vector<int> &patientsOrder,std::vector<std::string> &controlIDs, std::vector<std::string> &patientIDs)
{
   controlsOrder.clear();
   patientsOrder.clear();
   controlIDs.clear();
   patientIDs.clear();
   
   int subjectOrder=0;//id in the order of the graph loaded to memory
   
   std::ifstream file; 
   std::istringstream is; 
   std::string line, subjectID, source;

   ////load names of result files
   file.open(samplesFilename.c_str());
   if(file.fail())
   {
      std::cout << "File not found: "<<samplesFilename<<" !!!!\n";
      return;
   }	

   while(getline(file, line))
   {
      is.clear();
      is.str(line);
      is >> subjectID;
      
      if(subjectID.find(controlText)!=std::string::npos || controlText.compare("all")==0)
      {
         controlsOrder.push_back(subjectOrder);
         controlIDs.push_back(subjectID);
      }
      if(subjectID.find(patientText)!=std::string::npos || patientText.compare("all")==0)
      {
         patientsOrder.push_back(subjectOrder);
         patientIDs.push_back(subjectID);
      }
      subjectOrder++;
   }
   if(patientsOrder.size()==0 || controlsOrder.size()==0)
   {
      cerr<<"Encountered some issues with subject groups, as the number of patients:"<<patientsOrder.size()<<" and controls:"<<controlsOrder.size()<<endl;
      cerr<<"Most probably, this is due to incorrect subject groups: you entered control:<"<<controlText<<"> and patient:<"<<patientText<<">"<<endl;
      cerr<<"Error in Experiment::identifyGroups() function. Exiting!!!\n";
      exit(1);
   }
   file.close();
}
