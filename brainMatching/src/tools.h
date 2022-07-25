/* 
 * File:   tools.h
 * Author: yusuf
 *
 * Generated on May 10, 2017, 1:08 PM
 */

#ifndef TOOLS_H
#define TOOLS_H

#include <iostream>
#include <sstream>
#include <fstream>
#include "utility.h"
#include "geometry.h"
#include "dataset.h"

namespace Tools
{
	//applies consistency thresholding for a given connectivity matrix
	void thresholdConsistency(std::string controlsSubjectList,std::string allSubjectList, std::string streamlineFolder, float density, std::string outputFolder);
	void thresholdDensityByGroup(std::string controlsSubjectList,std::string allSubjectList, std::string streamlineFolder, float density, std::string outputFolder);
	void thresholdDensityPerSubject(std::string allSubjectList, std::string streamlineFolder, float density, std::string outputFolder);
	void thresholdMatrices(std::string allSubjectList, std::string streamlineFolder, std::string outputFolder, float threshold);

	void binarizeMatrices(std::string allSubjectList, std::string matrixFolder, std::string outputFolder, float threshold);

	void calculateDensityOfMatrices(std::string subjectList,std::string matrixFolder);
	void calculateStrengthOfMatrices(std::string subjectList,std::string matrixFolder,std::string hemisphereFile="");
	void checkGraphConnectivity(std::string featureName, std::string graphsFolder,bool printComponents=false,std::string type="graph",std::string samplesFile="");
	void checkMissingNode(std::string graphsFolder,int numNodes,std::string type="graph",std::string samplesFile="");
	
	void calculateMeanVector(std::string vectorsFolder, std::string outFile, bool append=false, int size=0, int skipLine=0);
	void calculateMeanMatrix(std::string matricesFolder, std::string outFile, bool append=false, int _row=0, int _col=0, int skipLine=0);
	void symmetrizeMatrices(std::string matricesFolder, std::string outFolder);
	void scaleUpMatricesMakingMinNonZeroEntryUnity(std::string matricesFolder, std::string outFolder);
	void subtractMatrices(std::string inFile1, std::string inFile2, std::string outFile);
	void reduceMatricesIntoFunctionalNetworks(std::string resultsFolder, std::string functionalMappingsFile, std::string outfile, bool normalize=false);
	void saveGraphsAsVectors(Dataset &dataset, std::string outputFile);
	void modifyMatrices(std::string matricesFolder, std::string outFolder);
}

#endif /* TOOLS_H */

