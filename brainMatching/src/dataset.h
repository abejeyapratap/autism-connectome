/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   dataset.h
 * Author: yusuf
 *
 * Created on April 10, 2019, 3:13 PM
 */

#ifndef DATASET_H
#define DATASET_H

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <algorithm> // std::find
#include <dirent.h>  //for readdir, opendir, closedir
#include "graph.h"
#include "edge.h"
#include "normalizer.h"
#include "parameterSet.h"
#include "initializerParameters.h"

class Dataset
{
	public:
		Dataset(const Dataset &dataset);
		Dataset(InitializerParameters &_initializerParameters,std::string folderPath, std::string objectNames);
		Dataset(InitializerParameters &_initializerParameters,std::string streamlineFolder, std::string fmriNetworkFolder, std::string extraFeaturesFolder, std::string objectNames);
		Dataset(){}
		~Dataset(){}

		//calculates path between nodes and scales/normalizes nodes/edges
		void preprocessGraphs(int order=-1);
		//shuffles graph edges in structure/function
		void shuffleGraphs(int seedSupplement, std::string type);
		
		///given a dataset, save connectomes in matrix form into file
		void saveConnectomeOfSingleSubject(int subjectOrder, std::string connectomeType, std::string outputPath);
		void saveAverageConnectomeAcrossSubjects(std::string connectomeType, std::string outputPath);
		void saveConnectomesOfAllSubjects(std::string connectomeType, std::string outputPath);
		
		//given a dataset, retrieve the connectome of one of the graphs
		void getConnectivityMatrixForGraph(int graphOrder, Edge::Feature edgeType, float **matrix, std::string functionalConnectivity="", bool nonnegativeConnectivity=true);
		
		//given a dataset, save one of the graphs into file in .grp format
		inline void saveGraph(int i, std::string filename){graphs[i].saveGraph(filename);}//for debug
		
		//normalize graphs in the dataset wrt graphs of other subjects
		inline void normalizeDataset(){normalizer.normalizeData(graphs);}
		inline void saveNormalizer(std::string filename){normalizer.save(filename);}
		inline void loadNormalizer(std::string filename){normalizer.load(filename);}
		inline Normalizer& getNormalizer(){return normalizer;}
		
		inline int getNumOfSubjects(){return graphs.size();}
		static int getNumOfSubjects(std::string folderPath, std::string objectNames);
		inline int getSizeOfAGraph(){return graphs.begin()->second.getNumNodes();}
		static int getSizeOfAGraph(std::string folderPath, std::string objectNames);
		
		inline void setInitializerParameters(InitializerParameters &_initializerParameters){initializerParameters = _initializerParameters;}
		inline InitializerParameters& getInitializerParameters(){return initializerParameters;}
		inline void printParameters(){initializerParameters.print("all");}
      
		inline std::map<int,Graph>::iterator getGraphIteratorBegin(){return graphs.begin();}
		inline std::map<int,Graph>::iterator getGraphIteratorEnd(){return graphs.end();}
		inline Graph& getGraph(int graphId){return graphs[graphId];}
			
		
		//load dataset from graph files
		void loadDataset(std::string folderPath, std::string objectNames, bool printWarning=false);
		void loadDataset(std::string folderPath, std::string objectNames, int rowStart, int rowEnd, int columnStart,int columnEnd, bool printWarning=false);
		void loadDatasetReduced(std::string folderPath, std::string objectNames, int graphToDrop, bool printWarning=false);
		//load dataset from separate connectome matrices
		void loadDataset(std::string streamlineFolder, std::string fmriNetworkFolder, std::string extraFeaturesFolder, std::string subjectList, bool printWarning=false);
      
		//save the graphs stored inside the dataset into a folder in the .grp format
		void saveDataset(std::string outputFolder);
		
	private:
		std::map<int,Graph> graphs;
		InitializerParameters initializerParameters;
		Normalizer normalizer;
};

#endif /* DATASET_H */

