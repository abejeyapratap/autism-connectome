/* 
 * File:   test.h
 * Author: yusuf
 *
 * Generated on May 1, 2018, 11:58 AM
 */

#ifndef TEST_H
#define TEST_H


#include <iostream> //cin,cout
#include <cstdlib> //atoi,stof, etc.
#include <string> //std::stof Note: stof requires C++11 as the compiler standard
#include <vector>
#include "graph.h"
#include "utility.h"

class Test
{
	public:
		Test();
		~Test(){};
		
		////////sandbox: functions for testing stuff
		static void shuffleGraph(std::string inFile, std::string outFile, int seedSupplement, int iterationBound, bool isFunctional=false);
		static void calculateCovarianceOfAGraph(std::string inFile);
		static void testRandom(float min, float max, int randCount, int binCount);
		static void testCommunicability(int numNodes,std::string inputMatrix, std::string referenceMatrix, std::string outputMatrix);
		static void testCode();
};

#endif /* TEST_H */

