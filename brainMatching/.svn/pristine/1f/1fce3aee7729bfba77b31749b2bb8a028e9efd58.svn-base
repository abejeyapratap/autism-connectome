/* 
 * File:   parameterSet.h
 * Author: yusuf
 *
 * Generated on November 12, 2016, 3:46 PM
 */

#ifndef PARAMETERSET_H
#define PARAMETERSET_H

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "utility.h"

class ParameterSet
{
	public:
		enum ParameterType{BETA_PARAM=0,ASSIGNMENT_PARAMS,SEPARATION_PARAMS,ASSIGNMENT_AND_SEPARATION_PARAMS,END_OF_PARAMS};
		
		static const int MAX_VALUE=-1,BETA=0;
		
		ParameterSet();
		ParameterSet(int _numAssParams, int _numSepParams);
		ParameterSet(std::vector<float> _parameters, int _numAssParams, int _numSepParams, float _maxValue=-1);
		ParameterSet(float _beta, std::vector<float> _assignmentParameters, std::vector<float> _separationParameters, float _maxValue=-1);
		ParameterSet(float *data, int dataColumn, int _numAssParams, int _numSepParams);
		ParameterSet(std::string line, int _numAssParams, int _numSepParams);
		ParameterSet(const ParameterSet &other);
		ParameterSet& operator=(const ParameterSet& other);
		
		virtual void print(std::string mode="", std::ostream &out=std::cout, std::string endOfLine="\n");
		
		void setParameter(int paramId, float value);

		float getParameter(int orderId);
		enum ParameterType getParameterType(int orderId);
		inline int getNumberOfParameters(){return parameters.size()+1;}
		
		void addValueToParameter(int paramId, float value);
		void addParameters(ParameterSet &other);
		void multiplyParameters(float value);
		
		void generateParameters(std::vector<float> &paramStart, std::vector<float> &paramEnd, std::vector<float> &paramStep);
		void roundParameters(int digit=1);
		
		bool areParametersInRange(std::vector<float> &paramStart, std::vector<float> &paramEnd);
		bool areParametersSumGreaterThanOne(enum ParameterType parameterType);
		bool areParametersSumLessThanOne(enum ParameterType parameterType);
		bool areParametersSumEqualToOne(enum ParameterType parameterType);

	protected:
		int numAssignmentParams, numSeparationParams;
		std::vector<float> parameters;//beta,assignmentParamters,separationParameters
		
		float maxValue;//maxValue in a given column of data. 
					   //This field is used for choosing the best set of parameters after the evaluation phase. 
					   //This is not part of the input parameters for the objective function
		
		const float epsilon = 0.001;
		
		void generateParametersLooper(ParameterSet parameterSet,int parameterId, std::vector<float> &paramStart, std::vector<float> &paramEnd, std::vector<float> &paramStep);
};

#endif /* PARAMETERSET_H */

