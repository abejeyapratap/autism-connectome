/* 
* File:   parameterSet.cpp
* Author: yusuf
*
* Generated on February 27, 2017, 6:53 PM
*/

#include "parameterSet.h"

//<editor-fold defaultstate="collapsed" desc=" Constructors:">
ParameterSet::ParameterSet()
{
   numAssignmentParams=-1;
   numSeparationParams=-1;
}

ParameterSet::ParameterSet(int _numAssParams, int _numSepParams): numAssignmentParams(_numAssParams),numSeparationParams(_numSepParams)
{
   //insert null values for beta, assignment parameters, and separation parameters
   for(int i=0;i<numAssignmentParams+numSeparationParams+1;i++)
      parameters.push_back(0);

   maxValue = -1;
}

ParameterSet::ParameterSet(std::vector<float> _parameters, int _numAssParams, int _numSepParams, float _maxValue) 
  : parameters(_parameters),numAssignmentParams(_numAssParams),numSeparationParams(_numSepParams),maxValue(_maxValue)
{
}

ParameterSet::ParameterSet(float _beta, std::vector<float> _assignmentParameters, std::vector<float> _separationParameters, float _maxValue) 
  : maxValue(_maxValue)
{
   parameters.push_back(_beta);
   for(std::vector<float>::iterator iter = _assignmentParameters.begin(); iter!=_assignmentParameters.end();iter++)
      parameters.push_back(*iter);
   for(std::vector<float>::iterator iter = _separationParameters.begin(); iter!=_separationParameters.end();iter++)
      parameters.push_back(*iter);
   numAssignmentParams = _assignmentParameters.size();
   numSeparationParams = _separationParameters.size();
}

ParameterSet::ParameterSet(float* data, int dataColumn, int _numAssParams, int _numSepParams)
{
   parameters.push_back(data[0]);
   for(int i=0;i<_numAssParams;i++)
      parameters.push_back(data[i+1]);
   for(int i=0;i<_numSepParams;i++)
      parameters.push_back(data[i+1+_numAssParams]);
   
   maxValue = data[dataColumn];
   
   numAssignmentParams = _numAssParams;
   numSeparationParams = _numSepParams;
}

ParameterSet::ParameterSet(std::string line,int _numAssParams, int _numSepParams)
{
   std::istringstream is;
   float parameter;
   is.clear();
   is.str(line);
   while(is>>parameter)
      parameters.push_back(parameter);
   
   numAssignmentParams = _numAssParams;
   numSeparationParams = _numSepParams;
}


ParameterSet::ParameterSet(const ParameterSet &other)
{
   parameters = other.parameters;
   maxValue = other.maxValue;
   numAssignmentParams = other.numAssignmentParams;
   numSeparationParams = other.numSeparationParams;
}

ParameterSet& ParameterSet::operator=(const ParameterSet& other)
{
   parameters = other.parameters;
   maxValue = other.maxValue;
   numAssignmentParams = other.numAssignmentParams;
   numSeparationParams = other.numSeparationParams;
   return *this;
}
//</editor-fold>

//<editor-fold defaultstate="collapsed" desc=" Accessors/mutators/helpers: setParameter(), getParameter(), print()">
void ParameterSet::setParameter(int paramId, float value)
{
   if(paramId == MAX_VALUE)
      maxValue = value;
   else
      parameters[paramId] = value;
}

void ParameterSet::addValueToParameter(int paramId, float value)
{
   if(paramId == MAX_VALUE)
      maxValue += value;
   else
      parameters[paramId] += value;
}

void ParameterSet::addParameters(ParameterSet& other)
{
   int numParams = parameters.size();
   if(numParams!=other.parameters.size())
   {
      std::cerr<<"number of parameters mismatch. Cannot execute addParameters function. Exiting!!\n";
      exit(1);
   }
   for(int i=0;i<numParams;i++)
      parameters[i]+=other.parameters[i];
   maxValue += other.maxValue;
}

void ParameterSet::multiplyParameters(float value)
{
   for(std::vector<float>::iterator iter=parameters.begin();iter!=parameters.end();iter++)
      *iter *= value;
   maxValue *= value;
}

float ParameterSet::getParameter(int paramId)
{
   if(paramId == MAX_VALUE)
      return maxValue;
   else
      return parameters[paramId];
}

enum ParameterSet::ParameterType ParameterSet::getParameterType(int orderId)
{
   if(orderId==0)
      return BETA_PARAM;
   else if(orderId>0 && orderId<=numAssignmentParams)
      return ASSIGNMENT_PARAMS;
   else if(orderId>numAssignmentParams && orderId<=numAssignmentParams+numSeparationParams)
      return SEPARATION_PARAMS;
   else
   {
      std::cerr<<"getParameterType() function is called for incorrect parameter order:"<<orderId<<" !! Exiting...\n";
      return END_OF_PARAMS;
   }
}

void ParameterSet::print(std::string mode, std::ostream &out, std::string endOfLine)
{
   for(std::vector<float>::iterator iter = parameters.begin();iter!=parameters.end();iter++)
      out<<*iter<<"\t";

   if(mode.compare("all")==0)
     out<<maxValue<<endOfLine;
}
//</editor-fold>

//<editor-fold defaultstate="collapsed" desc=" Parameter modifiers: roundParameters(), generateParameters()">
void ParameterSet::roundParameters(int digit)
{
   float roundBy = std::pow(10.0,digit);//round by 10^digit

   for(std::vector<float>::iterator iter = parameters.begin();iter!=parameters.end();iter++)
      *iter = ((float)Utility::round((*iter)*roundBy))/roundBy;
}

void ParameterSet::generateParameters(std::vector<float> &paramStart, std::vector<float> &paramEnd, std::vector<float> &paramStep)
{
   generateParametersLooper(*this,0,paramStart,paramEnd,paramStep);
}

void ParameterSet::generateParametersLooper(ParameterSet parameterSet,int parameterId, std::vector<float> &paramStart, std::vector<float> &paramEnd, std::vector<float> &paramStep)
{
   ParameterType parameterType = parameterSet.getParameterType(parameterId);
   
   //if summation of the parameters of this type up to this orderId is larger than 1, we should stop the recursion
   if(parameterSet.areParametersSumGreaterThanOne(parameterType)==true)
      return;

   //if we recently completed one of the assignment or separation bundles, they should sum up to one, 
   //otherwise, we should stop the recursion
   if((parameterType==SEPARATION_PARAMS && parameterSet.areParametersSumEqualToOne(ASSIGNMENT_PARAMS)==false) || 
      (parameterType==END_OF_PARAMS && parameterSet.areParametersSumEqualToOne(SEPARATION_PARAMS)==false))
      return;
   
   bool isParameterFixed = false;
   if(Utility::absoluteValue(paramStart[parameterId]-paramEnd[parameterId])<=epsilon)
      isParameterFixed = true;
   
   parameterSet.setParameter(parameterId,paramStart[parameterId]);
   
   if(isParameterFixed)
   {
      generateParametersLooper(parameterSet,parameterId+1,paramStart,paramEnd,paramStep);
   }
   else
   {
      while(parameterSet.getParameter(parameterId)<=paramEnd[parameterId]+paramStep[parameterType]/10.0)
      {
         generateParametersLooper(parameterSet,parameterId+1,paramStart,paramEnd,paramStep);
         parameterSet.addValueToParameter(parameterId,paramStep[parameterType]);
      }
   }
}

//</editor-fold>

//<editor-fold defaultstate="collapsed" desc=" Parameter range functions: areParametersInRange(),parametersSumGreaterThanOne(),parametersSumToOne ">		
bool ParameterSet::areParametersInRange(std::vector<float> &paramStart, std::vector<float> &paramEnd)
{
   int numParams = 1 + numAssignmentParams + numSeparationParams;	

   for(int i=0;i<numParams;i++)
   {
      if(parameters[i]<paramStart[i] || parameters[i]>paramEnd[i])
         return false;
   }

   return true;
}

bool ParameterSet::areParametersSumGreaterThanOne(enum ParameterType parameterType)
{
   float assignmentParametersSum=1, separationParametersSum=1;

   if(parameterType==ASSIGNMENT_PARAMS || parameterType==ASSIGNMENT_AND_SEPARATION_PARAMS)
   {
      assignmentParametersSum=0;
      for(int i=0;i<numAssignmentParams;i++)
         assignmentParametersSum += parameters[i+1];
   }
   else if(parameterType==SEPARATION_PARAMS || parameterType==ASSIGNMENT_AND_SEPARATION_PARAMS)
   {
      separationParametersSum=0;
      for(int i=0;i<numSeparationParams;i++)
         separationParametersSum += parameters[i+1+numAssignmentParams];
   }

   if(assignmentParametersSum>1+epsilon || separationParametersSum>1+epsilon )
      return true;
   else
      return false;
}

bool ParameterSet::areParametersSumLessThanOne(enum ParameterType parameterType)
{
   float assignmentParametersSum=1, separationParametersSum=1;

   if(parameterType==ASSIGNMENT_PARAMS || parameterType==ASSIGNMENT_AND_SEPARATION_PARAMS)
   {
      assignmentParametersSum=0;
      for(int i=0;i<numAssignmentParams;i++)
         assignmentParametersSum += parameters[i+1];
   }
   if(parameterType==SEPARATION_PARAMS || parameterType==ASSIGNMENT_AND_SEPARATION_PARAMS)
   {
      separationParametersSum=0;
      for(int i=0;i<numSeparationParams;i++)
         separationParametersSum += parameters[i+1+numAssignmentParams];
   }

   if(assignmentParametersSum<1-epsilon || separationParametersSum<1-epsilon )
      return true;
   else
      return false;
}

//function testing if the coefficients of the parameters for the assignment and separation cost each sum to one
bool ParameterSet::areParametersSumEqualToOne(enum ParameterType parameterType)
{
   float assignmentParametersSum=1, separationParametersSum=1;

   if(parameterType==ASSIGNMENT_PARAMS || parameterType==ASSIGNMENT_AND_SEPARATION_PARAMS)
   {
      assignmentParametersSum=0;
      for(int i=0;i<numAssignmentParams;i++)
         assignmentParametersSum += parameters[i+1];
   }
   if(parameterType==SEPARATION_PARAMS || parameterType==ASSIGNMENT_AND_SEPARATION_PARAMS)
   {
      separationParametersSum=0;
      for(int i=0;i<numSeparationParams;i++)
         separationParametersSum += parameters[i+1+numAssignmentParams];
   }

   if(Utility::absoluteValue(assignmentParametersSum-1)<epsilon && Utility::absoluteValue(separationParametersSum-1)<epsilon)
      return true;
   else
      return false;
}
//</editor-fold>

