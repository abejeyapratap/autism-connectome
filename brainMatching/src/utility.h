#ifndef _UTILITY_H
#define _UTILITY_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cmath>
#include <time.h>
#include <stdio.h>
#include <vector>
#include <map>
#include <sys/time.h> 
#include <dirent.h>//for listing contents of a folder
#include <cctype>//for tolower function
#include <algorithm>//for using std::sort()
#include <cstdarg> //for variadic functions, va_arg() etc
#include "geometry.h"//for the constants such as PI

//using namespace std;
//using std::string;
//using std::vector;
//using std::pair;

class Utility
{
public:
	Utility(){};
	~Utility(){};
	
	static constexpr float EPSILON = 1e-16;
	static constexpr float INF=std::numeric_limits<float>::max();//10000000;
	static constexpr float N_INF=-INF;

	// <editor-fold defaultstate="collapsed" desc=" Print and string related functions: stringify/numerify(), toLower/toUpper(), tokenizeSentence(), splitString(), printTab(), printStringstreamToFile()">
	template <typename Type>
	static inline std::string stringify(Type x)
	{
		std::ostringstream o;
		o << x;
		return o.str();
	}

	template <typename Type>
	static inline Type numerify ( const std::string &text )
	{
		std::stringstream ss(text);
		Type result;
		return ss >> result ? result : 0;
	}

	static inline void toLower( std::string &text)
	{
		for(int i=0;i<text.length();i++)
			text[i]=std::tolower(text[i]);
	}

	static inline void toUpper( std::string &text)
	{
		for(int i=0;i<text.length();i++)
			text[i]=std::toupper(text[i]);
	}

   //splits a given @text into a vector of strings named @tokens by separating from white spaces and new lines
	static inline void tokenizeSentence(std::string text, std::vector<std::string> &tokens)
	{
		std::istringstream iss(text);
		while(iss)
		{
			std::string token;
			iss >> token;
			tokens.push_back(token);
		}
		//if the last character of the sentence was a new line, we don't want to take it as a token
		if(tokens[tokens.size()-1].length()==0)
			tokens.erase(tokens.end()-1);
	}
	
   //give a @text and a list of delimiters stored in a single string @delims, splits the text and saves them into a string vector named @tokens
	static inline void splitString(const std::string& text, const std::string& delims, std::vector<std::string> &tokens)
	{
		std::size_t start = text.find_first_not_of(delims), end = 0;

		while((end = text.find_first_of(delims, start)) != std::string::npos)
		{
			tokens.push_back(text.substr(start, end - start));
			start = text.find_first_not_of(delims, end);
		}
		if(start != std::string::npos)
			tokens.push_back(text.substr(start));
	}
   
   //give a @text and a list of delimiters stored in a single string @delims, splits the text and 
   //returns the token in @orderNumber
   static inline std::string splitString(const std::string& text, const std::string& delims, int orderNumber)
	{
		std::size_t start = text.find_first_not_of(delims), end = 0;
		std::vector<std::string> tokens;

		while((end = text.find_first_of(delims, start)) != std::string::npos)
		{
			tokens.push_back(text.substr(start, end - start));
			start = text.find_first_not_of(delims, end);
		}
		if(start != std::string::npos)
			tokens.push_back(text.substr(start));
      
		if(orderNumber>=0)
		   return tokens[orderNumber];
		else
		{
		   int numTokens=tokens.size();
		   return tokens[numTokens+orderNumber];
		}  
	}

	static inline void printTab(int count)
	{
		for(int i=0;i<count;i++)
		std::cout<<"\t|";
	}

	static inline void loadFileIntoVectorOfIstringstream(std::string filename, std::vector<std::istringstream> &vss)
	{
		std::ifstream file; 
		std::string line;  

		file.open(filename.c_str());
		if(file.fail())
		{
		  std::cout << "File not found in loadFileIntoVectorOfIstringstream():" << filename << "\n";
		  return;
		}
		
		getline(file, line);	//ignore comment line
		while(getline(file, line))
		{
			vss.push_back(std::istringstream(line));
		}
	}
	
	static inline void printOstreamToFile(std::string filename, std::ostream &os)
	{
		std::ofstream file;
		file.open(filename.c_str());
		
		file<<os.rdbuf();//read content of the string stream into the file stream
		file.close();
	}
	
	static inline void appendOstreamToFile(std::string filename, std::ostream &os)
	{
		std::ofstream file;
		file.open(filename.c_str(),std::ios_base::app);
		
		file<<os.rdbuf();//read content of the string stream into the file stream
		file.close();
	}
	
	static inline void printStringToFile(std::string filename, std::string str)
	{
		std::ofstream file;
		file.open(filename.c_str());
		
		file<<str;
		file.close();
	}
	
	static inline void appendStringToFile(std::string filename, std::string str)
	{
		std::ofstream file;
		file.open(filename.c_str(),std::ios_base::app);
		
		file<<str;
		file.close();
	}
	
	//</editor-fold>

	//<editor-fold defaultstate="collapsed" desc=" c-style helper functions for reading command line parameters in main()">
	static char* getCmdOption(char ** begin, char ** end, const std::string & option)
	{
	   char ** itr = std::find(begin, end, option);
	   if (itr != end && ++itr != end)
	   {
		  return *itr;
	   }
	   return const_cast<char*>("");
	}
	
	static char* getCmdOptionWithOrder(char ** begin, char ** end, const std::string & option, int order)
	{
	   char ** itr = std::find(begin, end, option);
	   if (itr != end && itr+order != end)
	   {
		  return *(itr+order);
	   }
	   return const_cast<char*>("");
	}

	static bool cmdOptionExists(char** begin, char** end, const std::string& option)
	{
	   return std::find(begin, end, option) != end;
	}
	
	static std::string getCommandLine(int argc,char** argv)
	{
		std::string commandLine="";
		for(int i=0;i<argc;i++)
		   commandLine += std::string(argv[i]) + " ";
		commandLine += "\n";
		return commandLine;
	}

	//</editor-fold>
	
	// <editor-fold defaultstate="collapsed" desc=" Calculate Distance functions: calculateDistance()x3, calculateDotProduct()x2, calculateCosineDistance(), calculateMagnitude()">     
	//calculates Euclidean (i.e.,L2) distance between two vectors
	template <typename T>
	static T calculateL1Distance(const std::vector<T> &v1,const std::vector<T> &v2)
	{
		if(v1.size()!=v2.size())
		{
		   std::cerr<<"Utility::calculateDistance() cannot calculate distance btw vectors of different size\n";
		   return -1;
		}

		T distance=0;
		int size = v1.size();
		for(int i=0;i<size;i++)
			distance += std::abs(v1[i]-v2[i]);

		return distance;
	}
	
	template <typename T>
	static inline T calculateL2Distance(T x1, T y1, T x2, T y2)
	{
		T distance = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2);
		return sqrt(distance);
	}

	//calculates Euclidean (i.e.,L2) distance between two vectors
	template <typename T>
	static T calculateL2Distance(const std::vector<T> &v1,const std::vector<T> &v2)
	{
		if(v1.size()!=v2.size())
		{
		   std::cerr<<"Utility::calculateDistance() cannot calculate distance btw vectors of different size\n";
		   return -1;
		}

		T distance=0;
		int size = v1.size();
		for(int i=0;i<size;i++)
			distance += (v1[i]-v2[i])*(v1[i]-v2[i]);

		return sqrt(distance);
	}

	template <typename T>
	static T calculateL2Distance(const T* Mi, const T* Mj, int size)
	{
		T distance=0;
		for(int i=0;i<size;i++)
		   distance += (Mi[i]-Mj[i])*(Mi[i]-Mj[i]);

		return sqrt(distance);
	}
	
	template <typename T>
	static T calculateL2DistanceByIgnoringPairwiseRelations(const T* Mi, const T* Mj, int size, int i,int j)
	{
		T distance=0;
		for(int k=0;k<size;k++)
		{
			if(k!=i && k!=j)
				distance += (Mi[k]-Mj[k])*(Mi[k]-Mj[k]);
		}
		if(i==j)
			distance /= (float)(size-1);
		else
			distance /= (float)(size-2);
		
		return sqrt(distance);
	}
	
	static int calculateLevenshteinDistance(std::string s, std::string t)
	{
		// degenerate cases
		if (s == t) return 0;
		if (s.size() == 0) return t.size();
		if (t.size() == 0) return s.size();

		// create two work vectors of integer distances
		int v0length = t.size() + 1;
		int* v0 = new int[v0length];
		int* v1 = new int[v0length];

		// initialize v0 (the previous row of distances)
		// this row is A[0][i]: edit distance for an empty s
		// the distance is just the number of characters to delete from t
		for (int i = 0; i < v0length; i++)
			v0[i] = i;

		for (int i = 0; i < s.size(); i++)
		{
			// calculate v1 (current row distances) from the previous row v0

			// first element of v1 is A[i+1][0]
			//   edit distance is delete (i+1) chars from s to match empty t
			v1[0] = i + 1;

			// use formula to fill in the rest of the row
			for (int j = 0; j < t.size(); j++)
			{
				int cost = (s[i] == t[j]) ? 0 : 1;
				v1[j + 1] = minimum<int>(3,v1[j] + 1, v0[j + 1] + 1, v0[j] + cost);
			}

			// copy v1 (current row) to v0 (previous row) for next iteration
			for (int j = 0; j < v0length; j++)
				v0[j] = v1[j];
		}

		delete[] v0;
		delete[] v1;

		return v1[t.size()];
	}

	template <typename T>
	static int calculateLevenshteinDistance(std::vector<T> s, std::vector<T> t)
	{
		// degenerate cases
		if (s == t) return 0;
		if (s.size() == 0) return t.size();
		if (t.size() == 0) return s.size();

		// create two work vectors of integer distances
		int v0length = t.size() + 1;
		int* v0 = new int[v0length];
		int* v1 = new int[v0length];

		// initialize v0 (the previous row of distances)
		// this row is A[0][i]: edit distance for an empty s
		// the distance is just the number of characters to delete from t
		for (int i = 0; i < v0length; i++)
			v0[i] = i;

		for (int i = 0; i < s.size(); i++)
		{
			// calculate v1 (current row distances) from the previous row v0

			// first element of v1 is A[i+1][0]
			//   edit distance is delete (i+1) chars from s to match empty t
			v1[0] = i + 1;

			// use formula to fill in the rest of the row
			for (int j = 0; j < t.size(); j++)
			{
				int cost = (s[i] == t[j]) ? 0 : 1;
				v1[j + 1] = minimum<int>(3,v1[j] + 1, v0[j + 1] + 1, v0[j] + cost);
			}

			// copy v1 (current row) to v0 (previous row) for next iteration
			for (int j = 0; j < v0length; j++)
				v0[j] = v1[j];
		}

		delete[] v0;
		delete[] v1;

		return v1[t.size()];
	}

	template <typename T>
	static T calculateDotProduct(const std::vector<T> &v1,const std::vector<T> &v2)
	{
		if(v1.size()!=v2.size())
		{
		   std::cerr<<"Utility::calculateDotProduct() cannot calculate distance btw vectors of different size\n";
		   return -1;
		}
		
		T distance=0;
		int size = v1.size();
		for(int i=0;i<size;i++)
		   distance += v1[i]*v2[i];

		return distance;
	}
	
	template <typename T>
	static T calculateDotProduct(T* vec1, T* vec2, int size)
	{
		T result = 0;

		for(int i=0;i<size;i++)
			result += vec1[i]*vec2[i];
		return result;
	}


	template <typename T>
	static T calculateMagnitude(const std::vector<T> &v)
	{
		T distance = 0;
		for(typename std::vector<T>::iterator iter=v.begin();iter!=v.end();iter++)
			distance += (*iter)*(*iter);

		return (T)sqrt(distance);
	}

	template <typename T>
	static T calculateCosineDistance(const std::vector<T> &v1,const std::vector<T> &v2)
	{
		T dotProduct = calculateDotProduct<T>(v1,v2);
		T magnitudeOfV1 = calculateMagnitude<T>(v1);
		T magnitudeOfV2 = calculateMagnitude<T>(v2);

		return dotProduct/(magnitudeOfV1*magnitudeOfV2);
	}

	//</editor-fold>

	// <editor-fold defaultstate="collapsed" desc=" Math functions: sphericalToCartesian(), calculateDotProduct(), gcd(),lcm(), absoluteValue(), round(), roundUp(),log2(), log(), log10(), log_n(), square()">
	static void sphericalToCartesian(float phi, float theta, float rho, float &x, float &y, float &z)
	{
		x = rho*std::sin(phi)*std::cos(theta);
		y = rho*std::sin(phi)*std::sin(theta);
		z = rho*std::cos(phi);
	}

	//calculates greatest common divisor of two integers
	static inline int gcd(int a, int b)
	{
		for (;;)
		{
			if (a == 0) 
				return b;
			b %= a;
			if (b == 0) 
				return a;
			a %= b;
		}
	}

	//calculates least common multiple of two integers
	static inline int lcm(int a, int b)
	{
	  int temp = gcd(a, b);

	  return temp ? (a / temp * b) : 0;
	}

	template <typename Type>
	static inline Type absoluteValue(Type num)
	{
		return (num < 0 ? -num : num);
	}
	
	template <typename Type>
 	static inline Type sign(Type num)
 	{
 		if(num<0)
 			return -1;
 		else if(num>0)
 			return 1;
 		else 
 			return 0;
 	}
	
	template <typename T>
	static T average(int count, ...)
	{
		va_list ap;
		T tot = 0;
		va_start(ap, count); //Requires the last fixed parameter (to get the address)
		for(int j=0; j<count; j++)
			tot+=va_arg(ap, T); //Requires the type to cast to. Increments ap to the next argument.
		va_end(ap);
		return tot/count;
	}

	template <typename T>
	static T minimum(int count, ...)
	{
		va_list ap;

		va_start(ap, count); //Requires the last fixed parameter (to get the address)
		T min = va_arg(ap, T);
		for(int j=1; j<count; j++)
		{
			T temp = va_arg(ap, T);
			if(temp<min)
				min = temp;
		}
		va_end(ap);
		return min;
	}
	
	template <typename T>
	static T maximum(int count, ...)
	{
		va_list ap;

		va_start(ap, count); //Requires the last fixed parameter (to get the address)
		T max = va_arg(ap, T);
		for(int j=1; j<count; j++)
		{
			T temp = va_arg(ap, T);
			if(temp>max)
				max = temp;
		}
		va_end(ap);
		return max;
	}

	template <typename Type>
	static inline int round(Type num)
	{
		int temp;
		temp = (int)(num*10);
		temp %= 10;
		if(temp>=5)
			return (int)std::floor(num)+1;
		else
			return (int)std::floor(num);
	}

	template <typename Type>
	static inline int roundUp(Type num)
	{
		int temp;
		temp = (int)(num*10);
		temp %= 10;
		if(temp>0)
			return (int)std::floor(num)+1;
		else
			return (int)std::floor(num);
	}

	template <typename Type>
	static inline double log(Type num,float base)
	{
		return (double)std::log(num)/(double)std::log(base);
	}

	template <typename Type>
	static inline double log2(Type num)
	{
		return (double)std::log(num)/(double)std::log(2);
	}
	
	template <typename Type>
	static inline double log10(Type num)
	{
		return (double)std::log10(num);
	}
	
	template <typename Type>
	static inline double log_n(Type num)
	{
		return (double)std::log(num);
	}

	template <typename Type>
	static inline float square(Type num)
	{
		return num*num;
	}
	
	static inline int factorial(int num)
	{
		int val=1;
		for(int i=num;i>0;i--)
			val *= i;
		return val;
	}

	//</editor-fold>

	// <editor-fold defaultstate="collapsed" desc=" Statistics functions: calculateMean(), calculateStandardDeviation(), calculateCovariance(), calculateCorrelation(), fisherZTransform() ">
	template <typename Type>
	static inline float calculateMean(std::vector<Type> &vec)
	{
		float mean = 0;
		
		for(typename std::vector<Type>::iterator iter=vec.begin();iter!=vec.end();iter++)
			mean += *iter;
		mean /= (float)vec.size();
		
		return mean;
	}
	
	template <typename Type>
	static float calculateStandardDeviation(std::vector<Type> &vec, Type mean)
	{
		float stdDeviation = 0;
		
		for(typename std::vector<Type>::iterator iter=vec.begin();iter!=vec.end();iter++)
			stdDeviation += (*iter - mean)*(*iter - mean);
		stdDeviation /= (float)vec.size();
		stdDeviation = sqrt(stdDeviation);
		
		return stdDeviation;
	}
	
	template <typename Type>
	static float calculateStandardDeviation(std::vector<Type> &vec)
	{
		float mean=0, stdDeviation = 0;
		
		for(typename std::vector<Type>::iterator iter=vec.begin();iter!=vec.end();iter++)
			mean += *iter;
		mean /= (float)vec.size();
		
		for(typename std::vector<Type>::iterator iter=vec.begin();iter!=vec.end();iter++)
			stdDeviation += (*iter - mean)*(*iter - mean);
		stdDeviation /= (float)vec.size();
		stdDeviation = sqrt(stdDeviation);
		
		return stdDeviation;
	}
	
	//cov(x,y) = (sum_{i=1..n} (x_i-x_m)(y_y-y_m)) / n-1
	template <typename Type>
	static float calculateCovariance(std::vector<Type> &vec1, std::vector<Type> &vec2)
	{
		float mean1=0, mean2=0;
		float covariation = 0;
		
		int size = vec1.size();
		if(vec2.size()!=size)
		{
			std::cerr<<"Vector sizes mismatch for calculateCovariation()!! Exiting...\n";
			exit(1);
		}
		
		//calculate mean1
		for(typename std::vector<Type>::iterator iter=vec1.begin();iter!=vec1.end();iter++)
			mean1 += *iter;
		mean1 /= (float)size;
		
		//calculate mean2
		for(typename std::vector<Type>::iterator iter=vec2.begin();iter!=vec2.end();iter++)
			mean2 += *iter;
		mean2 /= (float)size;
		
		//calculate covariation
		for(int i=0;i<size;i++)
			covariation += (vec1[i]-mean1)*(vec2[i]-mean2);
		covariation /= (float)size;
		
		return covariation;
	}
	
	template <typename Type>
	static inline float calculateZScore(Type value, std::vector<Type> &distribution)
	{
		float mean = calculateMean(distribution);
		float std = calculateStandardDeviation(distribution,mean);
		float zScore=0;
		if(std!=0)
			zScore=((float)value-mean)/std;
		
		return zScore;
	}
	
	template <typename Type>
	static inline float calculateMeanSquaredError(std::vector<Type> &prediction, std::vector<Type> &truth)
	{
		float error = 0;
		int vectorSize = prediction.size();
		for(int i=0;i<vectorSize;i++)
		   error += (prediction[i]-truth[i])*(prediction[i]-truth[i]);
		error /= (float)vectorSize;
		return error;
	}
	
    template <typename Type>
	static inline float calculateMeanAbsoluteError(std::vector<Type> &prediction, std::vector<Type> &truth)
	{
		float error = 0;
		int vectorSize = prediction.size();
		for(int i=0;i<vectorSize;i++)
		   error += absoluteValue<float>(prediction[i]-truth[i]);
		error /= (float)vectorSize;
		return error;
	}
	
	//this function calculates the Pearson's correlation coefficient "r" for two given vectors
	//corr(x,y) = r(x,y) = cov(x,y) / s_x*s_y
	//s_x and s_y are the standard deviation of both random variables
	//s_x = sqrt( (sum_{i=1..n} (x_i-x_m)^2 )/n-1 )
	template <typename Type>
	static float calculateCorrelation(std::vector<Type> &vec1, std::vector<Type> &vec2)
	{
		float mean1=0, mean2=0;
		float covariation = 0;
		float stdDeviation1=0,stdDeviation2=0;
		float correlation;
		
		int size = vec1.size();
		if(vec2.size()!=size)
		{
			std::cerr<<"Vector sizes mismatch for calculateCorrelation()!! Exiting...\n";
			exit(1);
		}
		
		//calculate covariance but don't divide into number of samples
		for(typename std::vector<Type>::iterator iter=vec1.begin();iter!=vec1.end();iter++)
			mean1 += *iter;
		mean1 /= (float)size;
		
		for(typename std::vector<Type>::iterator iter=vec2.begin();iter!=vec2.end();iter++)
			mean2 += *iter;
		mean2 /= (float)size;
		
		for(int i=0;i<size;i++)
			covariation += (vec1[i]-mean1)*(vec2[i]-mean2);
		
		
		//calculate standard deviation for both random variables
		for(typename std::vector<Type>::iterator iter=vec1.begin();iter!=vec1.end();iter++)
			stdDeviation1 += (*iter - mean1)*(*iter - mean1);
		for(typename std::vector<Type>::iterator iter=vec2.begin();iter!=vec2.end();iter++)
			stdDeviation2 += (*iter - mean2)*(*iter - mean2);
		
		//calculate correlation
		correlation = covariation / sqrt(stdDeviation1*stdDeviation2);
		
		return correlation;
	}
	
   //given a correlation value r, calculates its Fisher's z-transform
   //formula is: z = artanh(r) =1/2*ln((1+r)/(1-r)) where r is the correlation between two vectors of data
	template <typename Type>
	static inline float fisherZTransform(Type r)
	{
		if(r>1)
		{
			std::cerr<<"attempting to appy Fisher's z-transform to an r value greater than 1:"<<r<<"\nExiting..."<<std::endl;
			exit(1);
		}
		float nominator = 1.0+r;
		float denominator = 1.0-r;
		denominator = (denominator < EPSILON ? denominator+EPSILON : denominator);
		
		return 0.5*std::log(nominator/denominator);
	}	
	// </editor-fold>
	
	// <editor-fold defaultstate="collapsed" desc=" Random number generators: geometricDistribution(), exponentialDistribution(), gaussianDistribution(), generateRandomNumber(),randDiscrete(), PRNG()">
	static double geometricDistribution(double p)
	{
		int count=0;
		while(true)
		{
		   count++;
		   if(PRNG() < p)
			  return (double)count;
		}
	}

	static double exponentialDistribution(double r, double n)
	{
		double z;//random radius z
		double p_z;//probability of z
		double rlnn = r*std::log(n);
		while(true)
		{
			//x is a random value in interval [0,rlnn)
			z = (double)(std::rand()%((int)(rlnn*1000)))/1000.0;
			//y is a random value in interval [1/(n-1)*r , n/(n-1)*r]
			//first pick a random number in interval [0 , 1/r]
			p_z = (double)(std::rand()%((int)(1000000.0/r)))/1000000.0;
			//then add 1/(n-1)*r to the generated intermediate y value to get final random y value
			p_z += 1.0/((n-1)*r);
			//if generated value fits our probability distribution, return it, otherwise try once more
			if(p_z < ((n*std::exp(-z/r))/((n-1)*r)) )
				return z;
			else 
				continue;       
		}
	}

	static double gaussianDistribution(double sigma, double r)
	{
		double sqrt2Pi = sqrt(2*Geometry::PI);
		double denominator = 1.0/(sigma*sqrt2Pi);
		double upperLimit = denominator;

		double randR, p_randR;

		while(true)
		{
			//randR is a random value in interval [0,r)
			randR = (double)(std::rand()%((int)(r*1000)))/1000.0;
			//p_randR is a random value for probability of randR, between [0,upperLimit] interval
			p_randR = (double)(std::rand()%((int)(1000000.0*upperLimit)))/1000000.0;;
			//if generated value fits our probability distribution, return it, otherwise try once more
			if(p_randR <= std::exp(-(((double)randR/sigma)*((double)randR/sigma))/2) * denominator )
				return randR;
			else 
				continue;       
		}
	}

	// a random number generator using rand() function of standard library
	static inline float generateRandomNumber(int seed, int interval)
	{
		float temp;
		int halfInterval = interval/2;

		std::srand ( seed );
		while(1)
		{
			temp = std::rand() % interval;
			temp -= halfInterval;
			if(temp!=0)
				break;
			else 
				continue;
		}
		return temp;
	}
	
	// a random number generator using rand() function of standard library
	static inline float generateRandomNumber(int seed, float rangeMin, float rangeMax)
	{
		float randomNumber;
		int range = 13849;//some arbitrary large range to randomly pick a number from

		std::srand ( (int)clock() + seed );
		while(1)
		{
			randomNumber = std::rand() % range;

			if(randomNumber!=0)
				break;
			else 
				continue;
		}
		
		randomNumber = rangeMin + ((float)randomNumber * (rangeMax-rangeMin))/(float)range;
		
		return randomNumber;
	}

	//returns a random number among the given set of values
	static float randDiscrete(int seed, std::vector<float> &values)
	{
		std::srand ( seed );
		int temp = std::rand() % values.size();
		return values[temp];
	}

	//another more accurate (more uniform) pseudo random number generator
	static inline double PRNG()
	{
		const int C = 25173;
		const int D = 13849;
		const int M = 32768;

		static int seedPRNG = 5323;

		seedPRNG = (C*seedPRNG+D)%M;

		return (double)seedPRNG/(double)M;  
	}

	static inline Vector generateRandomVector(int seed, int interval)
	{
		float temp[3];
		Vector returnVal;
		int halfInterval = interval/2;

		std::srand ( seed );
		for(int i=0;i<3;i++)
		{
			while(1)
			{
				temp[i] = std::rand() % interval;
				temp[i] -= halfInterval;
				if(temp[i]!=0)
					break;
				else 
					continue;
			}	

		}
		returnVal.set(temp[0]*0.001,temp[1]*0.001,temp[2]*0.001);
		return returnVal;
	}
	//</editor-fold>

	// <editor-fold defaultstate="collapsed" desc=" File manipulators: sumNumbers(), mergeFiles(), loadFileNames()--2 version, selectivelyCopyFiles(), writeTGA() ">
	
	//this function reads a given file and returns the content of the indicated line as a string
	static std::string readLineFromFile(std::string filename, int lineNumber)
	{
		std::ifstream file; 
		std::istringstream is; 
		std::string line;

		////load names of result files
		file.open(filename.c_str());
		if(file.fail())
		{
			std::cerr << "File not found readLineFromFile(): "<<filename<<" !!!!\n";
			exit(1);
		}	

		int count =0;
		bool found=false;
		while(getline(file, line))
		{
			if(count==lineNumber)
			{
				found=true;
				break;
			}
			else
				count++;
		}
		file.close();
		if(found==false)
			return "";
		else
			return line;
	}
	
	//reads a file that has a number at its each line, and returns the sum of those numbers
	static float sumNumbers(std::string filename)
	{
		std::ifstream inFile; 
		std::istringstream is;
		std::string line;
		float num;
		float sum=0;

		inFile.open(filename.c_str());
		if(inFile.fail())
		{
			std::cerr << "File not found sumNumbers():"<<filename<<std::endl;
			exit(1);
		}

		while(getline(inFile, line))
		{
			is.clear();
			is.str(line);
			is >> num;

			sum += num;
		}
		inFile.close();

		return sum;
	}

	//given a file, counts the number of lines in it
	static int countRows(std::string filename)
	{
		std::ifstream inFile; 
		std::string line;
		int count=0;

		inFile.open(filename.c_str());
		if(inFile.fail())
		{
			std::cerr << "File not found countRows():"<<filename<<std::endl;
			exit(1);
		}

		while(getline(inFile, line))
		{
			count++;
		}
		inFile.close();

		return count;
	}
	
	static int countRowsAfterSkipline(std::string filename, int skipline)
	{
		std::ifstream inFile; 
		std::string line;
		int count=0;

		inFile.open(filename.c_str());
		if(inFile.fail())
		{
			std::cerr << "File not found countRowsAfterSkipline():"<<filename<<std::endl;
			exit(1);
		}

		for(int i=0;i<skipline;i++)
			getline(inFile, line);//ignore skiplines
		
		while(getline(inFile, line))
			count++;

		inFile.close();

		return count;
	}
	
	static int countColumns(std::string filename)
	{
		std::ifstream inFile; 
		std::istringstream is;
		std::string line;
		int count=0;

		inFile.open(filename.c_str());
		if(inFile.fail())
		{
			std::cerr << "File not found in countColumns():"<<filename<<std::endl;
			exit(1);
		}

		getline(inFile, line);
		is.clear();
		is.str(line);
		
		std::string temp;

		while(is >> temp)
			count++;

		inFile.close();

		return count;
	}
	
	static int countColumnsAfterSkipline(std::string filename, int skipline)
	{
		std::ifstream inFile; 
		std::istringstream is;
		std::string line;
		int count=0;

		inFile.open(filename.c_str());
		if(inFile.fail())
		{
			std::cerr << "File not found in countColumnsAfterSkipline():"<<filename<<std::endl;
			exit(1);
		}

		for(int i=0;i<skipline;i++)
			getline(inFile, line);//ignore skiplines
		
		getline(inFile, line);
		is.clear();
		is.str(line);
		
		std::string temp;

		while(is >> temp)
			count++;

		inFile.close();

		return count;
	}
	
	//given a folderPath that includes the output files, and list of folder content in listFilename,
	//merges the files as a single file that is to be named outputFilename
	static void mergeFiles(std::string folderPath, std::string outputFilename, std::string listFilename="")
	{
		std::ifstream inFile; 
		std::istringstream is; 
		std::vector<std::string> resultFiles;
		std::string line, fileName, source; 
		std::cout<<folderPath<<"  "<<listFilename<<"  "<<outputFilename<<"\n";
		
		////load names of result files
		if(listFilename.compare("")==0)
			loadFilePathFromList(folderPath,listFilename,resultFiles);
		else
			loadFilePathFromFolder(folderPath,resultFiles);

		///now merge files
		std::ofstream outFile;
		outFile.open(outputFilename.c_str());

		for(int i=0;i<resultFiles.size();i++)
		{
			inFile.open(resultFiles[i].c_str());
			while(getline(inFile, line))
			{
				outFile << line << std::endl;
			}
			inFile.close();
		}

		outFile.close(); 
	}

	//given a list containing of filenames
	//this function returns the list of the files in a string vector
	static void loadFileNamesFromList(std::string fileList, std::vector<std::string> &fileNames)
	{
		std::ifstream file; 
		std::istringstream is; 
		std::string line, fileName, source;

		////load names of result files
		file.open(fileList.c_str());
		if(file.fail())
		{
			std::cerr << "File not found in loadFileNamesFromList():"<<fileList<<std::endl;
			exit(1);
		}	

		while(getline(file, line))
		{
			is.clear();
			is.str(line);
			is >> fileName;

			fileNames.push_back(fileName);
		}
		file.close();
	}

	//given a folderPath and a list of files that are included under this folder
	//this function returns the list of the files in a string vector by appending the folder path to each filename
	static void loadFilePathFromList(std::string folderPath, std::string fileList, std::vector<std::string> &fileNames)
   	{
		std::ifstream file; 
		std::istringstream is; 
		std::string line, fileName, source;

		////load names of result files
		file.open(fileList.c_str());
		if(file.fail())
		{
		   std::cerr << "File not found in loadFilePathFromList():"<<fileList<<std::endl;
			exit(1);
		}	

		while(getline(file, line))
		{
		   is.clear();
		   is.str(line);
		   is >> fileName;

		   source = folderPath;
		   source.append(fileName);
		   fileNames.push_back(source);
		}
		file.close();
	}

	//given a folderPath this function returns the list of the filenames in a 
	//string vector without appending the folder path to each filename
	static void loadFilenamesFromFolder(std::string folderPath, std::vector<std::string> &filenames)
	{
		unsigned char isFile =0x8, isFolder =0x4;
		DIR *Dir;
		struct dirent *DirEntry;
		std::string filePath;

      if( (Dir = opendir(folderPath.c_str())) == NULL)
      {
         std::cerr << "Error: Path not found in loadFilenamesFromFolder()... Path: "<<folderPath<<std::endl;
         exit(1);
      }
      
		while(DirEntry=readdir(Dir))
		{
			if ( DirEntry->d_type == isFile)
			{
				filePath = DirEntry->d_name;
				filenames.push_back(filePath);
			}
		}
		std::sort(filenames.begin(), filenames.end());
		closedir(Dir);
	}
	
	//given a folderPath this function returns the list of the files in a 
	//string vector by appending the folder path to each filename
	static void loadFilePathFromFolder(std::string folderPath, std::vector<std::string> &filenames)
	{
		unsigned char isFile =0x8, isFolder =0x4;
		DIR *Dir;
		struct dirent *DirEntry;
		std::string filePath;
      
      if( (Dir = opendir(folderPath.c_str())) == NULL)
      {
         std::cerr << "Error: Path not found in loadFilePathFromFolder()... Path: "<<folderPath<<std::endl;
         exit(1);
      }

		while(DirEntry=readdir(Dir))
		{
			if ( DirEntry->d_type == isFile)
			{
				filePath = folderPath;
				filePath.append(DirEntry->d_name);
				filenames.push_back(filePath);
			}
		}
		std::sort(filenames.begin(), filenames.end());
		closedir(Dir);
	}
	
	//given a file containing the list of filenames and 
	//a list of "object names" that are the substrings of the filenames 
	//this function returns the list of the filenames corresponding to the
	//object list in a string vector 
	static void loadFilePathFromListSelectively(std::string objectListFilename, std::string filenameListFilename, std::vector<std::string> &selectedFilenames)
	{
		std::vector<std::string> objectNames,allFilenames;
		loadFileNamesFromList(objectListFilename,objectNames);
		loadFileNamesFromList(filenameListFilename,allFilenames);

		int numAllFilenames = allFilenames.size();
		for(int i=0;i<numAllFilenames;i++)
		{
			int numRemainingObjects=objectNames.size();
			for(int j=0;j<numRemainingObjects;j++)
			{
				if(allFilenames[i].find(objectNames[j]) != std::string::npos)
				{
					selectedFilenames.push_back(allFilenames[i]);
					objectNames.erase(objectNames.begin()+j);
					break;
				}
			}
		}

		std::sort(selectedFilenames.begin(), selectedFilenames.end());
	}

	//given a list of "object names" that are the substrings of the filenames
	//and a folder path that includes the files that we want to obtain the paths for,
	//this function returns the list of the filenames corresponding to the
	//object list in a string vector 
	static void loadFilePathFromFolderSelectively(std::string objectListFilename, std::string folderPath, std::vector<std::string> &selectedFilenames)
	{
		std::vector<std::string> objectNames;
		loadFileNamesFromList(objectListFilename,objectNames);

		unsigned char isFile =0x8, isFolder =0x4;
		DIR *Dir;
		struct dirent *DirEntry;
      std::string filePath;
      
      if( (Dir = opendir(folderPath.c_str())) == NULL)
      {
         std::cerr << "Error: Path not found in loadFilePathFromFolderSelectively()... Path: "<<folderPath<<std::endl;
         exit(1);
      }
		
		while(DirEntry=readdir(Dir))
		{
			if ( DirEntry->d_type == isFile)
			{
				int objectCount = objectNames.size();
				for(int i=0;i<objectCount;i++)
				{
					if((std::string(DirEntry->d_name)).find(objectNames[i]) != std::string::npos)
					{
						filePath = folderPath;
						filePath.append(DirEntry->d_name);
						selectedFilenames.push_back(filePath);
						objectNames.erase(objectNames.begin()+i);
						break;
					}
				}
			}
		}
		if(objectNames.size()!=0)
			std::cerr<<"Warning: "<<objectNames.size()<<" file paths couldn't be loaded  as the files do not exist!!!\n";
		
		std::sort(selectedFilenames.begin(), selectedFilenames.end());
		closedir(Dir);
	}
	
	//given a list of "object names" that are the substrings of the filenames
	//and a folder path that includes the files that we want to obtain the paths for,
	//this function returns the list of the filenames corresponding to the
	//object list in a string vector 
	//while preserving the order of the files as they appear in the @objectListFilename
	static void loadFilePathFromFolderSelectivelyPreservingOrder(std::string objectListFilename, std::string folderPath, std::vector<std::string> &selectedFilenames, bool indicateFileNotFound=false, std::string fileNotFoundText="")
	{
		std::vector<std::string> objectNames,selectedFilenamesTemp;
		loadFileNamesFromList(objectListFilename,objectNames);

		unsigned char isFile =0x8;//isFolder =0x4;
		DIR *Dir;
		struct dirent *DirEntry;
		std::string filePath;
		
		for(std::vector<std::string>::iterator filename=objectNames.begin();filename!=objectNames.end();filename++)
		{
			if( (Dir = opendir(folderPath.c_str())) == NULL)
         {
            std::cerr << "Error: Path not found in loadFilePathFromFolderSelectivelyPreservingOrder()... Path: "<<folderPath<<std::endl;
            exit(1);
         }

         bool fileNotFound=true;

         while(DirEntry=readdir(Dir))
         {
            if ( DirEntry->d_type == isFile && (std::string(DirEntry->d_name)).find(*filename) != std::string::npos)
            {
               filePath = folderPath;
               filePath.append(DirEntry->d_name);
               selectedFilenames.push_back(filePath);
               fileNotFound=false;
               break;
            }
         }
         if(fileNotFound==true && indicateFileNotFound==true)
            selectedFilenames.push_back(fileNotFoundText);
         closedir(Dir);
		}
	}
	
	
	//given a string that contains the complete path of a file, 
	//this function prunes the trailing folder path and the extension at the end
	//and returns the filename
	static inline std::string pruneFilenameFromPath(std::string filePath)
	{
		std::string filename = filePath.substr(filePath.find_last_of("/")+1);
		return filename.substr(0,filename.find_last_of("."));
	}

	//-given the path of a source folder, a list of filenames to be copied from this folder, and the target folder to copy the reduced set,
	//       this function only copies the files that are included in the list to the target folder
	//-reducedSetList should contain one filename per line. Only the first word at each line is taken into account. The rest is ignored.
	//-extension of the filenames are ignored. You can take extensions into account by commenting the two line where we prune the extensions 
	//       from cellName and filename below.
	static void selectivelyCopyFiles(std::string srcDatasetPath, std::string reducedSetList, std::string targetDatasetPath)
	{
		unsigned char isFile =0x8, isFolder =0x4;
		DIR *Dir;
		struct dirent *DirEntry;

		std::ifstream file; 
		std::istringstream is; 
		std::string line, cellName;  

		file.open(reducedSetList.c_str());
		if(file.fail())
		{
			std::cerr << "File not found in selectivelyCopyFiles():"<<reducedSetList<<std::endl;
			exit(1);
		}	

		while(getline(file, line))
		{
		   is.clear();
		   is.str(line);
		   is >> cellName;
		   cellName.erase(cellName.find("."),4);//erase extension from the cellName

		   Dir = opendir(srcDatasetPath.c_str());
		   while(DirEntry=readdir(Dir))
		   {
				if ( DirEntry->d_type == isFile)
				{
					std::string filename = DirEntry->d_name;
					filename.erase(filename.find("."),4);//erase extension from the fileName
					if(filename.compare(cellName)==0)
					{
						std::cout<<filename<<std::endl;
						system(("cp "+srcDatasetPath+"/"+DirEntry->d_name+" "+targetDatasetPath+"/"+DirEntry->d_name).c_str());
						break;
					}
				}
		   }
		   closedir(Dir);
		}
		file.close();
	}

	static bool writeTGA(const char *file, short int width, short int height, unsigned char *outImage)
	{
		// To save a screen shot is just like reading in a image.  All you do
		// is the opposite.  Istead of calling fread to read in data you call
		// fwrite to save it.

		FILE *pFile;               // The file pointer.
		unsigned char uselessChar; // used for useless char.
		short int uselessInt;      // used for useless int.
		unsigned char imageType;   // Type of image we are saving.
		int index;                 // used with the for loop.
		unsigned char bits;    // Bit depth.
		long Size;                 // Size of the picture.
		int colorMode;
		unsigned char tempColors;

		// Open file for output.
		pFile = fopen(file, "wb");

		// Check if the file opened or not.
		if(!pFile) { fclose(pFile); return false; }

		// Set the image type, the color mode, and the bit depth.
		imageType = 2; colorMode = 3; bits = 24;

		// Set these two to 0.
		uselessChar = 0; uselessInt = 0;

		// Write useless data.
		fwrite(&uselessChar, sizeof(unsigned char), 1, pFile);
		fwrite(&uselessChar, sizeof(unsigned char), 1, pFile);

		// Now image type.
		fwrite(&imageType, sizeof(unsigned char), 1, pFile);

		// Write useless data.
		fwrite(&uselessInt, sizeof(short int), 1, pFile);
		fwrite(&uselessInt, sizeof(short int), 1, pFile);
		fwrite(&uselessChar, sizeof(unsigned char), 1, pFile);
		fwrite(&uselessInt, sizeof(short int), 1, pFile);
		fwrite(&uselessInt, sizeof(short int), 1, pFile);

		// Write the size that you want.
		fwrite(&width, sizeof(short int), 1, pFile);
		fwrite(&height, sizeof(short int), 1, pFile);
		fwrite(&bits, sizeof(unsigned char), 1, pFile);

		// Write useless data.
		fwrite(&uselessChar, sizeof(unsigned char), 1, pFile);

		// Get image size.
		Size = width * height * colorMode;

		// Now switch image from RGB to BGR.
		for(index = 0; index < Size; index += colorMode)
		{
			tempColors = outImage[index];
			outImage[index] = outImage[index + 2];
			outImage[index + 2] = tempColors;
		}

		// Finally write the image.
		fwrite(outImage, sizeof(unsigned char), Size, pFile);

		// close the file.
		fclose(pFile);

		return true;
	}

	//</editor-fold>

	// <editor-fold defaultstate="collapsed" desc=" Memory related functions: allocate2Dmemory(), allocate3Dmemory(), free2DMemory(), free3DMemory() ">
	template <typename Type>
	static Type** allocate2Dmemory(int i, int j)
	{
		Type** data = new Type*[i];
		for (int p=0;p<i;p++)
		{
			data[p] = new Type[j];
			for(int q=0;q<j;q++)
				data[p][q]=0;
		}
		return data;
	}
	
	template <typename Type>
	static void free2Dmemory(Type **data,int rowCount)
	{
		for(int i=0;i<rowCount;i++)
			delete[] data[i];
		delete[] data;
	}
	
	template <typename Type>
	static Type*** allocate3Dmemory(int i, int j, int k)
	{
		Type*** data = new Type**[i];
		for (int p=0;p<i;p++)
		{
			data[p] = new Type*[j];
			for(int q=0;q<j;q++)
			{
				data[p][q] = new Type[k];
				for(int r=0;r<k;r++)
					data[p][q][r]=0;
			}
		}
		return data;
	}
	
	template <typename Type>
	static void free3Dmemory(Type ***data,int i, int j)
	{
		for(int p=0;p<i;p++)
		{
			for(int q=0;q<j;q++)
			{
				delete[] data[p][q];
			}
			delete[] data[p];
		}
		delete[] data;
	}
	//</editor-fold>
	
	// <editor-fold defaultstate="collapsed" desc=" Vector utility functions: sortValuesOfVector(), calculateHistogram(), fillVector(), printVector()X2 , printVectorOfPairs(), print1DArray(), printMap()">
   //struct to be used for holding the elements of a matrix as value/index pairs.
	struct VectorValueIndex
	{
	   float value;
	   int i;
	   enum SortOrder{ASCENDING=0,DESCENDING};

	   VectorValueIndex(float _value,int _i){value=_value;i=_i;}

	   static inline bool greaterThanFunc(const VectorValueIndex &x, const VectorValueIndex &y){return (x.value>y.value);};
	   static inline bool lessThanFunc(const VectorValueIndex &x, const VectorValueIndex &y){return (x.value<y.value);};
	};
	
	//given a vector,sort them in ascending/descending order,and then return the sorted vector with the indices of the original vector.
	template <typename Type>
	static std::vector<VectorValueIndex> sortValuesOfVectorWithIndices(std::vector<Type> &vec, enum VectorValueIndex::SortOrder order)
	{
		std::vector<VectorValueIndex> sorted;
		int size = vec.size();
		for(int i=0;i<size;i++)
			sorted.push_back(VectorValueIndex(vec[i],i));

		if(order==VectorValueIndex::ASCENDING)
			std::sort(sorted.begin(),sorted.end(),VectorValueIndex::lessThanFunc);//greaterFunc is defined above
		else
			std::sort(sorted.begin(),sorted.end(),VectorValueIndex::greaterThanFunc);//greaterFunc is defined above
		
		return sorted;
	}
   
   //given a vector,sort them in ascending/descending order,and then return the vector.
	template <typename Type>
	static std::vector<Type> sortValuesOfVector(std::vector<Type> &vec, enum VectorValueIndex::SortOrder order)
	{
		std::vector<VectorValueIndex> sorted;
		int size = vec.size();
		for(int i=0;i<size;i++)
			sorted.push_back(VectorValueIndex(vec[i],i));

		if(order==VectorValueIndex::ASCENDING)
			std::sort(sorted.begin(),sorted.end(),VectorValueIndex::lessThanFunc);//greaterFunc is defined above
		else
			std::sort(sorted.begin(),sorted.end(),VectorValueIndex::greaterThanFunc);//greaterFunc is defined above
		
      std::vector<Type> sortedVector;
      sortedVector.reserve(size);
      for(std::vector<VectorValueIndex>::iterator iter=sorted.begin();iter!=sorted.end();iter++)
         sortedVector.push_back(iter->value);
      
		return sortedVector;
	}
   
   //given a vector of numbers, calculate a histogram
   //@vec: a vector of numbers
   //@numBins: number of bins. If not specified, number of bins will be set to the largest integer (rounded up) in the input vector
   //@ignoreZero: by default, zeros in the array are ignored in the calculation of the histogram. Set false to account for zeros.
   template <typename Type>
   static inline std::vector<int> calculateHistogram(std::vector<Type> &vec, int numBins=-1,bool ignoreZeros=true)
   {
      std::vector<Type> sortedVector = sortValuesOfVector<Type>(vec,VectorValueIndex::ASCENDING);
      sortedVector.pop_back();// for whatever reason, sort function is adding an extra element to the end of this sorted list. I couldn't figure out the reason and temporarily(!) putting this fix here.
      
      int size = sortedVector.size();
      float range = sortedVector[size]; //set the range as [0,maxValue]
      float binSize;
      if(numBins==-1)
      {
         numBins = Utility::roundUp(range);
         binSize = 1;
      }
      else
         binSize = float(range) / float(numBins);
      
      std::vector<int> histogram(numBins,0);
      
      for(typename std::vector<Type>::iterator iter=sortedVector.begin();iter!=sortedVector.end();iter++)
      {
         if(*iter < Utility::EPSILON && ignoreZeros==true)
            continue;
         if(*iter < Utility::EPSILON && ignoreZeros==false)
            histogram[0]++;
         else
         {
            int bin = int(*iter / binSize);
            histogram[bin]++;
         }
      }
      return histogram;
   }
	
	template <typename Type>
	static inline Type getMaxValueOfVector(Type *vec, int size)
	{
		Type max = vec[0];
		for(int i=0;i<size;i++)
		{
			if(max < vec[i])
				max = vec[i];
		}
		return max;
	}
	
	template <typename Type>
	static inline Type getMaxValueOfVector(std::vector<Type> &vec)
	{
		int size = vec.size();
		Type max = vec[0];
		for(int i=0;i<size;i++)
		{
			if(max < vec[i])
				max = vec[i];
		}
		return max;
	}
	
	template <typename Type>
	static inline Type getMinValueOfVector(Type *vec, int size)
	{
		Type min = vec[0];
		for(int i=0;i<size;i++)
		{
			if(min > vec[i])
				min = vec[i];
		}
		return min;
	}
	
	template <typename Type>
	static inline Type getMinValueOfVector(std::vector<Type> &vec)
	{
		int size = vec.size();
		Type min = vec[0];
		for(int i=0;i<size;i++)
		{
			if(min > vec[i])
				min = vec[i];
		}
		return min;
	}
	
	template <typename Type>
	static inline Type getMedianValueOfVector(std::vector<Type> &vec)
	{
		int size = vec.size();
		
		Type median = sortValuesOfVector<Type>(vec,VectorValueIndex::ASCENDING)[size/2];
		
		return median;
	}
	
	
	template <typename Type>
	static inline void fillVector(std::vector<Type> &vec, Type value)
	{
		for(typename std::vector<Type>::iterator iter=vec.begin();iter!=vec.end();iter++)
			*iter = value;
	}
	
	template <typename Type>
	static inline void fillVector(std::vector<Type> &vec, Type value,int size)
	{
		vec.clear();
		for(int i=0;i<size;i++)
			vec.push_back(value);
	}
	
	template <typename Type>
	static inline void fillVector(Type *vec, Type value, int size)
	{
		for(int i=0;i<size;i++)
			vec[i] = value;
	}
	
	template <typename Type>
	static inline void elementwiseAddVectors(std::vector<Type> vec1, std::vector<Type> vec2,std::vector<Type> &vecDst)
	{
		int size = vec1.size();
		if(size!=vec2.size())
		{
			std::cerr<<"attempting to elementwise add vectors of unequal size!!! Exiting...\n";
			exit(1);
		}
		
		for(int i=0;i<size;i++)
			vecDst[i] = vec1[i]+vec2[i];
		return;
	}
	
	template <typename Type>
	static inline std::vector<Type> elementwiseAddVectors(std::vector<Type> vec1, std::vector<Type> vec2)
	{
		int size = vec1.size();
		if(size!=vec2.size())
		{
			std::cerr<<"attempting to elementwise add vectors of unequal size!!! Exiting...\n";
			exit(1);
		}
		std::vector<Type> vecDst(size,0);

		for(int i=0;i<size;i++)
			vecDst[i] = vec1[i]+vec2[i];
		return vecDst;
	}
	
	template <typename Type1, typename Type2>
	static inline void multiplyVectorWithScalar(std::vector<Type1> &vec, Type2 scalar)
	{
		for(typename std::vector<Type1>::iterator iter=vec.begin();iter!=vec.end();iter++)
			*iter *= scalar;
	}
	
	template <typename Type>
	static void removeElementWithValueFromVector(std::vector<Type> &vec, Type value)
	{
		for(typename std::vector<Type>::iterator iter=vec.begin();iter!=vec.end();iter++)
		{
			if(*iter<=(value+EPSILON) && *iter>=(value-EPSILON))
			{
				vec.erase(iter);
				break;
			}
		}
		return;
	}
	
	template <typename Type>
	static int returnIndexOfMaxElement(std::vector<Type> &vec)
	{
		Type maxValue = vec[0];
		int maxIndex = 0;
		int index = 0;
		for(typename std::vector<Type>::iterator iter=vec.begin();iter!=vec.end();iter++)
		{
			if(*iter > maxValue)
			{
				maxValue = *iter;
				maxIndex = index;
			}
			index++;
		}
		return maxIndex;
	}
	
	//given a (potentially empty) vector @vec, this function randomly fills up the vector
	//with @dimension number of float values within the range [0,@upperRange]
	static void randomlySetVector(std::vector<float> &vec, int dimension, int upperRange)
	{
		std::srand((int)clock());
		vec.clear();
		vec.reserve(dimension);
		for(int i=0;i<dimension;i++)
			vec.push_back(std::rand()%upperRange + (rand()%100)/100.0);
	}
	
	static void randomlySetVector(std::vector<float> &vec, int dimension, float lowerRange, float upperRange)
	{
		std::srand((int)clock());
		vec.clear();
		vec.reserve(dimension);
		for(int i=0;i<dimension;i++)
			vec.push_back(generateRandomNumber(std::rand(),lowerRange,upperRange));
	}

	//given a (potentially empty) vector @vec, this function randomly fills up the vector
	//with @dimension number of values that are hard coded inside the vector @values within the function
	static void randomlySetVectorDiscrete(std::vector<float> &vec, std::vector<float> &values, int dimension)
	{
		int numOfValues = values.size(); //size of the array above
		std::srand((int)clock());
		vec.clear();
		vec.reserve(dimension);
		for(int i=0;i<dimension;i++)
			vec.push_back(values[std::rand()%numOfValues]);
	}
	
	//given a file that contains a vector in a single line, this function reads the contents and 
	//saves them into a vector
//	template <typename Type>
//	static inline void loadVectorFromFile(std::string filename, std::vector<Type> &vec)
//	{
//		std::ifstream file; 
//		std::istringstream is; 
//		std::string line;
//
//		////load names of result files
//		file.open(filename.c_str());
//		if(file.fail())
//		{
//			std::cout << "File not found!!!!\n";
//			return;
//		}	
//
//		getline(file, line);
//		is.clear();
//		is.str(line);
//
//		vec.clear();
//		vec.reserve(size);
//		Type tempVal;
//		
//		while(is >> tempVal)
//			vec.push_back(tempVal);
//		
//		file.close();
//	}
	
	//given a file containing values, one value at each line, this function reads them
	//and saves the values into a vector
	template <typename Type>
	static inline void loadVectorFromFile(std::string filename, std::vector<Type> &vec)
	{
		std::ifstream file; 
		std::istringstream is; 
		std::string line;
		Type tempVal;

		file.open(filename.c_str());
		if(file.fail())
		{
			std::cout << "File not found in loadVectorFromFile():"<<filename<<std::endl;
			return;
		}	
		
		//count rows and columns, this info will be used to read the vector either from rows or from columns
		int numRow=0,numCol=0;
		//count rows
		while(getline(file, line))
			numRow++;
		file.clear();//reset the file pointer to the beginning of the file
		file.seekg(0, std::ios::beg);
		
		//count columns
		getline(file, line);
		is.clear();
		is.str(line);
		while(is >> tempVal)
			numCol++;
		file.clear();//reset the file pointer to the beginning of the file
		file.seekg(0, std::ios::beg);
		
		if(numRow>1 && numCol>1)//if the file contains multiple rows and columns, then print error message
		{
			std::cerr<<"trying to load a matrix from file into a vector in loadVectorFromFile()!!! Exiting...\n";
			exit(1);
		}
		
		vec.clear();
		if(numRow>1)//read column wise
		{
			vec.reserve(numRow);
			while(getline(file, line))
			{
				is.clear();
				is.str(line);
				is >> tempVal;
				vec.push_back(tempVal);
			}
		}
		else if(numCol>1)//read row wise
		{
			vec.reserve(numCol);
			getline(file, line);
			is.clear();
			is.str(line);
			while(is >> tempVal)
				vec.push_back(tempVal);
		}
		file.close();
	}
	
	static inline bool doesVectorIncludeValue(std::vector<int> &vec, int value)
	{
		for(std::vector<int>::iterator iter=vec.begin();iter!=vec.end();iter++)
			if(*iter==value)
				return true;
			
		return false;
	}
	
	template <typename Type>
	static inline void loadVectorFromFileAfterSkipLine(std::string filename, std::vector<Type> &vec, int size, int skipLine)
	{
		std::ifstream file; 
		std::istringstream is; 
		std::string line;

		////load names of result files
		file.open(filename.c_str());
		if(file.fail())
		{
			std::cout << "File not found in loadVectorFromFileAfterSkipLine():"<<filename<<std::endl;
			return;
		}	

		for(int i=0;i<skipLine;i++)
			getline(file,line);
		
		getline(file, line);
		is.clear();
		is.str(line);

		vec.clear();
		vec.reserve(size);
		Type tempVal;
		
		for(int i=0;i<size;i++)
		{
			is >> tempVal;
			vec.push_back(tempVal);
		}

		file.close();
	}
	
	template <typename Type>
	static inline void saveVectorToFile(std::string filename, std::vector<Type> &vec, char separator='\t', char endOfLine='\n', int precision=2)
	{
		std::ofstream file;
		file.open(filename.c_str());
		file.setf(std::ios::fixed,std::ios::floatfield);
		file.precision(precision);
   
		for(typename std::vector<Type>::iterator iter=vec.begin();iter!=vec.end();iter++)
			file<<*iter<<separator;
		file<<endOfLine;
		
		file.close();
	}
	
	template <typename Type>
	static inline void appendVectorToFile(std::string filename, std::vector<Type> &vec, char separator='\t', char endOfLine='\n', int precision=2)
	{
		std::ofstream file;
		file.open(filename.c_str(),std::ios_base::app);
		file.setf(std::ios::fixed,std::ios::floatfield);
		file.precision(precision);
   
		for(typename std::vector<Type>::iterator iter=vec.begin();iter!=vec.end();iter++)
			file<<*iter<<separator;
		file<<endOfLine;
		
		file.close();
	}
	
	template <typename Type>
	static inline void printVector(Type *vec, int size, std::ostream &out=std::cout, char separator='\t', char endOfLine='\n', int precision=2)
	{
		out.precision(precision);
		for(int i=0;i<size;i++)
			out<<vec[i]<<separator;
		out<<endOfLine;
	}

	template <typename Type>
	static inline void printVector(std::vector<Type> vec, std::ostream &out=std::cout, char separator='\t', char endOfLine='\n', int precision=2)
	{
		out.precision(precision);
		for(int i=0;i<vec.size();i++)
			out<<vec[i]<<separator;
		out<<endOfLine;
	}
	
	template <typename Type1, typename Type2>
	static inline void printVectorOfPairs(std::vector<std::pair<Type1,Type2> > myVector, std::ostream &out=std::cout, char endOfLine='\n')
	{
		for(typename std::vector<std::pair<Type1,Type2> >::iterator iter=myVector.begin();iter!=myVector.end();iter++)
			out<<iter->first<<"\t"<<iter->second<<endOfLine;
	}
	
	template <typename Type>
	static inline void print1DArray(Type *arr, int size, std::ostream &out=std::cout, char separator='\t', char endOfLine='\n', int precision=2)
	{
		out.precision(precision);
		for(int i=0;i<size;i++)
			out<<arr[i]<<separator;
		out<<endOfLine;
	}
	
	template <typename Type1, typename Type2>
	static inline void printMap(std::map<Type1,Type2> myMap, std::ostream &out=std::cout, char endOfLine='\n')
	{
		for(typename std::map<Type1,Type2>::iterator iter=myMap.begin();iter!=myMap.end();iter++)
			out<<iter->first<<"\t"<<iter->second<<endOfLine;
	}
	
	// </editor-fold>
	
	// <editor-fold defaultstate="collapsed" desc=" Matrix utility functions: copyMatrix(), fillMatrix(), logScaleMatrix(), saveMatrixToFile(), loadMatrixFromFile(), printMatrix()x2, getMaxValueOfMatrix()x2, getMinValueOfMatrix()x2">
	//struct to be used for holding the elements of a matrix as value/index pairs.
	struct MatrixValueIndex
	{
	   float value;
	   int i,j;
	   enum SortOrder{ASCENDING=0,DESCENDING};

	   MatrixValueIndex(float _value,int _i, int _j){value=_value;i=_i;j=_j;}

	   static inline bool greaterThanFunc(const MatrixValueIndex &x, const MatrixValueIndex &y){return (x.value>y.value);};
	   static inline bool lessThanFunc(const MatrixValueIndex &x, const MatrixValueIndex &y){return (x.value<y.value);};
	};
	
	//given a square matrix, push the values contained in the upper triangular matrix into a vector of value/index structure, 
	//sort them in ascending/descending order,and then return the vector with the original indices of the sorted elements.
	template <typename Type>
	static std::vector<MatrixValueIndex> sortValuesOfUpperTriangularMatrixWithIndices(Type** matrix, int row, enum MatrixValueIndex::SortOrder order)
	{
		std::vector<MatrixValueIndex> sorted;
		for(int i=0;i<row;i++)
			for(int j=i+1;j<row;j++)
				sorted.push_back(MatrixValueIndex(matrix[i][j],i,j));

		if(order==MatrixValueIndex::ASCENDING)
			std::sort(sorted.begin(),sorted.end(),MatrixValueIndex::lessThanFunc);//greaterFunc is defined above
		else
			std::sort(sorted.begin(),sorted.end(),MatrixValueIndex::greaterThanFunc);//greaterFunc is defined above
		
		return sorted;
	}
   
   //given a square matrix, push the values contained in the upper triangular matrix into a vector of value/index structure, 
	//sort them in ascending/descending order,and then return the vector.
	template <typename Type>
	static std::vector<Type> sortValuesOfUpperTriangularMatrix(Type** matrix, int row, enum MatrixValueIndex::SortOrder order)
	{
		std::vector<MatrixValueIndex> sorted;
		for(int i=0;i<row;i++)
			for(int j=i+1;j<row;j++)
				sorted.push_back(MatrixValueIndex(matrix[i][j],i,j));

		if(order==MatrixValueIndex::ASCENDING)
			std::sort(sorted.begin(),sorted.end(),MatrixValueIndex::lessThanFunc);//greaterFunc is defined above
		else
			std::sort(sorted.begin(),sorted.end(),MatrixValueIndex::greaterThanFunc);//greaterFunc is defined above
		
		std::vector<Type> sortedVector;
		sortedVector.reserve(sorted.size());
		for(typename std::vector<MatrixValueIndex>::iterator iter=sorted.begin();iter!=sorted.end();iter++)
		   sortedVector.push_back(iter->value);
      
		return sortedVector;
	}
   
   //given a matrix, pushes its values into a vector, sorts them in ascending/descending order, and returns the sorted vector.
   template <typename Type>
	static std::vector<Type> sortValuesOfMatrix(Type** matrix, int row,int col, enum MatrixValueIndex::SortOrder order)
	{
		std::vector<MatrixValueIndex> sorted;
		for(int i=0;i<row;i++)
			for(int j=0;j<row;j++)
				sorted.push_back(MatrixValueIndex(matrix[i][j],i,j));

		if(order==MatrixValueIndex::ASCENDING)
			std::sort(sorted.begin(),sorted.end(),MatrixValueIndex::lessThanFunc);//greaterFunc is defined above
		else
			std::sort(sorted.begin(),sorted.end(),MatrixValueIndex::greaterThanFunc);//greaterFunc is defined above
		
		std::vector<Type> sortedVector;
		sortedVector.reserve(sorted.size());
		for(typename std::vector<MatrixValueIndex>::iterator iter=sorted.begin();iter!=sorted.end();iter++)
		   sortedVector.push_back(iter->value);
      
		return sortedVector;
	}
   
   //given a matrix of numbers, calculate a histogram
   //@matrix: a matrix of numbers
   //@numBins: number of bins. If not specified, number of bins will be set to the largest integer (rounded up) in the input vector
   //@ignoreZero: by default, zeros in the array are ignored in the calculation of the histogram. Set false to account for zeros.
   //@upperTriangular: by default, histogram is calculated for the entire matrix. Set true to upper triangular to calculate the histogram only for the upper triangular matrix. (note: matrix should be square for this option)
   template <typename Type>
   static inline std::vector<int> calculateHistogram(Type** matrix, int row,int col,int numBins=-1,bool ignoreZeros=true, bool upperTriangular=false)
   {
      if(upperTriangular==true)
      {
         if(row!=col)
         {
            std::cerr<<"calculateHistogram() is called for a non-square matrix with upperTriangular=true parameter!! Exiting..\n";
            exit(1);
         }
         std::vector<Type> upperTriangularMatrix = Utility::copyOffDiagonalUpperTriangularOfMatrixIntoVector(matrix,row);
         return Utility::calculateHistogram(upperTriangularMatrix,numBins,ignoreZeros);
      }
      else
      {
         std::vector<Type> vectorizedMatrix = Utility::copyMatrixTo1DVector<Type>(matrix,row,col);
         std::vector<Type> sortedVec = sortValuesOfVector<Type>(vectorizedMatrix,VectorValueIndex::ASCENDING);
         
         sortedVec.pop_back();// for whatever reason, sort function is adding an extra element to the end of this sorted list. I couldn't figure out the reason and temporarily(!) putting this fix here.
      
         int size = sortedVec.size();
         float range = sortedVec[size]; //set the range as [0,maxValue]

         float binSize;
         if(numBins==-1)
         {
            numBins = Utility::roundUp(range);
            binSize = 1;
         }
         else
            binSize = float(range) / float(numBins);

         std::vector<int> histogram(numBins,0);

         for(typename std::vector<Type>::iterator iter=sortedVec.begin();iter!=sortedVec.end();iter++)
         {
            if(*iter < Utility::EPSILON && ignoreZeros==true)
               continue;
            if(*iter < Utility::EPSILON && ignoreZeros==false)
               histogram[0]++;
            else
            {
               int bin = int(*iter / binSize);
               histogram[bin]++;
            }
         }
         return histogram;
      }
      
   }
	
	//given a square matrix, returns the ratio of number of non-diagonal entries that are non-zero/greater than zero/less-than zero to total number of non-diagonal entries
   //@type: takes values absolute/positive/negative, where the density for non-zero, positive, and negative values are calculated, respectively, Default behavior is "full"
	template <typename Type>
	static inline float calculateMatrixDensity(Type** matrix, int row, std::string type="absolute")
	{
		if(matrix==NULL)
		{
			std::cerr<<"calculateMatrixDensity() is called for an empty matrix!! Exiting..\n";
			exit(1);
		}

		int count = 0;
		for(int i=0;i<row;i++)
			for(int j=0;j<row;j++)
			{
				if(type=="absolute" && absoluteValue(matrix[i][j])>0 && i!=j)
					count++;
				else if(type=="positive" && matrix[i][j]>0 && i!=j)
				   count++;
				else if(type=="negative" && matrix[i][j]<0 && i!=j)
				   count++;
			}
		float density = ((float)count) / ((float) (row*(row-1)));
		return density;
	}
        
        //given a square matrix, returns the summation of edges in total connectome, inter hemispheric connections or intra hemispheric connections. Could give positive and negative connections separately
        //@connectivitySign: takes values positive/negative, where the connectivity strength for positive, and negative values are calculated, respectively, Default behavior is "positive"
        //@connectivityType: takes values total/interhemispheric/intrahemispheric, where the density for non-zero, positive, and negative values are calculated, respectively, Default behavior is "absolute"
        //@hemispheres: integer vector to indicate left/right hemispheres per node. -1 indicates cerebellum (or any non-left/right hemispheric ROI)
	template <typename Type>
	static inline float calculateMatrixConnectivityStrength(Type** matrix, int row, std::vector<int> hemispheres, std::string connectivityType="total", std::string connectivitySign="positive")
	{
		if(matrix==NULL)
		{
			std::cerr<<"calculateMatrixConnectivityStrength() is called for an empty matrix!! Exiting..\n";
			exit(1);
		}

		float connectivity=0;
		if(connectivitySign=="positive")
		{
			if(connectivityType=="total")
			{
				for(int i=0;i<row;i++)
					for(int j=0;j<row;j++)
						if(matrix[i][j]>0)
							connectivity+=matrix[i][j];
			}
			else if(connectivityType=="interhemispheric")
			{
				for(int i=0;i<row;i++)
					for(int j=0;j<row;j++)
						if(matrix[i][j]>0 && hemispheres[i]==hemispheres[j])
							connectivity+=matrix[i][j];
			}
			else if(connectivityType=="intrahemispheric")
			{
				for(int i=0;i<row;i++)
					for(int j=0;j<row;j++)
						if(matrix[i][j]>0 && hemispheres[i]!=hemispheres[j])
							connectivity+=matrix[i][j];
			}
		}
		else if(connectivitySign=="negative")
		{
			if(connectivityType=="total")
			{
				for(int i=0;i<row;i++)
					for(int j=0;j<row;j++)
						if(matrix[i][j]<0)
							connectivity+=matrix[i][j];
			}
			else if(connectivityType=="interhemispheric")
			{
				for(int i=0;i<row;i++)
					for(int j=0;j<row;j++)
						if(matrix[i][j]<0 && hemispheres[i]==hemispheres[j])
							connectivity+=matrix[i][j];
			}
			else if(connectivityType=="intrahemispheric")
			{
				for(int i=0;i<row;i++)
					for(int j=0;j<row;j++)
						if(matrix[i][j]<0 && hemispheres[i]!=hemispheres[j])
							connectivity+=matrix[i][j];
			}
		}


		connectivity /= 2.0;
		return connectivity;
	}
   
   //given a matrix, this function determines the positive/negativeness of the elements of the matrix.
   //@matrix: input matrix
   //@row: number of rows of the matrix
   //@col: number of columns of the matrix
   //@signMatrix: the output will be written into this matrix, -1 indicating negative, 0 indicating 0 and 1 indicating positive.
   template <typename Type>
   static inline void getSignsOfMatrixElements(Type** matrix, int row,int col,int **signMatrix)
   {
		for(int i=0;i<row;i++)
		{
			for(int j=0;j<col;j++)
			{
				if(matrix[i][j]<-Utility::EPSILON)
				   signMatrix[i][j]=-1;
				else if (matrix[i][j]>-Utility::EPSILON && matrix[i][j]<Utility::EPSILON)
				   signMatrix[i][j]=0;
				else if(matrix[i][j]>Utility::EPSILON)
				   signMatrix[i][j]=1;
				else
				{
				   std::cerr<<"Unexpected value in the matrix encountered in Utility::getSignsOfMatrixElements() function. Exiting!!!\n";
				   exit(1);
				}
			}
		}
   }
   
   template <typename Type>
	static inline void absoluteValue(Type** srcMatrix, Type** dstMatrix, int row,int col)
	{
		for(int i=0;i<row;i++)
		{
			for(int j=0;j<col;j++)
			{
				dstMatrix[i][j] = (srcMatrix[i][j] < 0 ? -srcMatrix[i][j] : srcMatrix[i][j]);
			}
      }
	}
	
	//given a matrix consisting of correlation values (i.e., r), calculates the Fisher's Z transform for each of its values and updates the matrix in place with calculated results
   //formula is: z = artanh(r) =1/2*ln((1+r)/(1-r)) where r is the correlation between two vectors of data
	template <typename Type>
	static inline void fisherZTransformMatrix(Type** matrix, int row, int col)
	{
		for(int i=0;i<row;i++)
		{
			for(int j=0;j<col;j++)
			{
				float nominator = (float) (1+matrix[i][j]);
				float denominator = 1-matrix[i][j];
				denominator = (absoluteValue(denominator) < EPSILON ? denominator+EPSILON : denominator);

				matrix[i][j] = 0.5*std::log(nominator/denominator);
			}
		}
	}

	//copies contents of the source matrix to the destination matrix
	template <typename Type>
	static inline void copyMatrix(Type** src, Type** dst, int row, int column)
	{
		if(src==NULL || dst==NULL)
		{
			std::cerr<<"copyMatrix() is called for an empty source/destination matrix!! Exiting..\n";
			exit(1);
		}

		for(int i=0;i<row;i++)
			for(int j=0;j<column;j++)
				dst[i][j]=src[i][j];
	}
	
	template <typename Type1, typename Type2>
	static inline void copyMatrix(Type1** src, Type2** dst, int row, int column)
	{
		if(src==NULL || dst==NULL)
		{
			std::cerr<<"copyMatrix() is called for an empty source/destination matrix!! Exiting..\n";
			exit(1);
		}
		
		for(int i=0;i<row;i++)
			for(int j=0;j<column;j++)
				dst[i][j] = (Type2)src[i][j];
	}
	
	template <typename Type>
	static inline void copyMatrix(Type*** src, Type*** dst, int axis1, int axis2, int axis3)
	{
		if(src==NULL || dst==NULL)
		{
			std::cerr<<"copyMatrix() is called for an empty source/destination matrix!! Exiting..\n";
			exit(1);
		}

		for(int i=0;i<axis1;i++)
			for(int j=0;j<axis2;j++)
				for(int k=0;k<axis3;k++)
					dst[i][j][k]=src[i][j][k];
	}
	
   //given a matrix, copy the elements of the matrix into a 1 dimensional array line by line, top to bottom
	template <typename Type>
	static inline Type* copyMatrixTo1DArray(Type** src, int rows, int cols) 
	{
		Type *array = new Type[rows*cols];

		for(int i=0;i<rows;i++)
		{
			for(int j=0;j<cols;j++)
				array[i*rows+j] = src[i][j];
		}
		return array;
	}
   
   //given a matrix, copy the elements of the matrix into a 1 dimensional vector line by line, top to bottom
   template <typename Type>
	static inline std::vector<Type> copyMatrixTo1DVector(Type** src, int rows, int cols) 
	{
		std::vector<Type> vec(rows*cols);

		for(int i=0;i<rows;i++)
		{
			for(int j=0;j<cols;j++)
				vec.push_back(src[i][j]);
		}
		return vec;
	}
	
	//given a 3 dimensional matrix, copy a row of a matrix at a given axis
	//for example,  if the axis is 0, then we return the elements contained in [i,j]th index of each matrix
	//that is, matrix[k][i][j] for k=0..sizeFixedAxis
	template <typename Type>
	static inline std::vector<Type> copyRowOfMatrixIntoVector(Type*** matrix, int sizeFixedAxis, int i, int j, int fixedAxis)
	{
		std::vector<Type> vec;
		vec.reserve(sizeFixedAxis);
		if(matrix==NULL)
		{
			std::cerr<<"copyRowOfMatrixIntoVector() is called for an empty source matrix!! Exiting..\n";
			exit(1);
		}

		if(fixedAxis==0)
		{
			for(int p=0;p<sizeFixedAxis;p++)
				vec.push_back(matrix[p][i][j]);
		}
		else if(fixedAxis==1)
		{
			for(int p=0;p<sizeFixedAxis;p++)
				vec.push_back(matrix[i][p][j]);
		}
		else if(fixedAxis==2)
		{
			for(int p=0;p<sizeFixedAxis;p++)
				vec.push_back(matrix[i][j][p]);
		}
		return vec;
	}
	
	//given a 2 dimensional square matrix, copy the diagonal of the matrix into a vector
	template <typename Type>
	static inline std::vector<Type> copyDiagonalOfMatrixIntoVector(Type** matrix, int size)
	{
		std::vector<Type> vec;
		vec.reserve(size);
		if(matrix==NULL)
		{
			std::cerr<<"copyDiagonalOfMatrixIntoVector() is called for an empty source matrix!! Exiting..\n";
			exit(1);
		}

		for(int p=0;p<size;p++)
			vec.push_back(matrix[p][p]);
		return vec;
	}
	
	//given a 2 dimensional square matrix, copy the off diagonal upper triangular of the matrix into a vector
	template <typename Type>
	static inline std::vector<Type> copyOffDiagonalUpperTriangularOfMatrixIntoVector(Type** matrix, int size)
	{
		std::vector<Type> vec;
		vec.reserve((size*(size-1))/2);
		if(matrix==NULL)
		{
			std::cerr<<"copyDiagonalOfMatrixIntoVector() is called for an empty source matrix!! Exiting..\n";
			exit(1);
		}

		for(int i=0;i<size;i++)
			for(int j=i+1;j<size;j++)
				vec.push_back(matrix[i][j]);
		return vec;
	}
	
	//given a 2D matrix, this function fills in the matrix with the same @value
	template <typename Type>
	static inline void fillMatrix(Type** dst, Type value, int row, int column)
	{
		if(dst==NULL)
		{
			std::cerr<<"fillMatrix() is called for an empty destination matrix!! Exiting..\n";
			exit(1);
		}

		for(int i=0;i<row;i++)
			for(int j=0;j<column;j++)
				dst[i][j]=value;
		return;
	}
	
	//given a 3D matrix, this function fills in the matrix with the same @value
	template <typename Type>
	static inline void fillMatrix(Type*** dst, Type value, int _i, int _j, int _k)
	{
		if(dst==NULL)
		{
			std::cerr<<"fillMatrix() is called for an empty destination matrix!! Exiting..\n";
			exit(1);
		}

		for(int i=0;i<_i;i++)
			for(int j=0;j<_j;j++)
				for(int k=0;k<_k;k++)
					dst[i][j][k]=value;
		return;
	}
	
	//given a square matrix, this function fills the diagonal of the matrix with the @value
	template <typename Type>
	static inline void fillDiagonalOfTheMatrix(Type** dst, Type value, int row)
	{
		if(dst==NULL)
		{
			std::cerr<<"fillDiagonalOfTheMatrix() is called for an empty destination matrix!! Exiting..\n";
			exit(1);
		}

		for(int i=0;i<row;i++)
				dst[i][i]=value;
		return;
	}
	
	//given a square matrix, this function fills the diagonal entries of the matrix with 
	//a random value in the range of [min,max] values within the column of the diagonal entry
	template <typename Type>
	static inline void fillDiagonalOfTheMatrixRandomlyColumnwise(Type** dst, int size, int seed=1)
	{
		if(dst==NULL)
		{
			std::cerr<<"fillDiagonalOfTheMatrix() is called for an empty destination matrix!! Exiting..\n";
			exit(1);
		}

		for(int i=0;i<size;i++)
		{
			dst[i][i]=INF;//while calculating min of the column, ignore the diagonal entry
			float min = Utility::getMinValueOfColumnOfMatrix(dst,i,size);
			dst[i][i]=N_INF;//while calculating max of the column, ignore the diagonal entry
			float max = Utility::getMaxValueOfColumnOfMatrix(dst,i,size);
			
			//generate a random value in the given range
			float randomValue = generateRandomNumber(i*seed,min,max);
			dst[i][i]=randomValue;
		}
				
		return;
	}
	
	//given a square matrix, this function fills the diagonal entries of the matrix with 
	//a random value in the range of [min,max] values within the column of the diagonal entry
	template <typename Type>
	static inline void fillDiagonalOfTheMatrixWithColumnwiseAverage(Type** dst, int size)
	{
		if(dst==NULL)
		{
			std::cerr<<"fillDiagonalOfTheMatrix() is called for an empty destination matrix!! Exiting..\n";
			exit(1);
		}

		for(int i=0;i<size;i++)
		{
			float average=0;
			for(int j=0;j<size;j++)
				average+=dst[i][j];
			average /= (float)size;
					
			dst[i][i]=average;
		}
				
		return;
	}
	
	//given a matrix @dst, this function replaces occurences of @value1 in the matrix with the @value2
	template <typename Type>
	static inline void replaceValueOfMatrix(Type** dst, Type value1, Type value2, int row, int column)
	{
		if(dst==NULL)
		{
			std::cerr<<"replaceValueOfMatrix() is called for an empty destination matrix!! Exiting..\n";
			exit(1);
		}

		for(int i=0;i<row;i++)
			for(int j=0;j<column;j++)
			{
				if(std::abs(value1-dst[i][j])<EPSILON)
					dst[i][j]=value2;
			}
		return;
	}
	
	//given a matrix, this function updates the values of the matrix with the log_2 of all of its elements 
	template <typename Type>
	static inline void logScaleMatrix(Type** dst, int row, int column)
	{
		if(dst==NULL)
		{
			std::cerr<<"logScaleMatrix() is called for an empty destination matrix!! Exiting..\n";
			exit(1);
		}

		for(int i=0;i<row;i++)
			for(int j=0;j<column;j++)
			{
				if(dst[i][j]<=1)
					dst[i][j] = 0;
				else
					dst[i][j] = log2(dst[i][j]);
			}
		return;
	}
	
	//given a matrix, this function updates the values of the matrix with the multiplicative inverse of all of its elements 
	template <typename Type>
	static inline void multiplicativeInverseMatrix(Type** dst, int row, int column)
	{
		if(dst==NULL)
		{
			std::cerr<<"multiplicativeInverseMatrix() is called for an empty destination matrix!! Exiting..\n";
			exit(1);
		}

		for(int i=0;i<row;i++)
			for(int j=0;j<column;j++)
			{
				if(absoluteValue(dst[i][j])<=EPSILON)
					dst[i][j] = INF;
				else if(absoluteValue(dst[i][j])==INF)
					dst[i][j] = 0;
				else
					dst[i][j] = 1.0/dst[i][j];
			}
		return;
	}
		
	//given a matrix, this function returns the sum of its elements
	template <typename Type>
	static inline Type sumElementsOfMatrix(Type** src1, int row, int column)
	{
		Type sum=0;
		if(src1==NULL)
		{
			std::cerr<<"sumElementsOfMatrix() is called for an empty source matrix!! Exiting..\n";
			exit(1);
		}

		for(int i=0;i<row;i++)
			for(int j=0;j<column;j++)
				sum+=src1[i][j];
		return sum;
	}
	
	//given a matrix, this function returns the sum of its elements
	template <typename Type>
	static inline float calculateMatrixDistanceL1(Type** src1, Type** src2, int row, int column)
	{
		Type sum=0;
		if(src1==NULL || src2==NULL)
		{
			std::cerr<<"calculateMatrixDistanceL1() is called for an empty source matrix!! Exiting..\n";
			exit(1);
		}

		for(int i=0;i<row;i++)
			for(int j=0;j<column;j++)
				sum+=absoluteValue(src1[i][j]-src2[i][j]);
		return sum;
	}
	
	//given a matrix, this function returns the sum of its elements
	template <typename Type>
	static inline float calculateMatrixDistanceL2(Type** src1, Type** src2, int row, int column)
	{
		Type sum=0;
		if(src1==NULL || src2==NULL)
		{
			std::cerr<<"calculateMatrixDistanceL1() is called for an empty source/destination matrix!! Exiting..\n";
			exit(1);
		}

		for(int i=0;i<row;i++)
			for(int j=0;j<column;j++)
				sum+=(src1[i][j]-src2[i][j])*(src1[i][j]-src2[i][j]);
		return sqrt(sum);
	}
	
	//given two matrices of the same size, this function multiplies the matrices elementwise 
	//and writes the result into the third matrix
	template <typename Type>
	static inline void elementwiseMultiplyMatrices(Type** src1, Type** src2, Type** dst, int row, int column)
	{
		if(src1==NULL || src2==NULL || dst==NULL)
		{
			std::cerr<<"elementwiseMultiplyMatrices() is called for an empty source/destination matrix!! Exiting..\n";
			exit(1);
		}

		for(int i=0;i<row;i++)
			for(int j=0;j<column;j++)
				dst[i][j]=src1[i][j]*src2[i][j];
		return;
	}
	
	//given two matrices of the same size, this function adds the matrices elementwise 
	//and writes the result into the third matrix
	template <typename Type>
	static inline void elementwiseAddMatrices(Type** src1, Type** src2, Type** dst, int row, int column)
	{
		if(src1==NULL || src2==NULL || dst==NULL)
		{
			std::cerr<<"elementwiseMultiplyMatrices() is called for an empty source/destination matrix!! Exiting..\n";
			exit(1);
		}

		for(int i=0;i<row;i++)
			for(int j=0;j<column;j++)
				dst[i][j]=src1[i][j]+src2[i][j];
		return;
	}
	
	//given two matrices of the same size, this function subtracts the first matrix from the second elementwise 
	//and writes the result into the third matrix
	template <typename Type>
	static inline void elementwiseSubtractMatrices(Type** src1, Type** src2, Type** dst, int row, int column)
	{
		if(src1==NULL || src2==NULL || dst==NULL)
		{
			std::cerr<<"elementwiseMultiplyMatrices() is called for an empty source/destination matrix!! Exiting..\n";
			exit(1);
		}

		for(int i=0;i<row;i++)
			for(int j=0;j<column;j++)
				dst[i][j]=src1[i][j]-src2[i][j];
		return;
	}
	
	//given a matrix and a scalar value, this function multiplies the elements of the @matrix
	//with the @scalar and writes the result into the input @matrix
	template <typename Type>
	static inline void multiplyMatrixWithScalar(Type** matrix, Type scalar, int row, int column)
	{
		if(matrix==NULL)
		{
			std::cerr<<"multiplyMatrixWithScalar() is called for an empty source matrix!! Exiting..\n";
			exit(1);
		}

		for(int i=0;i<row;i++)
			for(int j=0;j<column;j++)
				matrix[i][j]=matrix[i][j]*scalar;
		return;
	}
	
	//given a square matrix and a scalar value, this function multiplies the elements of the @matrix
	//except the diagonal entries with the @scalar and writes the result into the input @matrix
	template <typename Type>
	static inline void multiplyMatrixWithScalarExceptDiagonal(Type** matrix, Type scalar, int row)
	{
		if(matrix==NULL)
		{
			std::cerr<<"multiplyMatrixWithScalar() is called for an empty source matrix!! Exiting..\n";
			exit(1);
		}

		for(int i=0;i<row;i++)
			for(int j=0;j<row;j++)
			{
				if(i==j)
					continue;
				matrix[i][j]=matrix[i][j]*scalar;
			}
		return;
	}
	
	//given two matrices @src1 and @src2 where the number of columns of first matrix 
	//is equal to the number of rows to second matrix, this function performs matrix multiplication  
	//and writes the result into the third matrix @dst
	template <typename Type>
	static inline void multiplyMatrices(Type** src1, Type** src2, Type** dst, int r1, int c1, int c2)
	{
		if(src1==NULL || src2==NULL || dst==NULL)
		{
			std::cerr<<"multiplyMatrices() is called for an empty source/destination matrix!! Exiting..\n";
			exit(1);
		}

		for(int i = 0; i < r1; ++i)
			for(int j = 0; j < c2; ++j)
				dst[i][j]=0;
		
		for(int i=0;i<r1;i++)
			for(int j=0;j<c2;j++)
				for(int k=0;k<c1;k++)
					dst[i][j]+=src1[i][k]*src2[k][j];
		return;
	}
	
	//given a square matrix @src, this function raises the matrix to the indicated power  
	//and writes the result into the @dst matrix
	template <typename Type>
	static inline void matrixPower(Type** src, int power, int row)
	{
		if(src==NULL)
		{
			std::cerr<<"matrixPower() is called for an empty source matrix!! Exiting..\n";
			exit(1);
		}

		Type** dst1 = new Type*[row];
		Type** dst2 = new Type*[row];
		for (int i=0;i<row;i++)
		{
			dst1[i] = new Type[row];
			dst2[i] = new Type[row];
		}
		
		copyMatrix(src,dst1,row,row);
		
		for(int i=2;i<=power;i++)
		{
			multiplyMatrices(dst1,src,dst2,row,row,row);
			copyMatrix(dst2,dst1,row,row);
		}
		
		copyMatrix(dst1,src,row,row);
		
		for(int i=0;i<row;i++)
		{
			delete[] dst1[i];
			delete[] dst2[i];
		}
		delete[] dst1;
		delete[] dst2;
		
		return;
	}
	
	//given a square matrix @src, this function performs matrix exponentiation (i.e., e^X = \sum_{k=0}^{@upperLimit} (1/k!)X^k) 
	//and writes the result into the @src matrix
	//Note that, this function does not add e^0 into the sum, that is, it does not add identity matrix to the sum.
	template <typename Type>
	static inline void matrixExponent(Type** src, int upperLimit, int row)
	{
		if(src==NULL)
		{
			std::cerr<<"matrixExponent() is called for an empty source matrix!! Exiting..\n";
			exit(1);
		}

		//allocate memory
		Type** power = new Type*[row];
		Type** temp = new Type*[row];
		Type** accumulator = new Type*[row];
		for (int i=0;i<row;i++)
		{
			power[i] = new Type[row];
			temp[i] = new Type[row];
			accumulator[i] = new Type[row];
		}
		
		copyMatrix(src,power,row,row);
		copyMatrix(src,accumulator,row,row);
				
		for(int i=2;i<=upperLimit;i++)
		{
			multiplyMatrices(power,src,temp,row,row,row);
			copyMatrix(temp,power,row,row);//power matrix keeps the src matrix raised to power i
			float nFactorialInv = 1.0/(float)Utility::factorial(i);
			multiplyMatrixWithScalar<float>(temp,nFactorialInv,row,row);//divide power matrix by factorial for addition
			elementwiseAddMatrices(accumulator,temp,accumulator,row,row);
		}
		
		copyMatrix(accumulator,src,row,row);
		
		//free memory
		for(int i=0;i<row;i++)
		{
			delete[] power[i];
			delete[] temp[i];
			delete[] accumulator[i];
		}
		delete[] power;
		delete[] temp;
		delete[] accumulator;
		
		return;
	}
	
	
	//scale the values inside the matrix to fit in [0,1] interval where oldMin maps to 0 and oldMax maps to 1
	static void normalizeMatrix(float** matrix, int row, int column, float rangeMin=0, float rangeMax=1)
	{
		float minValue=INF,maxValue=N_INF;
		for(int i=0;i<row;i++)
			for(int j=0;j<row;j++)
			{
				if(matrix[i][j]>maxValue)
					maxValue=matrix[i][j];
				if(matrix[i][j]<minValue)
					minValue=matrix[i][j];
			}
		if(maxValue==minValue)
		{
			std::cerr<<"Divide by zero encountered in normalizeMatrix()!!! Exiting..."<<std::endl;
			exit(1);
		}
		float scaleRatio=1.0 / (maxValue-minValue);
		for(int i=0;i<row;i++)
			for(int j=0;j<row;j++)
			{
				matrix[i][j] = (matrix[i][j]-minValue)*scaleRatio;
			}
	}
	
	//given a matrix and a scalar value, this function adds the @scalar value 
	//to each element of the @matrix and writes the result into the input @matrix
	template <typename Type>
	static inline void addMatrixWithScalar(Type** matrix, Type scalar, int row, int column)
	{
		if(matrix==NULL)
		{
			std::cerr<<"addMatrixWithScalar() is called for an empty source matrix!! Exiting..\n";
			exit(1);
		}

		for(int i=0;i<row;i++)
			for(int j=0;j<column;j++)
				matrix[i][j]=scalar+matrix[i][j];
		return;
	}
	
	//given a square matrix and a scalar value, this function adds the @scalar value 
	//to each element of the @matrix except the diagonal entries and writes the result into the input @matrix
	template <typename Type>
	static inline void addMatrixWithScalarExceptDiagonal(Type** matrix, Type scalar, int row)
	{
		if(matrix==NULL)
		{
			std::cerr<<"addMatrixWithScalar() is called for an empty source matrix!! Exiting..\n";
			exit(1);
		}

		for(int i=0;i<row;i++)
			for(int j=0;j<row;j++)
			{
				if(i==j)
					continue;
				matrix[i][j]=scalar+matrix[i][j];
			}
		return;
	}
	
	//given a matrix and a scalar value, this function subtracts each element of the @matrix
	//from the @scalar value and writes the result into the input @matrix
	template <typename Type>
	static inline void subtractMatrixFromScalar(Type** matrix, Type scalar, int row, int column)
	{
		if(matrix==NULL)
		{
			std::cerr<<"subtractMatrixFromScalar() is called for an empty source matrix!! Exiting..\n";
			exit(1);
		}

		for(int i=0;i<row;i++)
			for(int j=0;j<column;j++)
				matrix[i][j]=scalar-matrix[i][j];
		return;
	}
	
	//given a matrix and a scalar value, this function subtracts each element of the @matrix
	//from the @scalar value and writes the result into the input @matrix
	template <typename Type>
	static inline void subtractMatrixFromScalar(Type** src, Type** dst, Type scalar, int row, int column)
	{
		if(src==NULL | dst==NULL)
		{
			std::cerr<<"subtractMatrixFromScalar() is called for an empty source/destination matrix!! Exiting..\n";
			exit(1);
		}

		for(int i=0;i<row;i++)
			for(int j=0;j<column;j++)
				dst[i][j]=scalar-src[i][j];
		return;
	}
	
	//given a square matrix and a scalar value, this function subtracts each element of the @matrix
	//except the diagonal entries from the @scalar value and writes the result into the input @matrix
	template <typename Type>
	static inline void subtractMatrixFromScalarExceptDiagonal(Type** matrix, Type scalar, int row)
	{
		if(matrix==NULL)
		{
			std::cerr<<"subtractMatrixFromScalar() is called for an empty source matrix!! Exiting..\n";
			exit(1);
		}

		for(int i=0;i<row;i++)
			for(int j=0;j<row;j++)
			{
				if(i==j)
					continue;
				matrix[i][j]=scalar-matrix[i][j];
			}
		return;
	}
	
	//given a matrix to be masked and a mask matrix which is assumed to consist of 0,1 values,
	//this function multiplies the two matrices elementwise.
	template <typename Type>
	static inline void maskMatrix(Type** src, int** mask, int row, int column)
	{
		if(src==NULL || mask==NULL)
		{
			std::cerr<<"maskMatrix() is called for an empty source/destination matrix!! Exiting..\n";
			exit(1);
		}

		for(int i=0;i<row;i++)
			for(int j=0;j<column;j++)
				src[i][j]=src[i][j]*mask[i][j];
		return;
	}
	
	//given a @matrix and a @threshold, this function replaces the elements of the matrix that are less than
	//the threshold with the provided @value
   //@absoluteValue: if true, then filter out the values of the matrix whose absolute value is less than the given threshold
	template <typename Type>
	static inline void filterOutElementsOfMatrixLessThanThreshold(Type** matrix, Type threshold, Type value, int row, int column,bool absoluteValue=false)
	{
		if(matrix==NULL)
		{
			std::cerr<<"filterOutElementsOfMatrixLessThanThreshold() is called for an empty source/destination matrix!! Exiting..\n";
			exit(1);
		}

      if(absoluteValue==false)
      {
         for(int i=0;i<row;i++)
            for(int j=0;j<column;j++)
               matrix[i][j] = matrix[i][j] < threshold ? value : matrix[i][j];
      }
      else
      {
         for(int i=0;i<row;i++)
            for(int j=0;j<column;j++)
               matrix[i][j] = Utility::absoluteValue(matrix[i][j]) < threshold ? value : matrix[i][j];
      }
		return;
	}
   
   //given a @matrix and a @threshold, this function replaces the elements of the matrix that are less than
	//the threshold with the provided @value
   //@absoluteValue: if true, then filter out the values of the matrix whose absolute value is less than the given threshold
	template <typename Type>
	static inline void filterOutElementsOfMatrixLessThanOrEqualToThreshold(Type** matrix, Type threshold, Type value, int row, int column,bool absoluteValue=false)
	{
		if(matrix==NULL)
		{
			std::cerr<<"filterOutElementsOfMatrixLessThanThreshold() is called for an empty source/destination matrix!! Exiting..\n";
			exit(1);
		}

		if(absoluteValue==false)
		{
		   for(int i=0;i<row;i++)
			  for(int j=0;j<column;j++)
				 matrix[i][j] = matrix[i][j] <= threshold ? value : matrix[i][j];
		}
		else
		{
		   for(int i=0;i<row;i++)
			  for(int j=0;j<column;j++)
				 matrix[i][j] = Utility::absoluteValue(matrix[i][j]) <= threshold ? value : matrix[i][j];
		}
		return;
	}
	
	//given a @matrix and a @threshold, this function replaces the elements of the matrix that are less than
	//the threshold with the provided @value
   //@absoluteValue: if true, then filter out the values of the matrix whose absolute value is greater than the given threshold
	template <typename Type>
	static inline void filterOutElementsOfMatrixGreaterThanThreshold(Type** matrix, Type threshold, Type value, int row, int column,bool absoluteValue=false)
	{
		if(matrix==NULL)
		{
			std::cerr<<"filterOutElementsOfMatrixLessThanThreshold() is called for an empty source/destination matrix!! Exiting..\n";
			exit(1);
		}

		if(absoluteValue==false)
		{
		   for(int i=0;i<row;i++)
			  for(int j=0;j<column;j++)
				 matrix[i][j] = matrix[i][j] > threshold ? value : matrix[i][j];
		}
		else
		{
		   for(int i=0;i<row;i++)
			  for(int j=0;j<column;j++)
				 matrix[i][j] = Utility::absoluteValue(matrix[i][j]) > threshold ? value : matrix[i][j];
		}
         
		return;
	}
	
   //given a @matrix and a @threshold, this function replaces the elements of the matrix that are less than
	//the threshold with the provided @value
   //@absoluteValue: if true, then filter out the values of the matrix whose absolute value is greater than the given threshold
	template <typename Type>
	static inline void filterOutElementsOfMatrixGreaterThanOrEqualToThreshold(Type** matrix, Type threshold, Type value, int row, int column,bool absoluteValue=false)
	{
		if(matrix==NULL)
		{
			std::cerr<<"filterOutElementsOfMatrixLessThanThreshold() is called for an empty source/destination matrix!! Exiting..\n";
			exit(1);
		}

		if(absoluteValue==false)
		{
		   for(int i=0;i<row;i++)
			  for(int j=0;j<column;j++)
				 matrix[i][j] = matrix[i][j] >= threshold ? value : matrix[i][j];
		}
		else
		{
		   for(int i=0;i<row;i++)
			  for(int j=0;j<column;j++)
				 matrix[i][j] = Utility::absoluteValue(matrix[i][j]) >= threshold ? value : matrix[i][j];
		}
         
		return;
	}
   
	//given a @matrix and a @threshold, this function replaces the elements of the matrix that are less than
	//the threshold with 0, and the ones greater than or equal to the @threshold with 1
	template <typename Type>
	static inline void binarizeMatrix(Type** matrix, Type threshold, int row, int column)
	{
		if(matrix==NULL)
		{
			std::cerr<<"binarizeMatrix() is called for an empty source matrix!! Exiting..\n";
			exit(1);
		}

		for(int i=0;i<row;i++)
			for(int j=0;j<column;j++)
				matrix[i][j] = matrix[i][j] >= threshold ? 1 : 0;
	}
	
	//given a square matrix, produces a symmetric matrix by summing the symmetric entries and replacing them with their mean
	template <typename Type>
	static inline void symmetrizeMatrixByMean(Type** dst, int row)
	{
		if(dst==NULL)
		{
			std::cerr<<"symmetrizeMatrixByMean() is called for an empty destination matrix!! Exiting..\n";
			exit(1);
		}

		for(int i=0;i<row;i++)
			for(int j=i+1;j<row;j++)
			{
				dst[i][j] = (dst[i][j] + dst[j][i])/2;
				dst[j][i] = dst[i][j];
			}
	}
	
	//given an upper triangular square matrix, fills the lower triangular part with the mirror of upper triangular
	//at the and giving a symmetric matrix
	template <typename Type>
	static inline void mirrorUpperTriangleOfTheMatrixToLowerTriangle(Type** src, int row)
	{
		if(src==NULL)
		{
			std::cerr<<"mirrorUpperTriangleOfTheMatrixToLowerTriangle() is called for an empty matrix!! Exiting..\n";
			exit(1);
		}

		for(int i=0;i<row;i++)
			for(int j=i+1;j<row;j++)
				src[j][i] = src[i][j];
	}
	
	//given a square matrix and a line to be removed, this function removes the row and column corresponding to this lineNum, 
	//and saves the result in a separate matrix
	template <typename Type>
	static inline void removeRowAndColumnOfMatrix(Type** src,Type** dst, int lineNum, int size)
	{
		if(src==NULL || dst==NULL)
		{
			std::cerr<<"removeRowAndColumnOfMatrix() is called for an empty matrix!! Exiting..\n";
			exit(1);
		}

		for(int i=0;i<size;i++)
		{
			for(int j=0;j<size;j++)
			{
				if(i<lineNum && j<lineNum)
					dst[i][j]=src[i][j];
				else if(i<lineNum && j>lineNum)
					dst[i][j-1]=src[i][j];
				else if(i>lineNum && j<lineNum)
					dst[i-1][j]=src[i][j];
				else if(i>lineNum && j>lineNum)
					dst[i-1][j-1]=src[i][j];
			}
		}
	}
	
	template <typename Type>
	static inline void saveMatrixToFile(std::string filename, Type **matrix, int row, int col, char separator='\t', char endOfLine='\n', int precision=2)
	{
		std::ofstream file;
		file.open(filename.c_str());
		file.setf(std::ios::fixed,std::ios::floatfield);
		file.precision(precision);
   
		for(int i=0;i<row;i++)
		{
			for(int j=0;j<col;j++)
            if(matrix[i][j]<EPSILON && matrix[i][j]>-EPSILON)
               file<<0<<separator;
            else
               file<<matrix[i][j]<<separator;
			file<<endOfLine;
		}
		
		file.close();
	}
	
	template <typename Type>
	static inline void appendMatrixToFile(std::string filename, Type **matrix, int row, int col, char separator='\t', char endOfLine='\n', int precision=2)
	{
		std::ofstream file;
		file.open(filename.c_str(),std::ios_base::app);
		file.setf(std::ios::fixed,std::ios::floatfield);
		file.precision(precision);
   
		for(int i=0;i<row;i++)
		{
			for(int j=0;j<col;j++)
				file<<matrix[i][j]<<separator;
			file<<endOfLine;
		}
		
		file.close();
	}
	
	template <typename Type>
	static inline void loadMatrixFromFile(std::string filename, Type **matrix, int row, int column)
	{
		std::ifstream file; 
		std::istringstream is; 
		std::string line;

		////load names of result files
		file.open(filename.c_str());
		if(file.fail())
		{
			std::cout << "File not found in loadMatrixFromFile():"<<filename<<std::endl;
			return;
		}	

		for(int i=0;i<row;i++)
		{
			getline(file, line);
			is.clear();
			is.str(line);

			for(int j=0;j<column;j++)
				is >> matrix[i][j];
		}
		file.close();
	}
	
	template <typename Type>
	static inline void loadMatrixFromFileAfterSkipLine(std::string filename, Type **matrix, int row, int column, int skipLine)
	{
		std::ifstream file; 
		std::istringstream is; 
		std::string line;

		////load names of result files
		file.open(filename.c_str());
		if(file.fail())
		{
			std::cout << "File not found!!!!\n";
			return;
		}	
		for(int i=0;i<skipLine;i++)
			getline(file,line);
		
		for(int i=0;i<row;i++)
		{
			getline(file, line);
			is.clear();
			is.str(line);
			for(int j=0;j<column;j++)
				is >> matrix[i][j];
		}
		file.close();
	}
	
	template <typename Type>
	static inline void printMatrix(Type **matrix, int size, std::ostream &out=std::cout, char separator='\t', char endOfLine='\n')
	{
		for(int i=0;i<size;i++)
		{
			for(int j=0;j<size;j++)
				out<<matrix[i][j]<<separator;
			out<<endOfLine;
		}
	}
	
	template <typename Type>
	static inline void printMatrix(Type **matrix, int row, int col, std::ostream &out=std::cout, char separator='\t', char endOfLine='\n')
	{
		for(int i=0;i<row;i++)
		{
			for(int j=0;j<col;j++)
				out<<matrix[i][j]<<separator;
			out<<endOfLine;
		}
	}
	
	template <typename Type>
	static inline Type getMaxValueOfRowOfMatrix(Type **matrix, int row, int columnCount)
	{
		Type max = matrix[0][0];
		for(int i=0;i<columnCount;i++)
		{
			if(max < matrix[row][i])
				max = matrix[row][i];
		}
		return max;
	}
	
	template <typename Type>
	static inline Type getMinValueOfRowOfMatrix(Type **matrix, int row, int columnCount)
	{
		Type min = matrix[0][0];
		for(int i=0;i<columnCount;i++)
		{
			if(min > matrix[row][i])
				min = matrix[row][i];
		}
		return min;
	}
	
	template <typename Type>
	static inline Type getMaxValueOfColumnOfMatrix(Type **matrix, int column, int rowCount)
	{
		Type max = matrix[0][0];
		for(int i=0;i<rowCount;i++)
		{
			if(max < matrix[i][column])
				max = matrix[i][column];
		}
		return max;
	}
	
	template <typename Type>
	static inline Type getMinValueOfColumnOfMatrix(Type **matrix, int column, int rowCount)
	{
		Type min = matrix[0][0];
		for(int i=0;i<rowCount;i++)
		{
			if(min > matrix[i][column])
				min = matrix[i][column];
		}
		return min;
	}
	
	template <typename Type>
	static inline Type getMaxValueOfMatrix(Type **matrix, int row, int column)
	{
		Type max = matrix[0][0];
		for(int i=0;i<row;i++)
		{
			for(int j=0;j<column;j++)
			{
				if(max < matrix[i][j])
					max = matrix[i][j];
			}
		}
		return max;
	}
	
	template <typename Type>
	static inline Type getMaxValueOfMatrixExceptDiagonal(Type **matrix, int row)
	{
		Type max = matrix[0][1];
		for(int i=0;i<row;i++)
		{
			for(int j=0;j<row;j++)
			{
				if(i==j)
					continue;
				if(max < matrix[i][j])
					max = matrix[i][j];
			}
		}
		return max;
	}
	
	template <typename Type>
	static inline Type getMaxValueOfMatrix(Type **matrix, int row, int column, int &maxRow, int &maxColumn)
	{
		Type max = matrix[0][0];
		maxRow = 0;
		maxColumn = 0;
		for(int i=0;i<row;i++)
		{
			for(int j=0;j<column;j++)
			{
				if(max < matrix[i][j])
				{
					max = matrix[i][j];
					maxRow = i;
					maxColumn = j;
				}
			}
		}
		return max;
	}
	
	template <typename Type>
	static inline Type getMinValueOfMatrix(Type **matrix, int row, int column)
	{
		Type min = matrix[0][0];
		for(int i=0;i<row;i++)
		{
			for(int j=0;j<column;j++)
			{
				if(min > matrix[i][j])
					min = matrix[i][j];
			}
		}
		return min;
	}
	
	template <typename Type>
	static inline Type getMinValueOfMatrix(Type **matrix, int row, int column, int &minRow, int &minColumn)
	{
		Type min = matrix[0][0];
		minRow = 0;
		minColumn = 0;
		for(int i=0;i<row;i++)
		{
			for(int j=0;j<column;j++)
			{
				if(min > matrix[i][j])
				{
					min = matrix[i][j];
					minRow = i;
					minColumn = j;
				}
			}
		}
		return min;
	}

	template <typename Type>
	static inline Type getNonZeroMinValueOfMatrix(Type **matrix, int row, int column)
	{
		int i=0,j=0;
		Type min = INF;
		
		for(int i=0;i<row;i++)
		{
			for(int j=0;j<column;j++)
			{
				if(min > matrix[i][j] && matrix[i][j]>=EPSILON)
					min = matrix[i][j];
			}
		}
		return min;
	}
	
	template <typename Type>
	static inline Type getMinValueOfMatrixGreaterThanThreshold(Type **matrix, float threshold, int row, int column)
	{
		int i=0,j=0;
		Type min = INF;
		
		for(int i=0;i<row;i++)
		{
			for(int j=0;j<column;j++)
			{
				if(min > matrix[i][j] && matrix[i][j]>=threshold)
					min = matrix[i][j];
			}
		}
		return min;
	}
	
	//</editor-fold>
	
	static long getTime()
	{
		long timeMilliseconds;
		struct timeval  time_data; /* seconds since 0 GMT */

		gettimeofday(&time_data,NULL);

		timeMilliseconds  = time_data.tv_usec;
		timeMilliseconds /= 1000;
		timeMilliseconds += time_data.tv_sec * 1000 ;

		return timeMilliseconds;
	}
};


#endif
