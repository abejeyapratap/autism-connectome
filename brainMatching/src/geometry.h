/* 
 * File:   geometry.h
 * Author: yusuf
 *
 * Generated on October 4, 2016, 5:44 PM
 */

#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <iostream>
#include <fstream>  //std::ofstream, std::ifstream
#include <sstream>
#include <cmath> //fabs()
#include <vector>
#include <limits> //std::numeric_limits<int>::max() and std::numeric_limits<double>::infinity();

//a three dimensional vector class
class Vector
{
    private:
		float x, y, z;

    public:
		Vector(){}
		
		Vector(float _x, float _y, float _z)
		{
			x = _x;
			y = _y;
			z = _z;
		}
		
		Vector(const Vector &other)
		{
			x = other.x;
			y = other.y;
			z = other.z;
		}
		
		Vector& operator=(const Vector& other)
		{
			x = other.x;
			y = other.y;
			z = other.z;
			return *this;
		}

		inline float getX(){return x;}
		
		inline float getY(){return y;}
		
		inline float getZ(){return z;}

		inline void set(float _x, float _y, float _z)
		{
			x = _x;
			y = _y;
			z = _z;
		}
		
		inline void set(Vector vec)
		{
			x = vec.x;
			y = vec.y;
			z = vec.z;
		}
		
		inline void randomlySetVector(int upperLimit)
		{
			srand((int)clock());
			x = rand()%upperLimit + (rand()%100)/100.0;
			y = rand()%upperLimit + (rand()%100)/100.0;
			z = rand()%upperLimit + (rand()%100)/100.0;
		}
		
		inline void add(float _x, float _y, float _z)
		{
			x += _x;
			y += _y;
			z += _z;
		}
		
		inline void add(Vector vec)
		{
			x += vec.x;
			y += vec.y;
			z += vec.z;
		}
		
		inline void addAndSet(Vector vec1, Vector vec2)
		{
			x = vec1.x + vec2.x;
			y = vec1.y + vec2.y;
			z = vec1.z + vec2.z;
		}
		
		inline void subtract(float _x, float _y, float _z)
		{
			x -= _x;
			y -= _y;
			z -= _z;
		}
		
		inline void subtract(Vector vec)
		{
			x -= vec.x;
			y -= vec.y;
			z -= vec.z;
		}
		
		inline Vector subtractAndSet(Vector vec1, Vector vec2)
		{
			x = vec1.x - vec2.x;
			y = vec1.y - vec2.y;
			z = vec1.z - vec2.z;

			return *this;
		}
		
		inline void multiply(float ratio)
		{
			x *= ratio;
			y *= ratio;
			z *= ratio;
		}
		
		inline void multiply(Vector ratio)
		{
			x *= ratio.x;
			y *= ratio.y;
			z *= ratio.z;
		}
		
		inline Vector multiplyAndReturn(float ratio)
		{
			Vector temp;
			temp.x = x * ratio;
			temp.y = y * ratio;
			temp.z = z * ratio;
			return temp;
		}
		
		inline void calculateVector(Vector P1, Vector P2)
		{//calculates vector that passes through P1,P2 points
			x = P2.x - P1.x;
			y = P2.y - P1.y;
			z = P2.z - P1.z;
		}
		
		inline void normalizeVector()
		{//calculates the unit vector
			float magnitude = this->calculateMagnitude();
			this->multiply(1.0/magnitude);
		}
		
		inline float calculateMagnitude()
		{
			return sqrt(x*x + y*y + z*z);
		}
		
		inline float calculateDistanceL1(Vector v1)
		{//for calculating distance btw two vertices
			return std::fabs(v1.x - x) + std::fabs(v1.y - y) + std::fabs(v1.z - z);
		}
		
		inline float calculateDistanceL2(Vector v1)
		{//for calculating distance btw two vertices
			return sqrt((v1.x - x)*(v1.x - x)+(v1.y - y)*(v1.y - y)+(v1.z - z)*(v1.z - z));
		}
		
		inline void averageVectors(Vector *vectorGroup, int count)
		{
			set(0,0,0);
			for(int i=0;i<count;i++)
			{
				x += vectorGroup[i].x;
				y += vectorGroup[i].y;
				z += vectorGroup[i].z;
			}
			this->multiply(1/(float)count);
		}
		
		inline void averageVectors(std::vector<Vector> &vectorGroup)
		{
			int count = vectorGroup.size();
			set(0,0,0);
			for(int i=0;i<count;i++)
			{
				x += vectorGroup[i].x;
				y += vectorGroup[i].y;
				z += vectorGroup[i].z;
			}
			this->multiply(1/(float)count);
		}
		
		inline void averageVectors(Vector vec1, Vector vec2)
		{
			x = (vec1.x + vec2.x) / 2;		
			y = (vec1.y + vec2.y) / 2;
			z = (vec1.z + vec2.z) / 2;
		}
		
		inline void calculateCrossProduct(Vector v1, Vector v2)
		{
			x = v1.y * v2.z - v1.z * v2.y;  
			y = v1.z * v2.x - v1.x * v2.z;
			z = v1.x * v2.y - v1.y * v2.x;
		}
		
		inline void calculateAngles(Vector v1)
		{//calculate angles of one vertex wrt origin on the 3D world coordinates
			x = std::atan2(v1.z,v1.y);
			y = std::atan2(v1.x,v1.z);
			z = std::atan2(v1.y,v1.x);
		}
		
		inline void calculateAngles(Vector v1, Vector v2)
		{//calculate angles of one vertex wrt another 
			x = std::atan2(v1.z-v2.z,v1.y-v2.y);
			y = std::atan2(v1.x-v2.x,v1.z-v2.z);
			z = std::atan2(v1.y-v2.y,v1.x-v2.x);
		}
		
		static inline float calculateAngle(Vector v1, Vector v2)
		{//given two vectors that lay on the same plane, find the angle between them
			//v1 . v2 = |v1||v2|cos(angle) => angle = acos((v1 . v2)/(|v1||v2|))
			float V1DotV2 = calculateDotProduct(v1,v2);
			float magnitudeV1 = v1.calculateMagnitude();
			float magnitudeV2 = v2.calculateMagnitude();
			float angle = std::acos(V1DotV2/(magnitudeV1*magnitudeV2));
			return angle;
		}
		
		static inline float calculateDotProduct(Vector v1, Vector v2)
		{
			return (v1.x * v2.x) + (v1.y * v2.y) + (v1.z * v2.z);
		}

		static inline float calculateDistance(Vector v1, Vector v2)
		{//for calculating distance btw two vertices
			return sqrt((v1.x - v2.x)*(v1.x - v2.x)+(v1.y - v2.y)*(v1.y - v2.y)+(v1.z - v2.z)*(v1.z - v2.z));
		}
		
		static inline float calculateMagnitude(Vector v1)
		{
			return sqrt(v1.x * v1.x + v1.y * v1.y + v1.z * v1.z);
		}
		
		static inline Vector calculateReflectionVector(Vector vector, Vector normal)
		{
			//R = V - ( 2 * (V [dot] N) * N )
			float dotProduct;
			Vector reflection;

			dotProduct = calculateDotProduct(vector, normal);//dotProduct = (v [dot] N)
			normal.multiply(2*dotProduct);
			reflection.subtractAndSet(vector, normal);

			return reflection;
		}

		inline std::string printToString(char separator)
		{
			std::ostringstream oss;
			oss<<x<<separator<<y<<separator<<z;
			return oss.str();
		}
		
		inline void print(std::ostream &out=std::cout, char separator='t',int precision=2)
		{
			out.precision(precision);
			out<<x<<separator<<y<<separator<<z<<separator;
		}

		/*static inline struct Vector calculateVectoralDistance(Vertex v1, Vertex v2)
		{
			Vector distance;
			distance.x = sqrt(((v1.y - v2.y)*(v1.y - v2.y))+((v1.z - v2.z)*(v1.z - v2.z)));
			distance.y = sqrt(((v1.x - v2.x)*(v1.x - v2.x))+((v1.z - v2.z)*(v1.z - v2.z)));
			distance.z = sqrt(((v1.x - v2.x)*(v1.x - v2.x))+((v1.y - v2.y)*(v1.y - v2.y)));
			return distance;
		}*/
};

class Triangle
{
	unsigned int corners[3];

	Triangle(){}

	Triangle(unsigned int a, unsigned int b, unsigned int c)
	{
		corners[0] = a;
		corners[1] = b;
		corners[2] = c;
	}
};

//a 3x3 matrix class
class Matrix
{
	public:
		Matrix(){};
		Matrix(Vector v1, Vector v2, Vector v3)
		{
			cell[0][0] = v1.getX();
			cell[0][1] = v1.getY();
			cell[0][2] = v1.getZ();
			cell[1][0] = v2.getX();
			cell[1][1] = v2.getY();
			cell[1][2] = v2.getZ();
			cell[2][0] = v3.getX();
			cell[2][1] = v3.getY();
			cell[2][2] = v3.getZ();			
		}
		~Matrix(){};

		inline void setRowwise(Vector v1, Vector v2, Vector v3)
		{
			cell[0][0] = v1.getX();
			cell[0][1] = v1.getY();
			cell[0][2] = v1.getZ();
			cell[1][0] = v2.getX();
			cell[1][1] = v2.getY();
			cell[1][2] = v2.getZ();
			cell[2][0] = v3.getX();
			cell[2][1] = v3.getY();
			cell[2][2] = v3.getZ();
		}
		inline void setColumnwise(Vector v1, Vector v2, Vector v3)
		{
			cell[0][0] = v1.getX();
			cell[1][0] = v1.getY();
			cell[2][0] = v1.getZ();
			cell[0][1] = v2.getX();
			cell[1][1] = v2.getY();
			cell[2][1] = v2.getZ();
			cell[0][2] = v3.getX();
			cell[1][2] = v3.getY();
			cell[2][2] = v3.getZ();
		}
		inline void set(Matrix mat)
		{
			for(int i=0;i<3;i++)
				for(int j=0;j<3;j++)
					cell[i][j] = mat.cell[i][j];
		}
		inline float calculateDeterminant()
		{
			return   ((cell[0][0]*cell[1][1]*cell[2][2]) + (cell[2][0]*cell[0][1]*cell[1][2])
				   + (cell[1][0]*cell[2][1]*cell[0][2]) - (cell[0][2]*cell[1][1]*cell[2][0]) 
				   - (cell[1][2]*cell[2][1]*cell[0][0]) - (cell[2][2]*cell[0][1]*cell[1][0]));
		}
		inline void multiplyMatricesAndSet(Matrix m1, Matrix m2)
		{//multiplies two 3x3 matrices and sets result to the matrix
			for(int i=0;i<3;i++)
				for(int j=0;j<3;j++)
				{
					cell[i][j] = m1.cell[i][0]*m2.cell[0][j] + m1.cell[i][1]*m2.cell[1][j] + m1.cell[i][2]*m2.cell[2][j];
				}
		}
		inline void multiplyWithVectorAndSetVector(Vector& position)
		{
			float x = position.getX()*cell[0][0] + position.getY()*cell[1][0] + position.getZ()*cell[2][0];
			float y = position.getX()*cell[0][1] + position.getY()*cell[1][1] + position.getZ()*cell[2][1];
			float z = position.getX()*cell[0][2] + position.getY()*cell[1][2] + position.getZ()*cell[2][2];
			position.set(x,y,z);
		}
		inline void calculateInverseMatrixAndSet(Matrix mat)
		{
			float determinant;

			determinant = mat.calculateDeterminant();

			cell[0][0] = (mat.cell[1][1]*mat.cell[2][2] - mat.cell[1][2]*mat.cell[2][1]) / determinant;
			cell[0][1] = (mat.cell[0][2]*mat.cell[2][1] - mat.cell[0][1]*mat.cell[2][2]) / determinant;
			cell[0][2] = (mat.cell[0][1]*mat.cell[1][2] - mat.cell[0][2]*mat.cell[1][1]) / determinant;
			cell[1][0] = (mat.cell[1][2]*mat.cell[2][0] - mat.cell[1][0]*mat.cell[2][2]) / determinant;
			cell[1][1] = (mat.cell[0][0]*mat.cell[2][2] - mat.cell[0][2]*mat.cell[2][0]) / determinant;
			cell[1][2] = (mat.cell[0][2]*mat.cell[1][0] - mat.cell[0][0]*mat.cell[1][2]) / determinant;
			cell[2][0] = (mat.cell[1][0]*mat.cell[2][1] - mat.cell[1][1]*mat.cell[2][0]) / determinant;
			cell[2][1] = (mat.cell[0][1]*mat.cell[2][0] - mat.cell[0][0]*mat.cell[2][1]) / determinant;
			cell[2][2] = (mat.cell[0][0]*mat.cell[1][1] - mat.cell[0][1]*mat.cell[1][0]) / determinant;
		}

	private:
		float cell[3][3];

};

class Geometry
{
	public:
		static constexpr float PI=3.1415926535897932384626433832795;
		static constexpr float RADIAN_TO_DEGREE=180.0/PI;
		static constexpr float DEGREE_TO_RADIAN=PI/180.0;
		static constexpr float INF=std::numeric_limits<float>::max();//10000000;
		static constexpr float N_INF=-INF;
		static const int NIL=-1;

		//given four points that form two roughly perpendicular lines to each other
		//this function calculates a local origin that is located at the center of four points.
		//outputs are written to the parameters: Vector& localOrigin 
		//inputs come with parameters:Vector upperPoint,Vector lowerPoint,Vector leftPoint,Vector rightPoint
		static void calculateLocalOrigin(Vector& localOrigin, Vector upperPoint,Vector lowerPoint,Vector leftPoint,Vector rightPoint)
		{
			float x = (upperPoint.getX() + lowerPoint.getX() + leftPoint.getX() + rightPoint.getX())/4;
			float y = (upperPoint.getY() + lowerPoint.getY() + leftPoint.getY() + rightPoint.getY())/4;
			float z = (upperPoint.getZ() + lowerPoint.getZ() + leftPoint.getZ() + rightPoint.getZ())/4;
			localOrigin.set(x,y,z);
		}

		//given four points that form two roughly perpendicular lines to each other
		//this function calculates a local coordinate system x,y,z but does not interfere with the origin
		//outputs are written to the parameters:Vector& xAxis,Vector& yAxis,Vector& zAxis
		//inputs come with parameters:Vector upperPoint,Vector lowerPoint,Vector leftPoint,Vector rightPoint
		static void calculateLocalAxis(Vector& xAxis,Vector& yAxis,Vector& zAxis, 
			Vector upperPoint,Vector lowerPoint,Vector leftPoint,Vector rightPoint)
		{
			xAxis.calculateVector(rightPoint, leftPoint);
			xAxis.normalizeVector();
			yAxis.calculateVector(lowerPoint, upperPoint);
			yAxis.normalizeVector();
			zAxis.calculateCrossProduct(xAxis, yAxis);
			zAxis.normalizeVector();
			xAxis.calculateCrossProduct(yAxis, zAxis);
			xAxis.normalizeVector();
		}

		//given a local coordinate system with three axis and the origin,
		//this function transforms the location of a given point to the local coordinates
		static void calculateLocalCoordinate(Vector& point, Vector xAxis,Vector yAxis,Vector zAxis, Vector localOrigin)
		{
			Matrix matrixA, matrixA_T;

			matrixA.setRowwise(xAxis,yAxis,zAxis);
			matrixA_T.calculateInverseMatrixAndSet(matrixA);

			Vector position;
			position.subtractAndSet(point,localOrigin);
			matrixA_T.multiplyWithVectorAndSetVector(position);
			point.set(position);
		}

		static inline float minimizeAngle(float val)
		{
			while(val >= 180.0)
				val -= 360.0;
			while(val <= -180.0)
				val += 360.0;
			return val;
		}
};

#endif /* GEOMETRY_H */

