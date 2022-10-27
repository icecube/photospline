#ifndef PHOTOSPLINE_CORE_BIVARIATESPLINE_H
#define PHOTOSPLINE_CORE_BIVARIATESPLINE_H

#include <vector>
#include "Eigen/Dense"

namespace photospline{
	enum splinetype{
		bspline, 
		bivariatespline
	};

	class bspline
	{
	public:
		double eval(int const splineid, double const x) const;	
	};
	class bivariatespline
	{
	private:
		class vertex
		{
		public:
			double x,y;
			vertex(float x,float y):x(x),y(y){}
		};
		class triangle
		{
		public:
			vertex v1, v2, v3;
			std::vector<vertex> neighbours;
			Eigen::Vector3f dim1;
			Eigen::Vector3f dim2;
			Eigen::Vector3f normal;
			Eigen::Vector3f point;
		};
		std::vector<triangle> triangles;

		static double det(vertex v1,vertex v2,vertex v3);
		static std::tuple<vertex, vertex,vertex> makerighthanded(vertex v1,vertex v2,vertex v3);
		static vertex projecttoplane(Eigen::Vector3f x, triangle t);
		static std::tuple<double, double, double> barycentriccoordinates(vertex x, vertex v1, vertex v2, vertex v3);
		static float chi(vertex x, vertex v1, vertex v2, vertex v3);
		static double M1(vertex x, vertex v1, vertex v2, vertex v3);
		static double M2(vertex x, vertex v1, vertex v2, vertex n1, vertex n2);
		static double M3(vertex x, vertex v1, vertex v2, vertex v3, vertex n1, vertex n2);
	public:
		double eval(int const splineid, double const theta, double const phi) const;	
	};
} //namespace photospline

#endif /* PHOTOSPLINE_CORE_BIVARIATESPLINE_H */

