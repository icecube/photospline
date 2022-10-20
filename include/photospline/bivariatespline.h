#ifndef PHOTOSPLINE_CORE_BIVARIATESPLINE_H
#define PHOTOSPLINE_CORE_BIVARIATESPLINE_H

#include <vector>
#include "Eigen/Dense"

namespace photospline{
	enum splinetype{
		bspline, 
		bivariatespline
	}

	class bspline
	{
	public:
		double eval(int const splineid, double const x) const;	
	}
	class bivariatespline
	{
	private:

		class vertex
		{
		public:
			double x,y;
			vertex(x,y):x(x),y(y){}
		}
		class triangle
		{
			vertex v1, v2, v3;
			std::vector<vertex> neighbours;
			Vector3f dim1;
			Vector3f dim2;
			Vector3f normal;
			Vector3f point;
		}
		std::vector<triangle> triangles;

		double det(vertex v1,vertex v2,vertex v3);
	public:
		double eval(int const splineid, double const theta, double const phi) const;	
	}
} //namespace photospline

#endif /* PHOTOSPLINE_CORE_BIVARIATESPLINE_H */

