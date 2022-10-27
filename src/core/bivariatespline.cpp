#include "photospline/bivariatespline.h"
#include <cmath>

namespace photospline
{
	double bivariatespline::det(vertex v1,vertex v2,vertex v3)
	{
		return v2.x*v3.y - v3.x*v2.y - v1.x*v3.y + v3.x*v1.y + v1.x*v3.y - v3.y*v1.x;
	}

	std::tuple<bivariatespline::vertex, bivariatespline::vertex, bivariatespline::vertex> bivariatespline::makerighthanded(vertex v1,vertex v2,vertex v3)
	{
		if (det(v1, v2, v3) < 0)
		{
			return std::tuple<bivariatespline::vertex, bivariatespline::vertex, bivariatespline::vertex>(v1, v3, v2);
		}
		else
		{
			return std::tuple<bivariatespline::vertex, bivariatespline::vertex, bivariatespline::vertex>(v1, v2, v3);
		}
	}

	bivariatespline::vertex bivariatespline::projecttoplane(Eigen::Vector3f x, triangle t)
	{
		double prod1 = t.point.dot(t.normal);
		double prod2 = x.dot(t.normal);
		double prod3 = prod1 / prod2;
		Eigen::Vector3f planeintersection = x * prod3;
		return vertex(t.dim1.dot(planeintersection), t.dim2.dot(planeintersection));
	}

	std::tuple<double, double, double> bivariatespline::barycentriccoordinates(bivariatespline::vertex x, bivariatespline::vertex v1, bivariatespline::vertex v2, bivariatespline::vertex v3)
	{
		double denominator = det(v1, v2, v3);
		return std::tuple<double, double, double>(det(x, v2, v3)/denominator, det(v1, x, v3)/denominator, det(v1, v2, x)/denominator);
	}

	float bivariatespline::chi(bivariatespline::vertex x, bivariatespline::vertex v1, bivariatespline::vertex v2, bivariatespline::vertex v3)
	{
		double x1, x2, x3;
		std::tie(x1, x2, x3) = barycentriccoordinates(x, v1, v2, v3);
		if (x1 < 0 || x2 < 0 || x3 < 0)
		{
			return 0;
		}
		else if (x1 == 0 || x2 == 0 || x3 == 0)
		{
			return 0.5;
		}
		else
		{
			return 1;
		}
	}

	double bivariatespline::M1(bivariatespline::vertex x, bivariatespline::vertex v1, bivariatespline::vertex v2, bivariatespline::vertex v3)
	{
		if (det(v1, v2, v3) < 0)
		{
			return chi(x, v1, v2, v3)/det(v1, v3, v2);
		}
		else
		{
			return chi(x, v1, v2, v3)/det(v1, v2, v3);	
		}
	}

	double bivariatespline::M2(bivariatespline::vertex x, bivariatespline::vertex v1, bivariatespline::vertex v2, bivariatespline::vertex n1, bivariatespline::vertex n2)
	{
		double result = 0;
		double x1, x2, x3;
		std::tie(x1, x2, x3) = barycentriccoordinates(x, v2, n1, n2);
		result += x1*M1(x, v1, n1, n2);
		result += x2*M1(x, v1, v2, n2);
		result += x3*M1(x, v1, v2, n1);
		return result;
	}

	double bivariatespline::M3(bivariatespline::vertex x, bivariatespline::vertex v1, bivariatespline::vertex v2, bivariatespline::vertex v3, bivariatespline::vertex n1, bivariatespline::vertex n2)
	{
		double result = 0;
		double x1, x2, x3;
		std::tie(x1, x2, x3) = barycentriccoordinates(x, v1, v2, v3);
		result += x1*M2(x, v2, v3, n1, n2);
		result += x2*M2(x, v1, v3, n1, n2);
		result += x3*M2(x, v1, v2, n1, n2);
		return result;
	}


	double bivariatespline::eval(int const splineid, double const theta, double const phi) const
	{
		Eigen::Vector3f tempx(std::sin(theta)*std::cos(phi), std::sin(theta)*std::sin(phi), std::cos(theta));
		triangle t = triangles[splineid];
		vertex x = projecttoplane(tempx, t);
		double result = 0;
		int size = t.neighbours.size();
		for(int i=0; i<size; i++)
		{
			auto n1 = t.neighbours[i];
			auto n2 =  t.neighbours[(i+1)%size];
			result += det(n1, n2, x)* M3(x, t.v1, t.v2, t.v3, n1, n2);
		}
		return result;

	}

}