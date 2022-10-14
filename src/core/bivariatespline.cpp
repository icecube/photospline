#include "photospline/bivariatespline.h"
#include <cmath>

double bivariatespline::det(vertex v1,vertex v2,vertex v3)
{
	return v2.x*v3.y - v3.x*v2.y - v1.x*v3.y + v3.x*v1.y + v1.x*v3.y - v3.y*v1.x;
}

std::tuple<int, int, int> bivariatespline::makerighthanded(vertex v1,vertex v2,vertex v3)
{
	if det(v1, v2, v3) < 0
	{
		return std::tuple(v1, v3, v2);
	}
	else
	{
		return std::tuple(v1, v2, v3);
	}
}

vertex projecttoplane(Vector3f x, triangle t)
{
	double prod1 = t.point.dot(t.normal);
	double prod2 = x.dot(t.normal);
	double prod3 = prod1 / prod2;
	Vector3f planeintersection = x * prod3;
	return vertex(t.dim1.dot(planeintersection), t.dim2.dot(planeintersection))
}

std::tuple<double, double, double> barycentriccoordinates(vertex x, vertex v1, vertex v2, vertex v3);
{
	double denominator = det(v1, v2, v3);
	return std::tuple(det(x, v2, v3)/denominator, det(v1, x, v3)/denominator, det(v1, v2, x)/denominator);
}

int chi(vertex x, vertex v1, vertex v2, vertex v3)
{
	auto x = barycentriccoordinates(x, v1, v2, v3);
	if (x[0] < 0 || x[1] < 0 || x[2] < 0)
	{
		return 0;
	}
	else
	{
		return 1;
	}
}

double M1(vertex x, vertex v1, vertex v2, vertex v3)
{
	if det(v1, v2, v3) < 0
	{
		return chi(x, v1, v2, v3)/det(v1, v3, v2);
	}
	else
	{
		return chi(x, v1, v2, v3)/det(v1, v2, v3);	
	}
}

double M2(vertex x, vertex v1, vertex v2, vertex n1, vertex n2)
{
	double result = 0;
	std::tuple<double, double, double> bcoords = barycentriccoordinates(x, v2, n1, n2);
	result += bcoords[0]*M1(x, v1, n1, n2);
	result += bcoords[1]*M1(x, v1, v2, n2);
	result += bcoords[2]*M1(x, v1, v2, n1);
	return result;
}

double M3(vertex x, vertex v1, vertex v2, vertex v3, vertex n1, vertex n2)
{
	double result = 0;
	std::tuple<double, double, double> bcoords = barycentriccoordinates(x, v1, v2, v3);
	result += bcoords[0]*M2(x, v2, v3, n1, n2);
	result += bcoords[1]*M2(x, v1, v3, n1, n2);
	result += bcoords[2]*M2(x, v1, v2, n1, n2);
	return result;
}


double bivariatespline::eval(int const splineid, double const theta, double const phi)
{
	Vector3f tempx(std::sin(theta)*std::cos(phi), std::sin(theta)*std::sin(phi), std::cos(theta));
	triangle t = trianles[splineid];
	vertex x = projecttoplane(tempx, t.normal, t.planepoint, t)
	double result = 0;
	int size = neighbours.size();
	for(i=0; i++; i<size)
	{
		n1 = neighbours[i];
		n2 =  neighbours[(i+1)%size]
		result += det(n1, n2, x)* M3(x, t.v1, t.v2, t.v3, n1, n2);
	}
	return result;

}

