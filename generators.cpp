#include "stdafx.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <random>
#include <ctime>
#include <fstream>
#include "generators.h"
# define M_PI           3.14159265358979323846

using namespace std;

UniformSimulator::UniformSimulator() {
}
vector<double> UniformSimulator::simulate(int n, int seed)
{
	double a = 16807;
	double m = pow(2, 31) - 1;
	double q = floor(m / a); double r = int(m) % int(a);
	double x = seed;
	double h = 1 / double(m);
	vector<double> rns;
	rns.push_back(x / m);
	for (int i = 0; i < n; i++)
	{
		int k = x / q;
		x = a* (x - k*q) - k*r;
		if (x < 0) x = x + m;
		//x = fmod(x*a,int(m)) ;
		rns.push_back((x + 0.5)*h);
	}
	return rns;
}
ExponentialSimulator::ExponentialSimulator()
{

}
vector<double> ExponentialSimulator::simulate(int n, double lambda, int seed)
{
	vector<double> rns = uni.simulate(n, seed);
	vector<double> rnsexp;
	rnsexp.resize(rns.size());
	transform(rns.begin(), rns.end(), rnsexp.begin(), [&](double u) {return -1 / lambda*log(u); });
	return rnsexp;
}

NormalSimulator::NormalSimulator()
{

}
vector<double> NormalSimulator::simulate(int n, double mean, double stddev, double rho, int seed)
{
	vector<double> rns = uni.simulate(n, seed);
	vector<double> rnsnorm;
	for (int i = 0; i < rns.size() / 2; i++)
	{
		double u1 = rns[i * 2]; double u2 = rns[i * 2 + 1];
		double R = -2 * log(u1); double V = 2 * M_PI*u2;
		double z1 = pow(R, 0.5)*cos(V); double z2 = pow(R, 0.5)*sin(V);
		rnsnorm.push_back(z1*stddev + mean); rnsnorm.push_back(z2*stddev + mean);
	}
	return rnsnorm;
}

BiNormalSimulator::BiNormalSimulator()
{

}
vector<vector<double>> BiNormalSimulator::simulate(int n, double mean, double stddev, double rho, int seed1, int seed2)
{
	vector<double> rns1 = norm.simulate(n , mean, stddev, 0, seed1);
	vector<double> rns2 = norm.simulate(n, mean, stddev, 0, seed2);
	vector<vector<double>> bivrns;
	for (int i = 0; i < n; i++)
	{
		double z1 = rns1[i]; double z2 = rns2[i];
		double y1 = z1;
		double y2 = rho*z1 + pow(1 - pow(rho, 2), 0.5)*z2;
		vector<double> e; e.push_back(y1); e.push_back(y2); bivrns.push_back(e);
	}
	return bivrns;
}
