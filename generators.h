#pragma once
#include <vector>
#include <ctime>
using namespace std;

class UniformSimulator
{
public:
	UniformSimulator();
	vector<double> simulate(int n, int seed = time(NULL));
};

class DistributionSimulator //uses uniform/Inverse transformation
{
public:
	DistributionSimulator() {};
	virtual vector<double> simulate(int n) {
		return vector<double>();
	};
protected:
	UniformSimulator uni;
};

class ExponentialSimulator :public DistributionSimulator //uses uniform/Inverse transformation
{
public:
	ExponentialSimulator();
	vector<double> simulate(int n, double lambda, int seed = time(NULL));
};

class NormalSimulator :public DistributionSimulator //uses uniform/Inverse transformation
{
public:
	NormalSimulator();
	vector<double> simulate(int n, double mean = 0, double stddev = 1, double rho = 0, int seed = time(NULL));
};

class BiNormalSimulator :public DistributionSimulator //uses uniform/Inverse transformation
{
public:
	BiNormalSimulator();
	vector<vector<double>> simulate(int n, double mean = 0, double stddev = 1, double rho = 0, int seed1 = time(NULL), int seed2=time(NULL));
private:
	NormalSimulator norm;
};