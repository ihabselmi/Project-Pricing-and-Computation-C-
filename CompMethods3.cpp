// CompMethods3.cpp : Defines the entry point for the console application.
//
#include "stdafx.h"
#include <iostream>
#include <algorithm>
#include "generators.h"
#include <ctime>
#include <fstream>
#include <string>
#include <math.h>
# define M_PI           3.14159265358979323846

using namespace std;

#include "stdafx.h"

double vecavg(vector<double> rns)
{
	double sum = 0;
	for (double r : rns) { sum += r; }
	return sum / rns.size();
}

double vecvar(vector<double> rns)
{
	double x_bar = vecavg(rns);
	double sum_x_sq_centered = 0;
	for (double n : rns)
	{
		sum_x_sq_centered += pow(n - x_bar, 2);
	}
	return sum_x_sq_centered / (rns.size() - 1);
}


void q1(int seed, int n)
{
	//p(y2 > 5)
	vector<double> y2;
	double step_size = 0.01;
	NormalSimulator norm;
	vector<double> rns = norm.simulate(n* ((int)2.0 / step_size), 0, 1);
	for (int i = 0; i < n; i++)
	{
		double yt = 3.0 / 4;
		for (double j = 0; j < 2; j+= step_size)
		{
			double dyt = ((2 / (1 + j))*yt + (1 + pow(j, 3)) / 3) * step_size + ( (1 + pow(j, 3)) / 3) *pow(step_size, 0.5)*(rns.back()) ;
			rns.pop_back();
			yt += dyt;
		}
		y2.push_back(yt);
	}
	double count = 0;
	for (double y : y2) { if (y > 5) count++; }
	cout << "Prob: " << count / y2.size() << endl;

	//X2 ^ (1/3);
	vector<double> x2;
	rns.clear();  rns = norm.simulate(n* ((int)2.0 / step_size), 0, 1);
	for (int i = 0; i < n; i++)
	{
		double xt = 1;
		for (double j = 0; j < 2; j += step_size)
		{
			double dxt = (0.2 - 0.5*xt) * step_size + (2.0/3) *pow(step_size, 0.5)*(rns.back());
			rns.pop_back();
			xt += dxt;
		}
		double xt_m = pow(abs(xt), 1.0 / 3); if (xt < 0) xt_m *= -1;
		x2.push_back(xt_m);
	}
	cout << "E1: " << vecavg(x2) << endl;

	//E y3
	vector<double> y3;
	rns.clear();  rns = norm.simulate(n* ((int)2.0 / step_size), 0, 1);
	for (int i = 0; i < n; i++)
	{
		double yt = 3.0/4;
		for (double j = 0; j < 2; j += step_size)
		{
			double dyt = ((2 / (1 + j))*yt + (1 + pow(j, 3)) / 3) * step_size + ((1 + pow(j, 3)) / 3) *pow(step_size, 0.5)*(rns.back());
			rns.pop_back();
			yt += dyt;
		}
		y3.push_back(yt);
	}
	cout << "E2: " << vecavg(y3) << endl;

	//E x2 * y2 * 1(X2 > 1)
	vector<double> e3;
	for (int i = 0; i < n; i++)
	{
		if (x2[i] > 1)
		{
			e3.push_back(x2[i] * y2[i]);
		}
		else {
			e3.push_back(0.0);
		}
	}
	cout << "E3: " << vecavg(e3) << endl;
}
void q2(int seed, int n) 
{
	double step_size = 0.01;
	NormalSimulator norm;
	vector<double> w = norm.simulate(n* ((int)3.0 / step_size), 0, 1, 0, seed);
	vector<double> z = norm.simulate(n* ((int)3.0 / step_size), 0, 1, 0, seed*2);
	//E (x+1)^(1/3)
	vector<double> x3; int k = 0;
	for (int i = 0; i < n; i++)
	{
		double xt = 1.0;
		for (double j = 0; j < 2; j += step_size)
		{
			double dxt = 0.25 * xt * step_size + (1.0/3) *xt* pow(step_size, 0.5)*(w[k]) - 0.75 * xt* pow(step_size, 0.5)*(z[k]);
			xt += dxt;
			k++;
		}
		x3.push_back(pow(1+xt, 1.0/3));
	}
	cout << "E1: "<< vecavg(x3) << endl;

	//E (1+y3)^(1/3)

	vector<double> y3;
	for (int i = 0; i < n; i++)
	{
		double yt = exp(-0.08 * 3 + (1.0 / 3)*(pow(3, 0.5))*(w[i]) + (3.0 / 4)*(pow(3, 0.5))*(z[i]));
		y3.push_back(pow(1+ yt, 1.0 / 3));
	}
	cout << "E2: " << vecavg(y3) << endl;
}

double q3_montecarlo(int seed, double S0, double T, double X, double r, double sig)
{
	int n = 100000;
	NormalSimulator norm;
	vector<double> rns = norm.simulate(n, 0, 1, 0, seed);
	vector<double> callprices_reduced;
	for (int i = 0; i < n; i++)
	{
		callprices_reduced.push_back(
			(max(S0 * exp(sig*pow(T,0.5)*rns[i] + (r-0.5*pow(sig,2)) * T) - X, 0.0)
			+ max(S0 * exp(sig*pow(T, 0.5)*-rns[i] + (r - 0.5*pow(sig, 2)) * T) - X, 0.0)) / 2
			);
	}
	cout << "Call price(Montecarlo Reduced Var): " << exp(-r*T)*vecavg(callprices_reduced) << endl;
	return exp(-r*T)*vecavg(callprices_reduced);
}
double Normal_CDF(const double & x)
{
	double k = 1.0 / (1.0 + 0.2316419 * x);
	double k_sum = k * (0.319381530 + k * (-0.356563782 + k * (1.781477937 + k * (-1.821255978 + 1.330274429 * k))));

	if (x >= 0.0)
	{
		return (1.0 - (1.0 / (pow(2 * M_PI, 0.5)))*exp(-0.5*x*x) * k_sum);
	}
	else
	{
		return 1.0 - Normal_CDF(-x);
	}
}
double q3_bs(double S0, double T, double X, double r, double sig)
{
	double d1 = (log(S0 / X) + (r + pow(sig, 2) / 2)*T) / (sig*pow(T, 0.5));
	double d2 = d1 - sig*pow(T, 0.5);
	double nd1 = Normal_CDF(d1); double nd2 = Normal_CDF(d2);
	cout << "Call price(BS): " << S0*nd1 - X*nd2*exp(-r * T) << endl;
	return S0*nd1 - X*nd2*exp(-r * T);
}
void q4(int seed1, int seed2, double T, double X)
{
	int n = 10000;
	double p = -0.6; double r = 0.03; double S0 = 100; double V0 = 0.05; double sig = 0.42; double alpha = 5.8; double beta = 0.0625;
	double step_size = 0.01;
	vector<vector<double>> w;
	BiNormalSimulator binorm; w = binorm.simulate(n*((int)T/step_size), 0, 1, p, seed1, seed2);
	int k = 0;
	vector<double> call_prices;
	for (int i = 0; i < n; i++)
	{
		double vt = V0;
		double st = S0;
		for (double j = 0; j < T; j += step_size)
		{
			double dst = r*st*step_size + pow(vt, 0.5)*st*pow(step_size, 0.5)*w[k][0];
			double dvt = alpha*(beta - vt)*step_size + sig*pow(vt, 0.5)*pow(step_size, 0.5)*w[k][1];
			st += dst;
			vt += dvt;
			if (vt < 0) { vt = 0; }
			k++;
		}
		call_prices.push_back(exp(-r*T)*max(st - X, 0.0));
	}
	cout << "C1: " << vecavg(call_prices) << endl;
}
vector<int> halton_factors(int x, int base)
{
	vector<int> factors;
	while (x > 0)
	{
		factors.push_back(x%base);
		x = (int)floor(float(x) / base);
	}
	//for (int n : factors) { cout << n; }; cout << endl;
	return factors;
}
vector<double> halton_gen(int base, int n)
{
	vector<double> halton_rns;
	for (int k = 1; k < (n + 1); k++)
	{
		vector<int> f = halton_factors(k, 2);
		double h_k = 0;
		for (int j = 0; j < f.size(); j++)
		{
			h_k += ((double)f[j]) / pow((double)base, j + 1);
		}
		halton_rns.push_back(h_k);
	}
	return halton_rns;
}
void q5(int seed, int n)
{

	vector<vector<double>> uni_lgm;
	UniformSimulator uni;
	//b
	vector<double> rns = halton_gen(2, 100);
	rns = halton_gen(7, 100);
	//c
	rns = halton_gen(2, 100);
	rns = halton_gen(4, 100);

	for (double d : rns) { cout << d << endl; }
}
int main()
{
	//q1(1000, 10000);
	//q2(1000, 10000);
	//q3_montecarlo(1000, 100, 2, 110, 0.05, 0.3);
	//q3_bs(100, 2, 110, 0.05, 0.3);
	//q4(1000, 2000, 2, 100);
	q5(1000, 1000);
	char c; cin >> c;
    return 0;
}

