#pragma once
#include <iostream>
#include <vector>

using namespace std;

double f(double x, double l, double Volleft, double Volright)
{
	if (x < l)
	{
		return Volleft;
	}
	else
	{
		return Volright;
	}
}

double Analitich_P(double x, double l, double Volleft, double Volright)
{
	double a = 1. / (l - 1 - Volright / Volleft * l);
	if (x < l)
	{
		return Volleft / Volright * a * x + 1;
	}
	else
	{
		return a * (x - 1);
	}
}

void vvodNachDan(double& P0, double& PN, double& h, int N)
{
	P0 = 1;
	PN = 0;
	h = 1.0 / N;
}

void vvodDan(double& L, double& k, double& m, double& p, double& mu)
{
	L = 100;
}


int progonka_3d(const vector<double>& c, const vector<double>& b, const vector<double>& a, const vector<double>& d, vector <double>& X, int n)
{
	//																			расчёт производится по принципуAX = D
	//========================================================================================================================================================================================================================   
	vector <double> Alpha(n, 0);
	vector <double> Betta(n, 0);
	//==========================================================================ПРЯМАЯ ПРОГОНКА
	Alpha[0] = -c[0] * 1. / b[0];
	Betta[0] = d[0] * 1. / b[0];
	for (int i = 1; i < n; i++)
	{
		double y = b[i] + a[i] * Alpha[i - 1];
		double t = d[i] - a[i] * Betta[i - 1];
		Alpha[i] = -c[i] * 1. / y;
		Betta[i] = t / y;
	}
	//==========================================================================ОБРАТНАЯ ПРОГОНКА				
	X[n - 1] = Betta[n - 1];

	for (int i = n - 2; i >= 0; i--)
	{
		X[i] = Alpha[i] * X[i + 1] + Betta[i];
	}
	return 0;
}


