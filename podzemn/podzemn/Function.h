#pragma once
#include <iostream>
#include <vector>
#include <iomanip>
#include <math.h>
#include <fstream>

using namespace std;
double f(
	double x,
	double l, 
	double Volleft,
	double Volright
);

double Analitich_P(
	double x,
	double l,
	double Volleft,
	double Volright
);
void vvodNach_bezrazmDan(
	double& P0,
	double& PN,
	double& h, 
	int N
);
void VvodANDVivod_razmernDan(
	double& L,
	double& k0, 
	double& m,
	double& delta_p,
	double& mu,
	const vector<double>& P,
	const vector<double>& x,
	double l,
	double Volleft,
	double Volright
);

void vvod_razmernDan(
	double& L,
	double& k0,
	double& m,
	double& p,
	double& mu
);

int progonka_3d(
	const vector<double>& c,
	const vector<double>& b,
	const vector<double>& a,
	const vector<double>& d,
	vector <double>& X,
	int n
);

void vivod_resaults(
	const vector<double>& x,
	const vector<double>& P,
	int N, 
	double volleft,
	double volright,
	double l
);

void vvodN(
	int &N
);
double viviod_newyazka(
	const vector<double>& x,
	const vector<double>& P,
	int N,
	double volleft,
	double volright,
	double l
);

void Solve(
	vector<double>& B,
	vector<double>& C,
	vector<double>& A,
	vector<double>& D,
	vector <double>& P,
	vector<double>& x,
	int N,
	double volleft,
	double volright,
	double l
);