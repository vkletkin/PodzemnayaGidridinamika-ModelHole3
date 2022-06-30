#include "Function.h"

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
		return Volright / Volleft * a * x + 1;
	}
	else
	{
		return a * (x - 1);
	}
}

void vvodNach_bezrazmDan(double& P0, double& PN, double& h, int N)
{
	P0 = 1;
	PN = 0;
	h = 1.0 / N;
}

void VvodANDVivod_razmernDan(double& L, double& k0, double& m, double& delta_p, double& mu, const vector<double>& P, const vector<double>& X, double l,double Volleft,double Volright)
{
	vvod_razmernDan(L, k0, m, delta_p, mu);

	vector<double> U_NUM(X.size() - 1);
	for (int i = 0; i < X.size() -  1; i++)
	{
		U_NUM[i] = -(P[i + 1] - P[i]) / (X[i + 1] - X[i]) * f(1. / 2 * (X[i] + X[i + 1]), l, Volleft, Volright);
	}

	double U0 = k0 / mu * delta_p / L;
	double U = U0 * U_NUM[0];// тут может быть интеграл
	double V = U / m;
	double t = L / V;  // тут мы должны использовать интеграл сфотографированно



	cout<< "скорость           U0[м/с] = " << U0<<endl 
		<< "скорость фильтрации U[м/с] = " << U << endl 
		<< "скорость истинная   V[м/с] = " << V << endl 
		<< "время               t[сек] = " << t << endl 
		<< "время             t[сутки] = " << int(t / (60 * 60 * 24.)) << endl;
	cout << fixed;


	for (int i = 0; i < U_NUM.size(); i++)
	{
		cout << i << "  P" << " = " << P[i] << "  " << "  U = " << U_NUM[i] << "  " << " f = " << f(1. / 2 * (X[i] + X[i + 1]), l, Volleft, Volright) << endl;
	}

}

void vvod_razmernDan(double& L, double& k0, double& m, double& delta_p, double& mu)
{
	L = 100;
	k0 =1e-12;
	m = 0.2;
	delta_p = 1e6;
	mu = 1e-3;
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

void vivod_resaults(const vector<double>& x, const vector<double>& P, int N, double volleft, double volright, double l )
{
	
	cout << " l = " << l << "   f(x=0)= " << volleft << "   f(x=1)= " << volright << endl;
	cout << "x       " << " f(x)      " << "P a       " << "P        " << "U a       " << "U     " << endl;
	cout << fixed;

	for (int i = 0; i < N; i++)
	{
		cout << x[i] << " " <<f(x[i], l, volleft, volright)
			<< " " << Analitich_P(x[i], l, volleft, volright) << " " << P[i]
			<< " " << -Analitich_P(x[i], l, volleft, volright) * f(x[i], l, volleft, volright) << " " << -P[i] * f(x[i], l, volleft, volright)<<endl;
	}
}

void vvodN(int& N)
{
	cout << "N(uzli) = ";
	cin >> N;
}


double viviod_newyazka(const vector<double>& x, const vector<double>& P, int N, double volleft, double volright, double l)
{
	double S = 0;
	for (int i = 0; i < N; i++)
	{
		S += pow(Analitich_P(x[i], l, volleft, volright) - P[i],2);
		//cout << S << endl;
	}
	S = S / N ;
	S = sqrt(S);
	//cout << "Невязка(" << N - 1 << ") = " << S << endl;
	return S;
}

void Solve( vector<double>& B,  vector<double>& C, vector<double>& A,  vector<double>& D, vector<double>& P, vector<double>& x, int N, double volleft, double volright, double l)
{
	double P0, PN, h;
	vvodNach_bezrazmDan(P0, PN, h, N);

	for (int i = 0; i < N + 1; i++)
	{
		x[i] = i * h;
	}


	for (int i = 1; i < N; i++)
	{
		A[i] = 1 / (x[i] - x[i - 1]) * f(1. / 2 * (x[i - 1] + x[i]), l, volleft, volright);
		B[i] = 1 / (x[i + 1] - x[i]) * f(1. / 2 * (x[i] + x[i + 1]), l, volleft, volright);
		C[i] = -(A[i] + B[i]);
		D[i] = 0;
	}

	// первая строчка матрицы
	A[0] = 0; B[0] = 0; C[0] = 1; D[0] = 1; P[0] = P0;
	// последняя строчка матрицы
	A[N] = 0; B[N] = 0; C[N] = 1; D[N] = 0; P[N] = PN;

	progonka_3d(
		B,											// выше диагонали
		C,											// cама диагональ
		A,											// ниже диагонали
		D,											// свободные значение
		P,											// искомый вектор
		N + 1);

	//vivod_resaults(x, P, N + 1, volleft, volright, l);

}
