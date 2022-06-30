
#include "Function.h"

int main()
{
	setlocale(LC_ALL, "Russian");

	double l = 0.75, volleft = 1, volright = 0.025/0.475;
	//double l = 0.5, volleft = 1, volright = 0.1;

	
	int N;
	vvodN(N);

	vector<double> B(N + 1, 0);											// выше диагонали
	vector<double> C(N + 1, 0);											// cама диагональ
	vector<double> A(N + 1, 0);											// ниже диагонали
	vector<double> D(N + 1, 0);											// свободные значение
	vector<double> P(N + 1, 0);											// искомый вектор
	vector<double> x(N + 1, 0);

	Solve(
		B,											// выше диагонали
		C,											// cама диагональ
		A,											// ниже диагонали
		D,											// свободные значение
		P,											// искомый вектор
		x,											// РАЗБИЕНИЕ
		N,
		volleft,
		volright,
		l);
	/*
	ofstream File("1.txt");
	for (int N = 1; N <=10; N++)
	{
		Solve(
			B,											// выше диагонали
			C,											// cама диагональ
			A,											// ниже диагонали
			D,											// свободные значение
			P,											// искомый вектор
			x,											// РАЗБИЕНИЕ
			N,
			volleft,
			volright,
			l);

		File << N << " " << viviod_newyazka(x, P, N + 1, volleft, volright, l) << endl;
	}
	File.close();
	*/
	double L, k0, m, delta_p, mu;
	VvodANDVivod_razmernDan(L, k0, m, delta_p, mu, P, x, l, volleft, volright);
}

