
#include "Source.cpp"

int main()
{
	setlocale(LC_ALL, "Russian");
	
	int N;
	cout << "N(uzli) = ";
	cin >> N;

	double P0, PN, h;
	vvodNachDan(P0, PN, h, N);

	vector<double> B(N + 1, 0);											// выше диагонали
	vector<double> C(N + 1, 0);											// cама диагональ
	vector<double> A(N + 1, 0);											// ниже диагонали
	vector<double> D(N + 1, 0);											// свободные значение
	vector<double> P(N + 1, 0);											// искомый вектор
	vector<double> x(N + 1, 0);


	for (int i = 0; i < N + 1; i++)
	{
		x[i] = i * h;
	}

	for (int i = 1; i < N; i++)
	{
		A[i] = 1 / (x[i] - x[i - 1]);
		B[i] = 1 / (x[i + 1] - x[i]);
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

	for (int i = 0; i < N + 1; i++)
	{
		cout << Analitich_P(x[i], 0.5, 1, 0) << "   " << P[i] << "  " << x[i] << "   " << endl;
	}

	/*
	double L, k, m, p, mu;
	vvodDan(L, k, m, p, mu);
	*/

	double U0 = 10e-12 * 10e6 / (10e-3 * 100.);
	double U = U0 * 1.;
	double V = U / 0.2;
	double t = 100. / V;
	cout << U0 << " " << U << " " << V << " " << t << "   " << t / (60 * 60 * 24.) << endl;
	cout << U;
}