#include "pch.h"
#include <iostream>
#include <fstream>
#include <conio.h>
#include <iomanip>
#include <algorithm>
#include <cmath>

using namespace std;

double Newton(double *yk1, double *yk, double Tk, double tau, int n);

const double MAXConst = 1000000000;

void ImplicitEulerMethod(double *u0, int n)
{
	ofstream cout("ImplicitEulerMethod.txt");

	cout << "n\n";
	cout << fixed << setprecision(9);

	double T = 1;
	double eps = 0.001;
	double tauMin = 0.01;
	double tauMax = 0.01;

	double Tk = 0, Tk1 = 0;
	double * yk = new double[n];
	double * yk_1 = new double[n];
	double * yk1 = new double[n];
	for (int i = 0; i < n; i++)
	{
		yk[i] = u0[i];
		yk_1[i] = u0[i];
		yk1[i] = u0[i];
	}

	for (int i = 0; i < n; i++)
	{
		if (i == 0)
		{
			cout << "      t" << setw(20);
			cout << "      u" << i + 1 << setw(19);
		}
		else if (i == n - 1)
			cout << "      u" << i + 1 << endl;
	}

	double tauK = tauMin, tauK_1 = tauMin, tauK1;
	double * epsK = new double[n];
	int steps = 0;
	while (Tk < T)
	{
		Tk1 = Tk + tauK;

		Newton(yk1, yk, Tk, tauK, n);


		for (int k = 0; k < n; k++)
			epsK[k] = -1 * (tauK / (tauK + tauK_1)) * (yk1[k] - yk[k] - tauK * (yk[k] - yk_1[k]) / tauK_1);

		for (int k = 0; k < n; k++)    // опт 
		{
			if (fabs(epsK[k]) > eps)
			{
				tauK /= 2;
				Tk1 = Tk;
				for (int i = 0; i < n; i++)
					yk1[i] = yk[i];
				continue;
			}
		}

		tauK1 = MAXConst;
		for (int i = 0; i < n; i++)              // 3 зон
		{
			if (fabs(epsK[i]) > eps)
				tauK1 = min(tauK1, tauK / 2);
			if (eps / 4 < fabs(epsK[i]) && fabs(epsK[i]) <= eps)
				tauK1 = min(tauK1, tauK);
			if (fabs(epsK[i]) <= eps / 4)
				tauK1 = min(tauK1, 2 * tauK);
		}

		if (tauK1 > tauMax)
			tauK1 = tauMax;

		for (int i = 0; i < n; i++)
		{
			if (i == 0)
				cout << Tk << setw(20);
			cout << setw(20) << yk[i];
			if (i == n - 1)
				cout << endl;
		}

		for (int i = 0; i < n; i++)
		{
			yk_1[i] = yk[i];
			yk[i] = yk1[i];
		}
		tauK_1 = tauK;
		tauK = tauK1;
		Tk = Tk1;
		if (Tk > T)
			Tk = T;

		steps++;
	}
	cout << "\\steps = " << steps << endl;
	return;
}
