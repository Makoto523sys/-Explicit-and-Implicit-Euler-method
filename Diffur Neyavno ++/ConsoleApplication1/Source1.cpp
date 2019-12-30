#include "pch.h"
#include <iostream>
#include <fstream>
#include <conio.h>
#include <iomanip>
#include <algorithm>
#include <cmath>

using namespace std;

const int itr_max = 100;

double** createMatrix(int x)
{
	double **A = new double*[x];
	for (int i = 0; i < x; i++)
		A[i] = new double[x + 1];
	return A;
}

void deleteMatrix(double **X, int x)
{
	for (int i = 0; i < x; i++)
		delete X[i];
	delete[] X;
}

void copyMatrix(double **X, double **copyX, int x)
{
	for (int i = 0; i < x; i++)
	{
		for (int j = 0; j < x + 1; j++)
			copyX[i][j] = X[i][j];
	}
}

double f1(double *uk1, double *uk, double t, double Tau)
{
	return uk1[0] - uk[0] - Tau * (-uk1[0] * uk1[1] + ((t < 1e-9) ? 0.0 : (sin(t) / t)));
}

double f2(double *uk1, double *uk, double t, double Tau)
{
	return uk1[1] - uk[1] - Tau * (-uk1[1] * uk1[1] + (3.125*t) / (1 + t * t));
}

bool Gauss(double *An, double **X, int x)
{
	for (int k = 0; k < x; k++)
	{
		double max = fabs(X[k][k]);
		int remeber = k;
		for (int i = k + 1; i < x; i++)
		{
			if (max < fabs(X[i][k]))
			{
				max = fabs(X[i][k]);
				remeber = i;
			}
		}

		if (fabs(max - 0) < 1e-6)
		{
			return 0;
		}

		if (k != remeber)
		{
			double *temp = X[k];
			X[k] = X[remeber];
			X[remeber] = temp;
		}

		double lead = X[k][k];
		for (int r = k; r < x + 1; r++)
		{
			X[k][r] /= lead;
		}

		for (int i = k + 1; i < x; i++)
		{
			double temp = X[i][k];
			for (int j = k; j < x + 1; j++)
			{
				X[i][j] -= X[k][j] * temp;
			}
		}
	}

	An[x - 1] = X[x - 1][x + 1 - 1];
	for (int i = x - 2; i >= 0; i--)
	{
		An[i] = X[i][x + 1 - 1];
		for (int j = i + 1; j < x + 1 - 1; j++)
		{
			An[i] -= X[i][j] * An[j];
		}
	}
	return 1;
}

typedef double(*pf)(double*, double*, double, double);
double Differential(pf f, double *uk1, double *uk, double t, double Tau, int n)
{
	double dx = 1e-9;
	double* D = new double[2];
	for (int i = 0; i < 2; i++)
	{
		if (i == n - 1)
		{
			D[i] = uk1[i] + dx;
			i++;
		}
		D[i] = uk1[i];
	}

	double F = f(uk1, uk, t, Tau);
	double dF = f(D, uk, t, Tau);

	delete[] D;
	return (dF - F) / dx;
}

double Newton(double *yk1, double *yk, double tk, double Tau, int n)
{
	double **F = createMatrix(n);
	double *An = new double[n];
	double e = 1e-9, b1 = 0, b2 = 0;

	int itr = 1;

	do {
		F[0][0] = Differential(f1, yk1, yk, tk, Tau, 1);
		F[0][1] = Differential(f1, yk1, yk, tk, Tau, 2);
		F[0][2] = -f1(yk1, yk, tk, Tau);
		F[1][0] = Differential(f2, yk1, yk, tk, Tau, 1);
		F[1][1] = Differential(f2, yk1, yk, tk, Tau, 2);
		F[1][2] = -f2(yk1, yk, tk, Tau);

		double **copyF = createMatrix(n);
		copyMatrix(F, copyF, n);

		if (!Gauss(An, copyF, n))
		{
			deleteMatrix(copyF, n);
			deleteMatrix(F, n);
			delete An;
			return 0;
		}

		yk1[0] += An[0];
		yk1[1] += An[1];

		if (fabs(f1(yk1, yk, tk, Tau)) > fabs(f2(yk1, yk, tk, Tau)))
			b1 = fabs(f1(yk1, yk, tk, Tau));
		else
			b1 = fabs(f2(yk1, yk, tk, Tau));

		for (int i = 0; i < n; i++)
		{
			if (fabs(An[i]) < 1)
				b2 = fabs(An[i]);
			else if (fabs(An[i]) >= 1)
				b2 = fabs(An[i] / yk1[i]);
		}

		deleteMatrix(copyF, n);
		itr++;
	} while ((b1 > e || b2 > e) && (itr < itr_max));

	deleteMatrix(F, n);
	delete An;

	return 1;
}