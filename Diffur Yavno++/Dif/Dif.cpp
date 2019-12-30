// 3 лаба Диффуры Явный  Эйлер 
#include "pch.h"
#include <iostream>
#include <math.h>
#include <vector>
#include <iomanip>
#include <fstream>

using namespace std;

double w;
double du(double u1, double u2, double t, int i);
double dw(double t);

int main()

{
	ofstream fout("output.txt");
	ofstream foutX("outputX.txt");
	ofstream foutY1("outputY1.txt");
	ofstream foutY2("outputY2.txt");
	fout << setw(15) << "numOfIt" << setw(15) << "y1" << setw(15) << "y2" << setw(15) << "tk" << endl;
	double tmax = 0.01, tow = 0.001, t = pow(10, -9), T = 1;
	int n = 2, numOfIt = 0;
	vector <double> y(n);
	y[0] = 0;
	y[1] = -0.412;
	vector <double> e(n);
	e[0] = 0.001 * 0.412;//0,001
	e[1] = 0.001 * 0.412;
	vector <double> f(n);
	vector <double> tk(n);
	while (t < T)
	{
		numOfIt++;
		for (int i = 0; i < n; i++)
		{
			f[i] = du(y[0], y[1], t, i + 1);   
			tk[i] = e[i] / (fabs(f[i]) + (e[i] / tmax));  //условия выбора шага  
			y[i] = y[i] + tow * f[i];
			if (tk[i] < tow)
				tow = tk[i];
		}
		t = t + tow;
		fout << setw(15) << numOfIt;
		for (int i = 0; i < n; i++)
			fout << setw(15) << y[i];
		foutY1 << y[0] << endl;
		foutY2 << y[1] << endl;
		foutX << t << endl;
		fout << setw(15) << t << endl;
	}
	return 0;
}

double du(double u1, double u2, double t, int i)   // задание функций 
{
	double du;
	if (i == 1)
		du = -u1 * u2 + (sin(t) / t);
	if (i == 2)
	{
		w = dw(t);
		du = -pow(u2, 2) + ((2.5 + (w / 40)) * t / (1 + pow(t, 2)));
	}
	return du;
}

double dw(double t)
{

	w = 25 + t * (48 - 25);
	return w;
}