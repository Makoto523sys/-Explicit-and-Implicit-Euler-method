// ConsoleApplication1.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include "pch.h"
#include <iostream>

void ImplicitEulerMethod(double *u0, int n);

int main()
{
	int n = 2;
	double * u = new double[n];
	u[0] = 0.0;   // начальные условия 
	u[1] = -0.412;

	ImplicitEulerMethod(u, n);
}
