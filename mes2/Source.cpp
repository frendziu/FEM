#include "Input.h"
#include "Element.h"
#include "Grid.h"
#include "element.h"
#include <iostream>
#include<vector>

#include<ctime>
#include<cstdlib>
#include<fstream>
#include<string>

#include<cmath>
#include<iomanip>
using namespace std;

const double eps = 1e-12;
bool gauss(int n, double ** AB, double * X)
{
	int i, j, k;
	double m, s;

	// eliminacja wspó³czynników

	for (i = 0; i < n - 1; i++)
	{
		for (j = i + 1; j < n; j++)
		{
			if (fabs(AB[i][i]) < eps) return false;
			m = -AB[j][i] / AB[i][i];
			for (k = i + 1; k <= n; k++)
				AB[j][k] += m * AB[i][k];
		}
	}

	// wyliczanie niewiadomych

	for (i = n - 1; i >= 0; i--)
	{
		s = AB[i][n];
		for (j = n - 1; j >= i + 1; j--)
			s -= AB[i][j] * X[j];
		if (fabs(AB[i][i]) < eps) return false;
		X[i] = s / AB[i][i];
	}
	return true;
}


int main()
{

	Input input;

	Grid grid;
	grid.displayGrid();

	double nH = Input::getN_H();
	double	nL = Input::getN_B();

	int size = (nH - 1) * (nL - 1);
	int size2 = nH * nL;
	/*element l;
	l.calculateMatrixH();
	l.calculateMatrixC();
	l.calculateMatrixHBC();
	l.calculateVectorP();*/

	//Element** elements = new Element*[size];
	//elements = grid.getElements();
	vector<Element*> elements = grid.getElements();
	vector<Node*> nodes = grid.getNodes();
	Node **tempNodes = new Node*[4];

	for (int i = 0; i < 9; i++)
	{
		cout << "Iteracja " << i + 1 << endl;
		int *nodesId = elements[i]->getNodesId();
		for (int i = 0; i < 4; i++)
		{
			tempNodes[i] = nodes[nodesId[i] - 1];

		}

		elements[i]->setXY(tempNodes);
		elements[i]->calculateMatrixH();
		elements[i]->calculateMatrixC();


		if (tempNodes[0]->onBound && tempNodes[1]->onBound)
			elements[i]->tabPow[0] = 1;
		else
			elements[i]->tabPow[0] = 0;
		if (tempNodes[1]->onBound && tempNodes[2]->onBound)
			elements[i]->tabPow[1] = 1;
		else
			elements[i]->tabPow[1] = 0;
		if (tempNodes[2]->onBound && tempNodes[3]->onBound)
			elements[i]->tabPow[2] = 1;
		else
			elements[i]->tabPow[2] = 0;
		if (tempNodes[3]->onBound && tempNodes[0]->onBound)
			elements[i]->tabPow[3] = 1;
		else
			elements[i]->tabPow[3] = 0;

		elements[i]->calculateMatrixHBC();
		elements[i]->matrix_final();
		elements[i]->calculateVectorP();
		cout << endl << endl;


	}

	double **AB, *X;
	int      n, i, j;
	cout << setprecision(3) << fixed;

	// odczytujemy liczbe niewiadomych
	n = size2;

	// tworzymy macierze AB i X
	AB = new double *[n];
	X = new double[n];
	for (i = 0; i < n; i++) AB[i] = new double[n + 1];

	double t_zero[16];
	for (int i = 0; i < 16; i++)
	{
		t_zero[i] = Input::getInitialTemperature();
	
	}


	grid.Global(elements, size, AB, t_zero);
	for (int k = 0; k < 10; k++)
	{
		grid.calculateGlobal(t_zero, AB);
		grid.showGlobalMatrix(size2);
		cout << endl << "X" << endl;
		if (gauss(n, AB, X))
		{
			for (i = 0; i < n; i++)
				cout << "x" << i + 1 << " = " << setw(9) << X[i]
				<< endl;
		}
		else
			cout << "DZIELNIK ZERO\n";

		for (int r = 0; r < 16; r++)
			t_zero[r] = X[r];

	}




	for (i = 0; i < n; i++) delete[] AB[i];
	delete[] AB;
	delete[] X;


	system("pause");
	return 0;
}
