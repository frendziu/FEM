#include "Grid.h"
#include "Input.h"
#include <iostream>


Grid::Grid()
{
	double verticalNodeNumber = Input::getN_B();
	double horizontalNodeNumber = Input::getN_H();
	double nodeHight = Input::getB() / (verticalNodeNumber - 1);
	double nodeWidth = Input::getH() / (horizontalNodeNumber - 1);

	bool onBound;
	for (unsigned long int i = 0; i < horizontalNodeNumber; i++)
	{
		for (unsigned long int j = 0; j < verticalNodeNumber; j++)
		{
			//sprawdzamy czy wezel globalny jest po lewej czy po prawie stronie siatki
			if (i == 0 || i == horizontalNodeNumber - 1)
			{
				onBound = true;
			}
		//jesli nie jest ani po prawej ani po lewej to sprawdzamy czy jest u gory albo u dolu siatki
			else if (j == 0 || j == verticalNodeNumber - 1)
			{
				onBound = true;
			}
			else
			{
				onBound = false;
			}
			gridNodes.push_back(new Node(i*nodeWidth, j*nodeHight, onBound));
		}
	}
	//gridElements = new Element *[9];
	int id = 0;
	for (unsigned long int i = 0; i < horizontalNodeNumber - 1; i++)
	{
		for (unsigned long int j = 1; j < verticalNodeNumber; j++)
		{
			int * nodesId = new int[4];
			nodesId[0] = j + verticalNodeNumber * i;
			nodesId[1] = j + verticalNodeNumber  * (i + 1);
			nodesId[2] = j + 1 + verticalNodeNumber * (i + 1);
			nodesId[3] = j + 1 + verticalNodeNumber * i;

			gridElements.push_back(new Element(++id, nodesId));
		}
	}
	globalH = new double *[gridNodes.size()];
	for (int i = 0; i < gridNodes.size(); i++)
		globalH[i] = new double[gridNodes.size()];
	for (int i = 0; i < gridNodes.size(); i++)
		for (int j = 0; j < gridNodes.size(); j++)
			globalH[i][j] = 0;

	globalC = new double *[gridNodes.size()];
	for (int i = 0; i < gridNodes.size(); i++)
		globalC[i] = new double[gridNodes.size()];
	for (int i = 0; i < gridNodes.size(); i++)
		for (int j = 0; j < gridNodes.size(); j++)
			globalC[i][j] = 0;

	globalHBC = new double *[gridNodes.size()];
	for (int i = 0; i < gridNodes.size(); i++)
		globalHBC[i] = new double[gridNodes.size()];
	for (int i = 0; i < gridNodes.size(); i++)
		for (int j = 0; j < gridNodes.size(); j++)
			globalHBC[i][j] = 0;

	globalHHBC = new double *[gridNodes.size()];
	for (int i = 0; i < gridNodes.size(); i++)
		globalHHBC[i] = new double[gridNodes.size()];
	for (int i = 0; i < gridNodes.size(); i++)
		for (int j = 0; j < gridNodes.size(); j++)
			globalHHBC[i][j] = 0;

	globalP = new double[gridNodes.size()];
	for (int i = 0; i < gridNodes.size(); i++)
		globalP[i] = 0;

	globalP_final = new double[gridNodes.size()];
	for (int i = 0; i < gridNodes.size(); i++)
		globalP_final[i] = 0;
}


Grid::~Grid()
{

}


void Grid::displayGrid()
{
	for (unsigned long int i = 0; i < 9; i++)
	{
		std::cout << gridElements[i]->getId();
		int * nodesId = gridElements[i]->getNodesId();
		std::cout << " --> ";
		std::cout << nodesId[0] << ", ";
		std::cout << nodesId[1] << ", ";
		std::cout << nodesId[2] << ", ";
		std::cout << nodesId[3] << std::endl;
	}
}

vector<Node*> Grid::getNodes()
{
	return gridNodes;
}
//Element** Grid::getElements()
//{
//	return gridElements;
//}

vector<Element*> Grid::getElements()
{
	return gridElements;
}

void Grid::calculateGlobal(double *t_zero, double**AB)
{
	for (int i = 0; i < 16; i++)
		globalP_final[i] = 0;
	for (int i = 0; i < 16; i++)
		for (int j = 0; j < 16; j++)
		{
			globalP_final[i] += (globalC[i][j] / step_time)*t_zero[j];
		}
	for (int i = 0; i < 16; i++)
		globalP_final[i] += globalP[i];

	for (int i = 0; i < 16; i++)
	{
		for (int j = 0; j < 17; j++)
		{
			if (j == 16)
				AB[i][j] = globalP_final[i];
			else
				AB[i][j] = globalHHBC[i][j];
		}
	}
}

void Grid::Global(vector<Element*> elements, double size, double** AB, double* t_zero)
{
	for (int k = 0; k < size; k++)
	{
		globalP_final[k] = 0;
		for (int r = 0; r < size; r++)
		{
			globalH[k][r] = 0;
			globalC[k][r] = 0;
			globalHBC[k][r] = 0;
			globalHHBC[k][r] = 0;
		}
	}
	for (int i = 0; i < size; i++)
	{
		double **localH = elements[i]->calculateMatrixH();
		double **localC = elements[i]->calculateMatrixC();
		double **localHBC = elements[i]->calculateMatrixHBC();
		double *localP = elements[i]->calculateVectorP();
		int *nodeId = elements[i]->getNodesId();

		for (int j = 0; j < 4; j++)
		{
			globalP[nodeId[j] - 1] += localP[j];
			for (int k = 0; k < 4; k++)
			{
				globalH[nodeId[j] - 1][nodeId[k] - 1] += localH[j][k];
				globalC[nodeId[j] - 1][nodeId[k] - 1] += localC[j][k];
				globalHBC[nodeId[j] - 1][nodeId[k] - 1] += localHBC[j][k];
			}
		}
	}
	for (int i = 0; i < 16; i++)
		for (int j = 0; j < 16; j++)
		{
			globalHHBC[i][j] += globalH[i][j] + globalHBC[i][j] + (globalC[i][j] / step_time);
		}


}

void Grid::showGlobalMatrix(double size)
{
	cout << endl << "Global C" << endl;
	for (int i = 0; i < size; i++)
	{
	for (int j = 0; j < size; j++)
	cout << globalC[i][j] << " ";
	cout << endl;
	}

	cout << endl << endl << "Global H" << endl;
	for (int i = 0; i < size; i++)
	{
	for (int j = 0; j < size; j++)
	cout << globalH[i][j] << " ";
	cout << endl;
	}

	cout << endl << endl << "Global [H] + [C]/dT" << endl;
	for (int i = 0; i < size; i++)
	{
	for (int j = 0; j < size; j++)
	cout << globalHHBC[i][j] << " ";
	cout << endl;
	}


	cout << endl << endl << "Vector P Global" << endl;
	for (int i = 0; i < size; i++)
		cout << globalP_final[i] << "  ";
}