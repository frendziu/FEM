#pragma once
#include "Node.h"
#include "Element.h"
#include <vector>

using namespace std;
class Grid
{

public:
	double **globalH;
	double **globalC;
	double **globalHBC;
	double **globalHHBC;
	double *globalP;
	double *globalP_final;
	double step_time = Input::getSimulationStepTime();
	Grid();
	~Grid();
	vector<Node*> gridNodes;
	vector<Element*> gridElements;
	//Element** gridElements;
	void displayGrid();
	vector<Node*> getNodes();
	//Element** getElements();
	vector<Element*> getElements();

	void Global(vector<Element*>, double, double**, double*);
	void showGlobalMatrix(double);
	void calculateGlobal(double*, double**);
};
