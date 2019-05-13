#include "Input.h"

#include <string>

Input::Input()
{
	file.open("input.csv", std::ios::in);
	if (file.good() == true)
	{
		std::string line;

		getline(file, line, '\n');
		initialTemperature = stod(line);


		getline(file, line, '\n');
		simulationTime = stod(line);


		getline(file, line, '\n');
		simulationStepTime = stod(line);

		getline(file, line, '\n');
		ambientTemperature = stod(line);

		getline(file, line, '\n');
		alfa = stod(line);

		getline(file, line, '\n');
		H = stod(line);

		getline(file, line, '\n');
		B = stod(line);

		getline(file, line, '\n');
		n_h = stod(line);

		getline(file, line, '\n');
		n_b = stod(line);


		getline(file, line, '\n');
		specificHeat = stod(line);

		getline(file, line, '\n');
		conductivity = stod(line);


		getline(file, line, '\n');
		density = stod(line);

		getline(file, line, '\n');
		iterations = stod(line); //iteracje dla Gaussa

		file.close();
	}
}


Input::~Input()
{
}


double Input::getInitialTemperature()
{
	return initialTemperature;
}

double Input::getSimulationTime()
{
	return simulationTime;
}
double Input::getSimulationStepTime()
{
	return simulationStepTime;
}
double Input::getAmbientTemperature()
{
	return ambientTemperature;
}
double Input::getAlfa()
{
	return alfa;
}
double Input::getH()
{
	return H;
}
double Input::getB()
{
	return B;
}
double Input::getN_H()
{
	return n_h;
}
double Input::getN_B()
{
	return n_b;
}
double Input::getSpecificHeat()
{
	return specificHeat;
}
double Input::getConductivity()
{
	return conductivity;
}
double Input::getDensity()
{
	return density;
}

double Input::getIterationsNumber()
{
	return iterations;
}
double Input::initialTemperature;
double Input::simulationTime;
double Input::simulationStepTime;
double Input::ambientTemperature;
double Input::alfa;
double Input::H;
double Input::B;
double Input::n_h;
double Input::n_b;
double Input::specificHeat;
double Input::conductivity;
double Input::density;
double Input::iterations;