#pragma once
#include <fstream>

class Input
{
	std::fstream file;
	static double initialTemperature;
	static double simulationTime;
	static double simulationStepTime;
	static double ambientTemperature;
	static double alfa;
	static double H;
	static double B;
	static double n_h;
	static double n_b;
	static double specificHeat;
	static double conductivity;
	static double density;
	static double iterations;
public:
	Input();
	~Input();

	static double getInitialTemperature();
	static double getSimulationTime();
	static double getSimulationStepTime();
	static double getAmbientTemperature();
	static double getAlfa();
	static double getH();
	static double getB();
	static double getN_H();
	static double getN_B();
	static double getSpecificHeat();
	static double getConductivity();
	static double getDensity();
	static double getIterationsNumber();
};