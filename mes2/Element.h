#pragma once
#include "Node.h"
#include "Input.h"
#include<array>
#include <iostream>
using namespace std;
class Element
{
	
public:
	int id;
	Element(int id, int * nodesId);
	Element();
	~Element();

	unsigned long int getId();

	int * getNodesId();
	
	int * nodesId;
	double conductionRatio;


	////////////////////////////
	int i, j;

	double k = Input::getConductivity();
	double c = Input::getSpecificHeat();
	double ro = Input::getDensity();
	double alfa = Input::getAlfa();
	double talfa = Input::getAmbientTemperature();

	double tabX[4];
	double tabY[4];

	double val = 1 / sqrt(3);

	double tabKSI[4] = { val*(-1), val, val, val*(-1) };
	double tabETA[4] = { val*(-1), val*(-1), val, val }; //punkty do calkowania po objetosci

	/*double tabKSI[4] = { -0.577350269189626, 0.577350269189626, 0.577350269189626, -0.577350269189626 };
	double tabETA[4] = { -0.577350269189626, -0.577350269189626, 0.577350269189626, 0.577350269189626 };*/

	double tabL[4];
	double tabPow[4];

	double N1pow1[4] = { 0.788675, 0.211325, 0, 0 }; //N z punktow do calkowania po powierzchni
	double N2pow1[4] = { 0.211325, 0.788675, 0 ,0 };

	double N1pow2[4] = { 0, 0.788675, 0.211325,  0 };
	double N2pow2[4] = { 0 , 0.211325, 0.788675, 0 };

	double N1pow3[4] = { 0, 0 , 0.788675, 0.211325 };
	double N2pow3[4] = { 0 ,0, 0.211325, 0.788675 };

	double N1pow4[4] = { 0.211325, 0, 0, 0.788675 };
	double N2pow4[4] = { 0.788675, 0 ,0,  0.211325 };

	double tabN1[4];
	double tabN2[4];
	double tabN3[4];
	double tabN4[4];
	double tabXp[4];
	double tabYp[4];
	double dNdksi[4][4];
	double dNdeta[4][4];
	double J[4][4];
	double detJ[4];
	double J1[4][4];
	double dNdx[4][4];
	double dNdy[4][4];

	double dNdxdNdxp1[4][4];
	double dNdxdNdxp2[4][4];
	double dNdxdNdxp3[4][4];
	double dNdxdNdxp4[4][4];

	double dNdydNdyp1[4][4];
	double dNdydNdyp2[4][4];
	double dNdydNdyp3[4][4];
	double dNdydNdyp4[4][4];


	double Kp1[4][4];
	double Kp2[4][4];
	double Kp3[4][4];
	double Kp4[4][4];

	double macierzH[4][4];


	double tabC1pc[4][4];
	double tabC2pc[4][4];
	double tabC3pc[4][4];
	double tabC4pc[4][4];

	double macierzC[4][4];


	double matrix_HBC[4][4];
	
	double matrixH_final[4][4];
	double tabJ[4];

	double pc1pow1[4][4];
	double pc2pow1[4][4];
	double sumpow1[4][4];


	double pc1pow2[4][4];
	double pc2pow2[4][4];
	double sumpow2[4][4];


	double pc1pow3[4][4];
	double pc2pow3[4][4];
	double sumpow3[4][4];


	double pc1pow4[4][4];
	double pc2pow4[4][4];
	double sumpow4[4][4];

	//wektor
	double wektorP[4];
	
	double wek1p1[4];
	double wek1p2[4];
	double sumawek1[4];
	double wek2p1[4];
	double wek2p2[4];
	double sumawek2[4];
	double wek3p1[4];
	double wek3p2[4];
	double sumawek3[4];
	double wek4p1[4];
	double wek4p2[4];
	double sumawek4[4];



	void setXY(Node **node);
	double * calculateVectorP();
	double ** calculateMatrixC();
	double ** calculateMatrixH();
	double ** calculateMatrixHBC();
	double ** matrix_final();

	
};
