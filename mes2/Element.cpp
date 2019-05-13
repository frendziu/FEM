#include "Element.h"
#include "Input.h"


Element::Element(int id, int * nodesId)
{
	this->id = id;
	this->nodesId = nodesId;
	//this->conductionRatio = Input::getConductivity();
	
}


Element::~Element()
{

}
Element::Element()
{

}

unsigned long int Element::getId()
{
	return id;
}

int * Element::getNodesId()
{
	return nodesId;
}


void Element::setXY(Node **node)
{
	tabX[0] = node[0]->x, tabX[1] = node[1]->x, tabX[2] = node[2]->x, tabX[3] = node[3]->x;
	tabY[0] = node[0]->y, tabY[1] = node[1]->y, tabY[2] = node[2]->y, tabY[3] = node[3]->y;
}

double ** Element::calculateMatrixH()
{

	for (i = 0; i < 4; i++)
	{
		tabN1[i] = 0.25 * (1 - tabKSI[i]) * (1 - tabETA[i]);
	}

	for (i = 0; i < 4; i++)
	{
		tabN2[i] = 0.25 * (1 + tabKSI[i]) * (1 - tabETA[i]);
	}

	for (i = 0; i < 4; i++)
	{
		tabN3[i] = 0.25 * (1 + tabKSI[i]) * (1 + tabETA[i]);
	}

	for (i = 0; i < 4; i++)
	{
		tabN4[i] = 0.25 * (1 - tabKSI[i]) * (1 + tabETA[i]);
	}

	for (i = 0; i < 4; i++)
	{
		tabXp[i] = (tabN1[i] * tabX[0]) + (tabN2[i] * tabX[1]) + (tabN3[i] * tabX[2]) + (tabN4[i] * tabX[3]);
	}

	for (i = 0; i < 4; i++)
	{
		tabYp[i] = (tabN1[i] * tabY[0]) + (tabN2[i] * tabY[1]) + (tabN3[i] * tabY[2]) + (tabN4[i] * tabY[3]);
	}

	for (i = 0; i < 4; i++)
	{
		dNdeta[i][0] = -0.25*(1 - tabKSI[i]);
		dNdeta[i][1] = -0.25*(1 + tabKSI[i]);
		dNdeta[i][2] = 0.25*(1 + tabKSI[i]);
		dNdeta[i][3] = 0.25*(1 - tabKSI[i]);
	}


	for (i = 0; i < 4; i++)
	{
		dNdksi[i][0] = -0.25*(1 - tabETA[i]);
		dNdksi[i][1] = 0.25*(1 - tabETA[i]);
		dNdksi[i][2] = 0.25*(1 + tabETA[i]);
		dNdksi[i][3] = -0.25*(1 + tabETA[i]);
	}

	for (j = 0; j < 4; j++) //j - punkt calkowania
	{
		J[j][0] = dNdksi[j][0] * tabX[0] + dNdksi[j][1] * tabX[1] + dNdksi[j][2] * tabX[2] + dNdksi[j][3] * tabX[3];
		J[j][1] = dNdksi[j][0] * tabY[0] + dNdksi[j][1] * tabY[1] + dNdksi[j][2] * tabY[2] + dNdksi[j][3] * tabY[3];
		J[j][2] = dNdeta[j][0] * tabX[0] + dNdeta[j][1] * tabX[1] + dNdeta[j][2] * tabX[2] + dNdeta[j][3] * tabX[3];
		J[j][3] = dNdeta[j][0] * tabY[0] + dNdeta[j][1] * tabY[1] + dNdeta[j][2] * tabY[2] + dNdeta[j][3] * tabY[3];
	}

	for (i = 0; i < 4; i++)
	{
		detJ[i] = J[i][0] * J[i][3] - J[i][1] * J[i][2]; //wyznacznik Jakobianu, s³u¿y do przechodzenia pomiêdzy uk³adami, a jego wyznacznacznik mowi jak duzy jest blad
	}


	for (j = 0; j < 4; j++)
	{
		J1[j][0] = J[j][0] / detJ[0];
		J1[j][1] = J[j][1] / detJ[1];
		J1[j][2] = J[j][2] / detJ[2];
		J1[j][3] = J[j][3] / detJ[3];
	}

	for (j = 0; j < 4; j++)
	{
		dNdx[j][0] = J1[j][0] * dNdksi[j][0] + J1[j][1] * dNdeta[j][0];
		dNdx[j][1] = J1[j][0] * dNdksi[j][1] + J1[j][1] * dNdeta[j][1];
		dNdx[j][2] = J1[j][0] * dNdksi[j][2] + J1[j][1] * dNdeta[j][2];
		dNdx[j][3] = J1[j][0] * dNdksi[j][3] + J1[j][1] * dNdeta[j][3];
	}

	double dNdy[4][4];

	for (j = 0; j < 4; j++)
	{
		dNdy[j][0] = J1[j][2] * dNdksi[j][0] + J1[j][3] * dNdeta[j][0];
		dNdy[j][1] = J1[j][2] * dNdksi[j][1] + J1[j][3] * dNdeta[j][1];
		dNdy[j][2] = J1[j][2] * dNdksi[j][2] + J1[j][3] * dNdeta[j][2];
		dNdy[j][3] = J1[j][2] * dNdksi[j][3] + J1[j][3] * dNdeta[j][3];
	}

	for (j = 0; j < 4; j++)
	{
		dNdxdNdxp1[j][0] = dNdx[0][0] * dNdx[0][j];
		dNdxdNdxp1[j][1] = dNdx[0][1] * dNdx[0][j];
		dNdxdNdxp1[j][2] = dNdx[0][2] * dNdx[0][j];
		dNdxdNdxp1[j][3] = dNdx[0][3] * dNdx[0][j];

	}

	for (j = 0; j < 4; j++)
	{
		dNdxdNdxp2[j][0] = dNdx[1][0] * dNdx[1][j];
		dNdxdNdxp2[j][1] = dNdx[1][1] * dNdx[1][j];
		dNdxdNdxp2[j][2] = dNdx[1][2] * dNdx[1][j];
		dNdxdNdxp2[j][3] = dNdx[1][3] * dNdx[1][j];

	}

	for (j = 0; j < 4; j++)
	{
		dNdxdNdxp3[j][0] = dNdx[2][0] * dNdx[2][j];
		dNdxdNdxp3[j][1] = dNdx[2][1] * dNdx[2][j];
		dNdxdNdxp3[j][2] = dNdx[2][2] * dNdx[2][j];
		dNdxdNdxp3[j][3] = dNdx[2][3] * dNdx[2][j];

	}

	for (j = 0; j < 4; j++)
	{
		dNdxdNdxp4[j][0] = dNdx[3][0] * dNdx[3][j];
		dNdxdNdxp4[j][1] = dNdx[3][1] * dNdx[3][j];
		dNdxdNdxp4[j][2] = dNdx[3][2] * dNdx[3][j];
		dNdxdNdxp4[j][3] = dNdx[3][3] * dNdx[3][j];

	}

	for (j = 0; j < 4; j++)
	{
		dNdydNdyp1[j][0] = dNdy[0][0] * dNdy[0][j];
		dNdydNdyp1[j][1] = dNdy[0][1] * dNdy[0][j];
		dNdydNdyp1[j][2] = dNdy[0][2] * dNdy[0][j];
		dNdydNdyp1[j][3] = dNdy[0][3] * dNdy[0][j];

	}

	for (j = 0; j < 4; j++)
	{
		dNdydNdyp2[j][0] = dNdy[1][0] * dNdy[1][j];
		dNdydNdyp2[j][1] = dNdy[1][1] * dNdy[1][j];
		dNdydNdyp2[j][2] = dNdy[1][2] * dNdy[1][j];
		dNdydNdyp2[j][3] = dNdy[1][3] * dNdy[1][j];

	}

	for (j = 0; j < 4; j++)
	{
		dNdydNdyp3[j][0] = dNdy[2][0] * dNdy[2][j];
		dNdydNdyp3[j][1] = dNdy[2][1] * dNdy[2][j];
		dNdydNdyp3[j][2] = dNdy[2][2] * dNdy[2][j];
		dNdydNdyp3[j][3] = dNdy[2][3] * dNdy[2][j];

	}

	for (j = 0; j < 4; j++)
	{
		dNdydNdyp4[j][0] = dNdy[3][0] * dNdy[3][j];
		dNdydNdyp4[j][1] = dNdy[3][1] * dNdy[3][j];
		dNdydNdyp4[j][2] = dNdy[3][2] * dNdy[3][j];
		dNdydNdyp4[j][3] = dNdy[3][3] * dNdy[3][j];

	}


	for (j = 0; j < 4; j++)
	{
		Kp1[j][0] = k*(dNdxdNdxp1[j][0] + dNdydNdyp1[j][0])* detJ[0];
		Kp1[j][1] = k*(dNdxdNdxp1[j][1] + dNdydNdyp1[j][1])* detJ[0];
		Kp1[j][2] = k*(dNdxdNdxp1[j][2] + dNdydNdyp1[j][2])* detJ[0];
		Kp1[j][3] = k*(dNdxdNdxp1[j][3] + dNdydNdyp1[j][3])* detJ[0];
	}

	for (j = 0; j < 4; j++)
	{
		Kp2[j][0] = k*(dNdxdNdxp2[j][0] + dNdydNdyp2[j][0])* detJ[1];
		Kp2[j][1] = k*(dNdxdNdxp2[j][1] + dNdydNdyp2[j][1])* detJ[1];
		Kp2[j][2] = k*(dNdxdNdxp2[j][2] + dNdydNdyp2[j][2])* detJ[1];
		Kp2[j][3] = k*(dNdxdNdxp2[j][3] + dNdydNdyp2[j][3])* detJ[1];
	}

	for (j = 0; j < 4; j++)
	{
		Kp3[j][0] = k*(dNdxdNdxp3[j][0] + dNdydNdyp3[j][0])* detJ[2];
		Kp3[j][1] = k*(dNdxdNdxp3[j][1] + dNdydNdyp3[j][1])* detJ[2];
		Kp3[j][2] = k*(dNdxdNdxp3[j][2] + dNdydNdyp3[j][2])* detJ[2];
		Kp3[j][3] = k*(dNdxdNdxp3[j][3] + dNdydNdyp3[j][3])* detJ[2];
		
	}

	for (j = 0; j < 4; j++)
	{
		Kp4[j][0] = k*(dNdxdNdxp4[j][0] + dNdydNdyp4[j][0])* detJ[3];
		Kp4[j][1] = k*(dNdxdNdxp4[j][1] + dNdydNdyp4[j][1])* detJ[3];
		Kp4[j][2] = k*(dNdxdNdxp4[j][2] + dNdydNdyp4[j][2])* detJ[3];
		Kp4[j][3] = k*(dNdxdNdxp4[j][3] + dNdydNdyp4[j][3])* detJ[3];

		

	}

	for (j = 0; j < 4; j++)
	{
		macierzH[j][0] = Kp1[j][0] + Kp2[j][0] + Kp3[j][0] + Kp4[j][0];
		macierzH[j][1] = Kp1[j][1] + Kp2[j][1] + Kp3[j][1] + Kp4[j][1];
		macierzH[j][2] = Kp1[j][2] + Kp2[j][2] + Kp3[j][2] + Kp4[j][2];
		macierzH[j][3] = Kp1[j][3] + Kp2[j][3] + Kp3[j][3] + Kp4[j][3];
	}

	cout << "Macierz H: " << endl;
	for (j = 0; j < 4; j++)
	{

		for (i = 0; i < 4; i++)
		{
			cout << macierzH[j][i] << "   ";
		}
		cout << endl;
	}
	double **tab3 = new double*[4];
	for (int i = 0; i < 4; i++)
		tab3[i] = macierzH[i];

	return tab3;
}

double ** Element::calculateMatrixC()
{
	for (int i = 0; i < 4; i++)
	{
		tabC1pc[i][0] = tabN1[i] * tabN1[0] * detJ[0] * c * ro;
		tabC1pc[i][1] = tabN1[i] * tabN1[1] * detJ[0] * c*ro;
		tabC1pc[i][2] = tabN1[i] * tabN1[2] * detJ[0] * c*ro;
		tabC1pc[i][3] = tabN1[i] * tabN1[3] * detJ[0] * c*ro;

	}

	for (int i = 0; i < 4; i++)
	{
		tabC2pc[i][0] = tabN2[i] * tabN2[0] * detJ[1] * c * ro;
		tabC2pc[i][1] = tabN2[i] * tabN2[1] * detJ[1] * c*ro;
		tabC2pc[i][2] = tabN2[i] * tabN2[2] * detJ[1] * c*ro;
		tabC2pc[i][3] = tabN2[i] * tabN2[3] * detJ[1] * c*ro;

	}

	for (int i = 0; i < 4; i++)
	{
		tabC3pc[i][0] = tabN3[i] * tabN3[0] * detJ[2] * c * ro;
		tabC3pc[i][1] = tabN3[i] * tabN3[1] * detJ[2] * c*ro;
		tabC3pc[i][2] = tabN3[i] * tabN3[2] * detJ[2] * c*ro;
		tabC3pc[i][3] = tabN3[i] * tabN3[3] * detJ[2] * c*ro;

	}

	for (int i = 0; i < 4; i++)
	{
		tabC4pc[i][0] = tabN4[i] * tabN4[0] * detJ[3] * c * ro;
		tabC4pc[i][1] = tabN4[i] * tabN4[1] * detJ[3] * c*ro;
		tabC4pc[i][2] = tabN4[i] * tabN4[2] * detJ[3] * c*ro;
		tabC4pc[i][3] = tabN4[i] * tabN4[3] * detJ[3] * c*ro;

	}

	for (j = 0; j < 4; j++)
	{
		macierzC[j][0] = tabC1pc[j][0] + tabC2pc[j][0] + tabC3pc[j][0] + tabC4pc[j][0];
		macierzC[j][1] = tabC1pc[j][1] + tabC2pc[j][1] + tabC3pc[j][1] + tabC4pc[j][1];
		macierzC[j][2] = tabC1pc[j][2] + tabC2pc[j][2] + tabC3pc[j][2] + tabC4pc[j][2];
		macierzC[j][3] = tabC1pc[j][3] + tabC2pc[j][3] + tabC3pc[j][3] + tabC4pc[j][3];
	}

	cout << "Matrix C: " << endl;
	for (j = 0; j < 4; j++)
	{

		for (i = 0; i < 4; i++)
		{
			cout << macierzC[j][i] << "   ";
		}
		cout << endl;
	}
	double **tab1 = new double*[4];
	for (int i = 0; i < 4; i++)
		tab1[i] = macierzC[i];

	return tab1;
}

double ** Element::calculateMatrixHBC()
{
	for (int i = 0; i <4; i++)
	{
		tabL[i] = sqrt(pow((tabX[(i + 1) % 4] - tabX[i]), 2) + pow((tabY[(i + 1) % 4] - tabY[i]), 2)); //dlugosc boku
	}
	for (i = 0; i < 4; i++)
	{
		tabJ[i] = tabL[i] / 2.0; //jakobian 1D
	}

	for (i = 0; i < 4; i++)
	{
		pc1pow1[i][0] = N1pow1[0] * N1pow1[i] * alfa;
		pc1pow1[i][1] = N1pow1[1] * N1pow1[i] * alfa;
		pc1pow1[i][2] = N1pow1[2] * N1pow1[i] * alfa;
		pc1pow1[i][3] = N1pow1[3] * N1pow1[i] * alfa;
	}


	for (i = 0; i < 4; i++)
	{
		pc2pow1[i][0] = N2pow1[0] * N2pow1[i] * alfa;
		pc2pow1[i][1] = N2pow1[1] * N2pow1[i] * alfa;
		pc2pow1[i][2] = N2pow1[2] * N2pow1[i] * alfa;
		pc2pow1[i][3] = N2pow1[3] * N2pow1[i] * alfa;
	}

	for (i = 0; i < 4; i++)
	{
		sumpow1[i][0] = (pc1pow1[i][0] + pc2pow1[i][0])*tabJ[0];
		sumpow1[i][1] = (pc1pow1[i][1] + pc2pow1[i][1])*tabJ[0];
		sumpow1[i][2] = (pc1pow1[i][2] + pc2pow1[i][2])*tabJ[0];
		sumpow1[i][3] = (pc1pow1[i][3] + pc2pow1[i][3])*tabJ[0];
	}

	for (i = 0; i < 4; i++)
	{
		pc1pow2[i][0] = N1pow2[0] * N1pow2[i] * alfa;
		pc1pow2[i][1] = N1pow2[1] * N1pow2[i] * alfa;
		pc1pow2[i][2] = N1pow2[2] * N1pow2[i] * alfa;
		pc1pow2[i][3] = N1pow2[3] * N1pow2[i] * alfa;
	}


	for (i = 0; i < 4; i++)
	{
		pc2pow2[i][0] = N2pow2[0] * N2pow2[i] * alfa;
		pc2pow2[i][1] = N2pow2[1] * N2pow2[i] * alfa;
		pc2pow2[i][2] = N2pow2[2] * N2pow2[i] * alfa;
		pc2pow2[i][3] = N2pow2[3] * N2pow2[i] * alfa;
	}

	for (i = 0; i < 4; i++)
	{
		sumpow2[i][0] = (pc1pow2[i][0] + pc2pow2[i][0])*tabJ[1];
		sumpow2[i][1] = (pc1pow2[i][1] + pc2pow2[i][1])*tabJ[1];
		sumpow2[i][2] = (pc1pow2[i][2] + pc2pow2[i][2])*tabJ[1];
		sumpow2[i][3] = (pc1pow2[i][3] + pc2pow2[i][3])*tabJ[1];
	}

	for (i = 0; i < 4; i++)
	{
		pc1pow3[i][0] = N1pow3[0] * N1pow3[i] * alfa;
		pc1pow3[i][1] = N1pow3[1] * N1pow3[i] * alfa;
		pc1pow3[i][2] = N1pow3[2] * N1pow3[i] * alfa;
		pc1pow3[i][3] = N1pow3[3] * N1pow3[i] * alfa;
	}


	for (i = 0; i < 4; i++)
	{
		pc2pow3[i][0] = N2pow3[0] * N2pow3[i] * alfa;
		pc2pow3[i][1] = N2pow3[1] * N2pow3[i] * alfa;
		pc2pow3[i][2] = N2pow3[2] * N2pow3[i] * alfa;
		pc2pow3[i][3] = N2pow3[3] * N2pow3[i] * alfa;
	}

	for (i = 0; i < 4; i++)
	{
		sumpow3[i][0] = (pc1pow3[i][0] + pc2pow3[i][0])*tabJ[2];
		sumpow3[i][1] = (pc1pow3[i][1] + pc2pow3[i][1])*tabJ[2];
		sumpow3[i][2] = (pc1pow3[i][2] + pc2pow3[i][2])*tabJ[2];
		sumpow3[i][3] = (pc1pow3[i][3] + pc2pow3[i][3])*tabJ[2];
	}
	for (i = 0; i < 4; i++)
	{
		pc1pow4[i][0] = N1pow4[0] * N1pow4[i] * alfa;
		pc1pow4[i][1] = N1pow4[1] * N1pow4[i] * alfa;
		pc1pow4[i][2] = N1pow4[2] * N1pow4[i] * alfa;
		pc1pow4[i][3] = N1pow4[3] * N1pow4[i] * alfa;
	}


	for (i = 0; i < 4; i++)
	{
		pc2pow4[i][0] = N2pow4[0] * N2pow4[i] * alfa;
		pc2pow4[i][1] = N2pow4[1] * N2pow4[i] * alfa;
		pc2pow4[i][2] = N2pow4[2] * N2pow4[i] * alfa;
		pc2pow4[i][3] = N2pow4[3] * N2pow4[i] * alfa;
	}

	for (i = 0; i < 4; i++)
	{
		sumpow4[i][0] = (pc1pow4[i][0] + pc2pow4[i][0])*tabJ[3];
		sumpow4[i][1] = (pc1pow4[i][1] + pc2pow4[i][1])*tabJ[3];
		sumpow4[i][2] = (pc1pow4[i][2] + pc2pow4[i][2])*tabJ[3];
		sumpow4[i][3] = (pc1pow4[i][3] + pc2pow4[i][3])*tabJ[3];
	}

	for (i = 0; i < 4; i++)
	{
		matrix_HBC[i][0] = tabPow[0] * sumpow1[i][0] + tabPow[1] * sumpow2[i][0] + tabPow[2] * sumpow3[i][0] + tabPow[3] * sumpow4[i][0];
		matrix_HBC[i][1] = tabPow[0] * sumpow1[i][1] + tabPow[1] * sumpow2[i][1] + tabPow[2] * sumpow3[i][1] + tabPow[3] * sumpow4[i][1];
		matrix_HBC[i][2] = tabPow[0] * sumpow1[i][2] + tabPow[1] * sumpow2[i][2] + tabPow[2] * sumpow3[i][2] + tabPow[3] * sumpow4[i][2];
		matrix_HBC[i][3] = tabPow[0] * sumpow1[i][3] + tabPow[1] * sumpow2[i][3] + tabPow[2] * sumpow3[i][3] + tabPow[3] * sumpow4[i][3];
	}

	cout << "MATRIX H - warunki brzegowe " << endl;
	for (j = 0; j < 4; j++)
	{

		for (i = 0; i < 4; i++)
		{
			cout << matrix_HBC[j][i] << "   ";
		}
		cout << endl;
	}

	double ** tab = new double*[4];
	for (int i = 0; i < 4; i++)
		tab[i] = matrix_HBC[i];
	return tab;

}

double ** Element::matrix_final()
{
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			matrixH_final[i][j] = (macierzH[i][j] + matrix_HBC[i][j]);

	cout << endl << "Matrix H final" << endl;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			cout << matrixH_final[i][j] << " ";
		cout << endl;
	}

	double ** tab = new double*[4];
	for (int i = 0; i < 4; i++)
		tab[i] = matrixH_final[i];
	return tab;
}

double *Element::calculateVectorP()
{
	for (i = 0; i < 4; i++)
	{
		wek1p1[i] = alfa * N1pow1[i] * talfa * tabJ[i];
	}

	for (i = 0; i < 4; i++)
	{
		wek1p2[i] = alfa * N2pow1[i] * talfa * tabJ[i];
	}

	for (i = 0; i < 4; i++)
	{
		sumawek1[i] = wek1p1[i] + wek1p2[i];
	}

	for (i = 0; i < 4; i++)
	{
		wek2p1[i] = alfa * N1pow2[i] * talfa * tabJ[i];
	}

	for (i = 0; i < 4; i++)
	{
		wek2p2[i] = alfa * N2pow2[i] * talfa * tabJ[i];
	}

	for (i = 0; i < 4; i++)
	{
		sumawek2[i] = wek2p1[i] + wek2p2[i];
	}

	for (i = 0; i < 4; i++)
	{
		wek3p1[i] = alfa * N1pow3[i] * talfa * tabJ[i];
	}

	for (i = 0; i < 4; i++)
	{
		wek3p2[i] = alfa * N2pow3[i] * talfa * tabJ[i];
	}

	for (i = 0; i < 4; i++)
	{
		sumawek3[i] = wek3p1[i] + wek3p2[i];
	}

	for (i = 0; i < 4; i++)
	{
		wek4p1[i] = alfa * N1pow4[i] * talfa * tabJ[i];
	}

	for (i = 0; i < 4; i++)
	{
		wek4p2[i] = alfa * N2pow4[i] * talfa * tabJ[i];
	}

	for (i = 0; i < 4; i++)
	{
		sumawek4[i] = wek4p1[i] + wek4p2[i];
	}




	for (i = 0; i < 4; i++)
	{
		wektorP[i] = tabPow[0] * sumawek1[i] + tabPow[1] * sumawek2[i] + tabPow[2] * sumawek3[i] + tabPow[3] * sumawek4[i];
	}

	cout << "Wektor P" << endl;
	for (i = 0; i < 4; i++)
	{
		cout << wektorP[i] << " ";

	}
	cout << endl;

	double *tab2 = new double[4];
	for (int i = 0; i < 4; i++)
		tab2[i] = wektorP[i];
	return tab2;
}


