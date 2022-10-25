// Lorenz.h
// header for Lorenz.cpp

#pragma once

#define NUMSTAGES 5

class Lorenz
{
public:
	Lorenz(double rho, double sigma, double beta, double timeStep);
	~Lorenz();
	void EvaluateTimestep();
	void EvaluateTimeDerivs();
	double OutputNusselt();

private:
	//model parameters
	//all the following model parameters are >0
    //see https://www2.physics.ox.ac.uk/sites/default/files/profiles/read/lect6-43147.pdf
	double mRho;    //Rayleigh number
	double mSigma;  //Prandtl number
	double mBeta;   //coupling strength

	//initial data
	double mX0;
	double mY0;
	double mZ0;

	//computational parameters
	double mDt;  //for one-way coupling, does not really matter if timesteps not the same between models

	//model state storage appropriate to LSRK
	double mX;
	double mY;
	double mZ;
	double mXDot;
	double mYDot;
	double mZDot;
	double mXS;
	double mYS;
	double mZS;

};
