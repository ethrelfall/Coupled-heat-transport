// Diffusion.h
// header for Diffusion.cpp

#pragma once

#include <iostream>

#define FINITEDIFFRES 200

#pragma warning(disable:4996)  //needed in VS2022 to stop security check error

class Diffusion
{
public:
	Diffusion(double DiffusionCoefficient, double timeStep);
	~Diffusion();
	void EvaluateTimestep(double LHS_val);
	void OutputData();

private:
	//model parameters
	double mD;  //diffusion coefficient

	//computational parameters
	int mResolution;
	double mDt;
	double mDx;

	//dimensionless group appearing in matrix used for time-evolution
	double mAlpha;

	double mTime;

	//coefficients used to solve tridiagonal matrix
	double* mpCdash;
	double* mpDdash;

	//storage for model state
	double* mpData;

	//I/O
	FILE * mpOutFile;
};
