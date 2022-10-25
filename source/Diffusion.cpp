// Diffusion.cpp
// solves 1D time-dependent diffusion problem
// Neumann BC on left side, homogeneous (i.e. zero) Dirichlet on right side
// Finite-difference, second-order in time and space, Crank-Nicolson for time stepping (unconditionally stable)
// Ed Threlfall October 2022

#include "Diffusion.h"

Diffusion::Diffusion(double DiffusionCoefficient, double timeStep) : mD(DiffusionCoefficient), mDt(timeStep)
{
	mResolution = FINITEDIFFRES;
	mDx = 1.0 / mResolution;

    mAlpha = 0.5*mD * mDt / (mDx * mDx);  //Crank-Nicolson method gives the factor of 0.5

    mTime = 0.0;

	mpCdash = new double[mResolution];
	mpDdash = new double[mResolution];
	mpData  = new double[mResolution];

    for (int i=0; i<mResolution; i++)
    {
        mpCdash[i] = 0.0;
        mpDdash[i] = 0.0;
        mpData[i]  = 0.0;
    }

    mpOutFile = fopen("diffusion_heat_flux_vs_t.txt", "w");
    fprintf(mpOutFile, "t,L,R\n");
}

Diffusion::~Diffusion()
{
	delete[] mpCdash;
	delete[] mpDdash;
	delete[] mpData;

    fclose(mpOutFile);
}
void Diffusion::EvaluateTimestep(double LHS_val)
{
    //Gauss elimination to solve tridiagonal system
    //see e.g. The Nature of Mathematical Modelling, N. Gershenfeld, p.83
    //uses Crank-Nicolson to improve time accuracy to second-order

    //mpData[0] = LHS_val;  //for Dirichlet BC LHS; not needed if Neumann LHS 
    mpData[mResolution - 1] = 0.0;              //Dirichlet BC RHS (zero-temperature)

    //For Crank-Nicolson data needs to be transformed into data + alpha/2 * matrix * data
    double* pDataCopy = new double[mResolution];
#if 1
    for (int i=0; i<mResolution; i++)
    {
        pDataCopy[i] = mpData[i];
    }

    //Neumann BC LHS.  Was derived using central difference at point zero
    //and T_{1}-T_{-1} = dx * gradient to derive value for putative T_{-1}
    mpData[0] = pDataCopy[0] + 2.0*mAlpha * (-pDataCopy[0]+pDataCopy[1]-2.0*mDx*LHS_val);

    for (int i = 1; i < mResolution-1; i++)
    {
        mpData[i] = mpData[i] + mAlpha * (pDataCopy[i+1] - 2.0 * pDataCopy[i]+ pDataCopy[i-1]);
    }
#endif

    //for Neumann BC LHS
    mpCdash[0] = -2.0*mAlpha / (1.0 + 2.0*mAlpha);
    mpDdash[0] = mpData[0] / (1.0 + 2.0*mAlpha);

    //for Dirichlet BC LHS
    //mpCdash[0] = 0.0;
    //mpDdash[0] = mpData[0];

    for (int i = 0; i < mResolution - 2; i++)
    {
        double recipDenom = 1.0/(1.0 + 2.0 * mAlpha + mAlpha * mpCdash[i]);
        mpCdash[i + 1] = -mAlpha * recipDenom;
        mpDdash[i + 1] = (mpData[i + 1] + mAlpha * mpDdash[i]) * recipDenom;
    }

    mpCdash[mResolution - 1] = 0;
    mpDdash[mResolution - 1] = mpData[mResolution - 1];

    for (int i = 0; i < mResolution - 1; i++)
    {
        int idx = mResolution - i - 2;
        mpData[idx] = mpDdash[idx] - mpCdash[idx] * mpData[idx + 1];
    }

    mTime += mDt;

    delete[] pDataCopy;
}

void Diffusion::OutputData()
{
    //output gradients at LHS and RHS
    fprintf(mpOutFile, "%.6e, %.6e, %.6e\n", mTime, (mpData[0] - mpData[1]) / mDx, (mpData[mResolution - 2] - mpData[mResolution - 1]) / mDx);
}
