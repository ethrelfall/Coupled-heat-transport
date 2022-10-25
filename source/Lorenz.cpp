// Lorenz.cpp
// Generates time series for Lorenz model
// Ed Threlfall October 2022

#include "Lorenz.h"

Lorenz::Lorenz(double rho, double sigma, double beta, double timeStep) : mRho(rho), mSigma(sigma), mBeta(beta), mDt(timeStep)
{
    mX0 = 0.0;
    mY0 = 1.0;
    mZ0 = 0.0;

    mX = mX0;
    mY = mY0;
    mZ = mZ0;

    mXS = mXDot = 0.0;
    mYS = mYDot = 0.0;
    mZS = mZDot = 0.0;
}

Lorenz::~Lorenz()
{

}

void Lorenz::EvaluateTimestep()
{
    double A[NUMSTAGES];
    double B[NUMSTAGES];

    //scheme 2 Carpenter - Kennedy (NUMSTAGES must be 5)
    //see Fourth-Order 2N-Storage Runge-Kutta Schemes, NASA Technical Memorandum 109112
    A[0] = 0.0;
    A[1] = -0.4801594388478;
    A[2] = -1.4042471952;
    A[3] = -2.016477077503;
    A[4] = -1.056444269767;

    B[0] = 0.1028639988105;
    B[1] = 0.7408540575767;
    B[2] = 0.7426530946684;
    B[3] = 0.4694937902358;
    B[4] = 0.1881733382888;

    for (int iStage = 0; iStage < NUMSTAGES; iStage++)
    {
        EvaluateTimeDerivs();

        mXS = A[iStage] * mXS + mDt * mXDot;
        mYS = A[iStage] * mYS + mDt * mYDot;
        mZS = A[iStage] * mZS + mDt * mZDot;
        mX += B[iStage] * mXS;
        mY += B[iStage] * mYS;
        mZ += B[iStage] * mZS;
    }
}

void Lorenz::EvaluateTimeDerivs()
{
    mXDot = mSigma * (mY - mX);
    mYDot = mX * (mRho - mZ) - mY;
    mZDot = mX * mY - mBeta * mZ;
}

double Lorenz::OutputNusselt()
{
    return 1.0 + 2.0 * mY / mRho - (1.0+2.0*mY0/mRho);  //more interesting proxy, offset to start at zero
    //return 1.0 + 2.0 * mZ / gRho  - (1.0+2.0*mZ0/mRho);  //actual Nusselt number, offset to start at zero
}
