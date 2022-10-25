// CoupledTurbulenceDiffusion1D.cpp
// Simple test of one way coupling from Lorenz model output to 1D diffusion
// contains `main' function
// Requires input file "inputs.txt" with parameters rho, sigma, beta, diffusion coeff
// e.g. 
// 28.000000 10.000000 2.666667 1.0
// Ed Threlfall October 2022

#include <iostream>

#include "Lorenz.h"
#include "Diffusion.h"

int main()
{
    std::cout << "1-way coupled Lorenz model and diffusion time series generator.\n";

    //read parameters
    double rho, sigma, beta;  //Lorenz model parameters
    double diffusionCoeff;

    FILE* paramFile = fopen("inputs.txt", "r");
    fscanf(paramFile, "%lf %lf %lf %lf", &rho, &sigma, &beta, &diffusionCoeff);
    fclose(paramFile);

    double timeStep = 0.001;
    int numSteps = 100000;
    int outputFrequency = 10;

    Lorenz lorenz = Lorenz(rho, sigma, beta, timeStep);

    Diffusion diffusion = Diffusion(diffusionCoeff, timeStep);

    for (int i = 0; i < numSteps; i++)
    {
        if (!(i % outputFrequency))
        {
            diffusion.OutputData();
        }

        lorenz.EvaluateTimestep();
        diffusion.EvaluateTimestep(lorenz.OutputNusselt());
    }

    return 1;
}

