#include <iostream>
#include <string>
#include <iomanip>
#include <sstream>
#include <vector>
#include <fstream>
#include "parameters3.h"
#include "stdafx.h"


using namespace std;
#pragma once

const int Nn = (Nx + 1) * (Ny + 1) * (Nz + 1); // Number of total grid points
complex<Float> *psi = new complex<Float>[Nn];
Float *phi = new Float[Nn];  
Float dx, dy, dz;

// Here map all the 1D data back to 3D by spicifying the counting method
#define ijk(i,j,k) ((((i) * (Ny + 1)) + (j)) * (Nz + 1) + (k)) // Start from 0, increase by z then y then x. So (1,0,0) is the (Nz+1)*(Ny+1)-th point
#define psi(i,j,k) (psi[ijk(i,j,k)])
#define phi(i,j,k) (phi[ijk(i,j,k)])


void readpsiphi(string file){
    // define variables
	string fx, fy, fz, fpsi, fphi; //variables from file are here
	vector<float>V_x;
	vector<float>V_y;
	vector<float>V_z;
	vector<float>V_psi;
    vector<float>V_phi;
    //number of lines
    int n = 0;

    ifstream coeff(file); // opening the file
    if (coeff.is_open()) //if the file is opened
    {
        string line;
        getline(coeff,line);

        while (!coeff.eof()) //While the end is not reached
        {
            //I have 5 sets {x, y, z, psi, phi} so use 5 getlines
            getline(coeff, fx, ' ');
            V_x.push_back(stof(fx));
            getline(coeff, fy, ' ');
            V_y.push_back(stof(fy));
            getline(coeff, fz, ' ');
            V_z.push_back(stof(fz));
            getline(coeff, fpsi, ' ');
            V_psi.push_back(stof(fpsi));
            getline(coeff, fphi, '\n'); //new line after psi value
            V_phi.push_back(stof(fphi));

            n += 1; //increment number of lines
        }
        coeff.close(); //closing the file
        cout << "Finished reading! Number of entries: " << n << endl;    
    }
    else cout << "Unable to open file";

    for (int i = 0; i <= Nx; i++)
	{
		for (int j = 0; j <= Ny; j++)
			for (int k = 0; k <= Nz; k++)
		{
			phi(i, j, k) = V_phi[ijk(i, j, k)];
			psi(i, j, k) = V_psi[ijk(i, j, k)];
        }
    }
}

int main(){
    int i, j, k;
    dx = (xr - xl) / Nx;
	dy = (yr - yl) / Ny;
	dz = (zr - zl) / Nz;
    FILE *fileini;
    // Initial conditions
    printf("Setting initial conditions...\n");
    fflush(stdout);

	// Zero arrays
	memset(phi, 0, sizeof(phi));
	memset(psi, 0, sizeof(psi));

    // Read initial psi and phi from file
    readpsiphi("psi_phi.dat");

    fileini = fopen("psi_phi_new.dat", "w");

    for (i = 0; i <= Nx; i++)
    {
        for (j = 0; j <= Ny; j++)
            for (k = 0; k <= Nz; k++)
                fprintf(fileini, "%lg %lg %lg %lg %lg\n", xl + i * dx, yl + j * dy,
                                    zl + k * dz, psi(i, j, k), phi(i, j, k));
        fprintf(fileini, "\n");  // For Gnuplot
    }
    fclose(fileini);

    
}