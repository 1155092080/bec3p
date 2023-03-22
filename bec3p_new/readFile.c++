#include <iostream>
#include <string>
#include <iomanip>
#include <sstream>
#include <vector>
#include <fstream>
#include "stdafx.h"


using namespace std;
#define Float double
#define Nx 0//120
#define Ny 0//120
#define Nz 7//120

const int Nn = 8;//(Nx + 1) * (Ny + 1) * (Nz + 1); // Number of total grid points
complex<Float> *psi = new complex<Float>[Nn];
Float *phi = new Float[Nn];  
Float dx, dy, dz;
// Physical size of simulation volume in units of [L]
const Float xl = -3.0f, yl = -3.0f, zl = -3.0f;
const Float xr = 3.0f, yr = 3.0f, zr = 3.0f;

// Here map all the 1D data back to 3D by spicifying the counting method
#define ijk(i,j,k) ((((i) * (Ny + 1)) + (j)) * (Nz + 1) + (k)) // Start from 0, increase by z then y then x. So (1,0,0) is the (Nz+1)*(Ny+1)-th point
#define psi(i,j,k) (psi[ijk(i,j,k)])
#define phi(i,j,k) (phi[ijk(i,j,k)])


void readpsiphi(string file){
    // define variables
	string fx, fy, fz, fpsi, fphi; //variables from file are here
	vector<Float>V_x;
	vector<Float>V_y;
	vector<Float>V_z;
	vector<Float>V_psi;
    vector<Float>V_phi;
    //number of lines
    int n = 1;

    ifstream coeff(file); // opening the file
    if (coeff.is_open()) //if the file is opened
    {
        string line;
        cout << "file opened!" <<endl;
        while (!coeff.eof()) //While the end is not reached
        {
            //I have 5 sets {x, y, z, psi, phi} so use 5 getlines
            getline(coeff, fx, '\t');
            cout << n <<endl;
            V_x.push_back(stod(fx));
            getline(coeff, fy, ' ');
            V_y.push_back(stod(fy));
            getline(coeff, fz, ' ');
            V_z.push_back(stod(fz));
            getline(coeff, fpsi, ' ');
            V_psi.push_back(stod(fpsi));
            getline(coeff, fphi, '\n'); //new line after psi value
            V_phi.push_back(stod(fphi));

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

//*********************************************************************
// Read a file as the initial state
//*********************************************************************
void readdouble(string file){
    Float *f_x = new Float[Nn];
    Float *f_y = new Float[Nn];
    Float *f_z = new Float[Nn];
	complex<Float> *f_psi = new complex<Float>[Nn];
    Float *f_phi = new Float[Nn];
    ifstream ifs(file, ios::in); // opening the file
    if (!ifs.is_open())
    {
        cout << "open file fail!" << endl;
    } 
    else
    {
        cout << "open file successful!" << endl;
        for (int i = 0; i < Nn; i++)
        {
            ifs >> f_x[i] >> f_y[i] >> f_z[i] >> f_psi[i] >> f_phi[i];
        }
        ifs.close();
        cout << "Finished reading! Number of entries: " << Nn << endl;
    }
    for (int i = 0; i <= Nx; i++)
	{
		for (int j = 0; j <= Ny; j++)
			for (int k = 0; k <= Nz; k++)
		{
			phi(i, j, k) = f_phi[ijk(i, j, k)];
			psi(i, j, k) = f_psi[ijk(i, j, k)];
        }
    }

}
int main(){
    int i, j, k;
    dx = (xr - xl) / 1;
	dy = (yr - yl) / 1;
	dz = (zr - zl) / 1;
    FILE *fileini;
    // Initial conditions
    printf("Setting initial conditions...\n");
    fflush(stdout);

	// Zero arrays
    memset(phi, 0, sizeof(phi));
    memset(psi, 0, sizeof(psi));

    // Read initial psi and phi from file
    cout << "Start reading" << endl;
    readdouble("test.dat");

    fileini = fopen("psi_phi_new.dat", "w");

    for (i = 0; i <= Nx; i++)
    {
        for (j = 0; j <= Ny; j++)
            for (k = 0; k <= Nz; k++)
                fprintf(fileini, "%e %e %e %e %e\n", xl + i * dx, yl + j * dy,
                                    zl + k * dz, psi(i, j, k), phi(i, j, k));
        fprintf(fileini, "\n");  // For Gnuplot
    }
    fclose(fileini);

    
}
