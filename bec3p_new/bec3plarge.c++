// BEC3P.cpp : Self-gravitating Bose-Einstein condensate
//
// This version dated 2013/08/05.

#include "stdafx.h"
#include "parameters3large.h"
#include "mkpath.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <limits>
#include <cmath>

#ifdef USECL
#include <CL/opencl.h>
#include <SDKCommon.hpp>
#include <SDKFile.hpp>
#endif

#pragma warning(disable:4996)

using namespace std;
using std::string;
const int Nn = (Nx + 1) * (Ny + 1) * (Nz + 1); // Number of total grid points

#define NRMN 100 // Check recent NRMN iterations to see if the loop is oscillatory
#define NRMC 1000 // Maximum iterations for each direction evolution
#define NRMP 10000 // Maximum iterations for Poisson evolution
#define NTOL 1e-7 // Relative error of oscillatory loop checking. Should be smaller than tolGPE

complex<Float> eye, dt;
Float t, dx, dy, dz, idx2, idy2, idz2, dxl, dyl, dzl, idxl2, idyl2, idzl2;
string path = prefix;

// Here flatten all the 3D data into 1D by ijk(i,j,k) function
complex<Float> *psi = new complex<Float>[Nn];  
complex<Float> *psi_n = new complex<Float>[Nn];
Float *dens = new Float[Nn];
Float *phi = new Float[Nn];
Float *phiU = new Float[Nn];
Float *UU = new Float[Nn];
Float *res = new Float[Nn];
Float *phiBary = new Float[Nn];
Float *xgrid = new Float[Nx + 1];
Float *ygrid = new Float[Ny + 1];
Float *zgrid = new Float[Nz + 1];
vector<Float>r_interp;
vector<Float>psi_interp;

// Here map all the 1D data back to 3D by spicifying the counting method
#define ijk(i,j,k) ((((i) * (Ny + 1)) + (j)) * (Nz + 1) + (k)) // Start from 0, increase by z then y then x. So (1,0,0) is the (Nz+1)*(Ny+1)-th point
#define psi(i,j,k) (psi[ijk(i,j,k)])
#define psi_n(i,j,k) (psi_n[ijk(i,j,k)])
#define density(i,j,k) (dens[ijk(i,j,k)])
#define phi(i,j,k) (phi[ijk(i,j,k)])
#define phiU(i,j,k) (phiU[ijk(i,j,k)])
#define U(i,j,k) (UU[ijk(i,j,k)])
#define res(i,j,k) (res[ijk(i,j,k)])
#define phiBary(i,j,k) (phiBary[ijk(i,j,k)])

// Forward declarations
Float init(int, int, int);
Float fermi(Float, int, int, int);
Float BaryU(int, int, int);
Float DMiniphi(int i, int j, int k);	
Float get_normsimp();
void get_cond(Float &, Float &, Float &, Float, Float &);

void get_phi();
void get_Vtr();
void get_U(Float);
void calc_rhs_x(complex<Float>, complex<Float>,
				complex<Float>, complex<Float>,
				complex<Float>, complex<Float>,
				complex<Float>, complex<Float>);
void solve_x(complex<Float>, complex<Float>,
			 complex<Float>, complex<Float>,
			 complex<Float>, complex<Float>,
			 complex<Float>, complex<Float>);
void calc_rhs_y(complex<Float>, complex<Float>,
				complex<Float>, complex<Float>,
				complex<Float>, complex<Float>,
				complex<Float>, complex<Float>);
void solve_y(complex<Float>, complex<Float>,
			 complex<Float>, complex<Float>,
			 complex<Float>, complex<Float>,
			 complex<Float>, complex<Float>);
void calc_rhs_z(complex<Float>, complex<Float>, complex<Float>,
				complex<Float>, complex<Float>, complex<Float>);
void solve_z(complex<Float>, complex<Float>, complex<Float>,
				complex<Float>, complex<Float>, complex<Float>);
Float energy(Float, FILE*);
void movie(int);
void thomas(complex<Float> *, complex<Float> *, complex<Float> *,
			complex<Float> *, int m);
void get_density();
void readdouble(std::string file);
Float interpini(Float r);
void readinterp(string file);
template<typename T>
int nearestNeighbourIndex(std::vector<T>& x, T& value);
template<typename T>
std::vector<T> interp1(std::vector<T>& x, std::vector<T>& y, std::vector<T>& x_new);

//**********************************************************************
void loopdetect(Float *nrma, Float norm, const char *pdir, int &nrmc)
{
	int i, j;
	int nrmn = nrmc++; // Current iteration number

	if (nrmc > NRMC) // If current iteration number exceeds the maximum, exit and show divergence
	{
		fprintf(stderr, "No convergence after %d "
				 "iterations in the %s direction.\n", nrmc, pdir);
		exit(1);
	}

	if (nrmn >= NRMN) // If current iteration number exceeds the nrma memory, shift the array
	{
		memmove(nrma, &(nrma[1]), sizeof(nrma) - sizeof(nrma[0]));
		nrmn = NRMN - 1;
	}
	nrma[nrmn] = norm;

	for (i = nrmn - 1; i >= 0; i--)
	{
		//cout << "Loop detect " << fabs(norm-nrma[i]) <<endl;
		if (fabs(norm - nrma[i]) < fabs(NTOL * norm))
		{
			fprintf(stderr, "Oscillatory loop in the %s direction is "
						 "detected, solution is not convergent.\n", pdir);
			for (j = i; j <= nrmn; j++)
				fprintf(stderr, "norm[%d]=%lg\n", j, nrma[j]);
			exit(1);
		}
	}
}

//**********************************************************************
int _tmain(int argc, _TCHAR* argv[])
{
	complex<Float> xi;
	Float omega, gamma, mu;
	Float norm, norm0;
	Float E0, E1 = 0, E2 = 0;
	int i, j, k, itime, ktime;
	FILE *fileini, *fileerg, *file21, *file22, *file23, *file24, *file_current;
	Float nrma[NRMN];
	int nrmc;
	complex<Float> foo1X, foo1XS, foo1Y, foo1YS, foo1ZS;
	complex<Float> foo2XS, foo2YS;
	complex<Float> foo3X, foo3XS, foo3Y, foo3YS, foo3Z, foo3ZS;
	complex<Float> foo4;
	complex<Float> foo5X, foo5Y, foo5Z;
    complex<Float> foo1Xl, foo1XSl, foo1Yl, foo1YSl, foo1ZSl;
	complex<Float> foo2XSl, foo2YSl;
	complex<Float> foo3Xl, foo3XSl, foo3Yl, foo3YSl, foo3Zl, foo3ZSl;
	complex<Float> foo4l;
	complex<Float> foo5Xl, foo5Yl, foo5Zl;
	double norm_ini;
	int mkdirretval;
    mkdirretval=light::mkpath(path.c_str());
    std::cout << mkdirretval << '\n';
	string filepath;
	bool imagt = true;
#ifdef BARY
printf("DM Simulation with baryonic matter potential...\n");
fflush(stdout);
#endif

#ifdef USECL
	initializeCL();
#endif

	// Fixed grad parameters for fine part

	eye = complex<Float>(0, 1);
	dx = (xr - xl) / Nxf;
	dy = (yr - yl) / Nyf;
	dz = (zr - zl) / Nzf;

	idx2 = 1 / dx / dx;
	idy2 = 1 / dy / dy;
	idz2 = 1 / dz / dz;

    // Fixed grad parameters for large part

    dxl = (xl - xll) / Nxl;
	dyl = (yl - yll) / Nyl;
	dzl = (zl - zll) / Nzl;

	idxl2 = 1 / dxl / dxl;
	idyl2 = 1 / dyl / dyl;
	idzl2 = 1 / dzl / dzl;

	norm0 = N;
     
	// Initial conditions
printf("Setting initial conditions...\n");
fflush(stdout);

	// Zero arrays
	memset(phi, 0, sizeof(phi));
	memset(psi, 0, sizeof(psi));
	memset(psi_n, 0, sizeof(psi));
	memset(UU, 0, sizeof(UU));
	memset(res, 0, sizeof(res));
	memset(dens, 0, sizeof(dens));
#ifdef REALTIME // Real-time evolution directly from initial conditions
	dt = complex<Float>(tau, 0.);
	imagt = false;
#else
	dt = complex<Float>(0., -tau);
#endif
	omega = omega0;
	gamma = 0;
#ifdef GRAV
	mu = 0 ;
#else
	mu = 0;//0.5 / SQ(R);
#endif

	xi = (eye + gamma) / (gamma * gamma + 1);

printf("Setting up helper variables...\n");
fflush(stdout);

	// Helper variables for fine part that need to be computed only once
	foo1X = eye * omega * xi * dt / (4 * dx);
	foo1XS = -xi * dt * idx2 / (Float)4;
	foo1Y = eye * omega * xi * dt / (4 * dy);
	foo1YS = -xi * dt * idy2 / (Float)4;
	foo1ZS = -xi * dt * idz2 / (Float)4;
	foo2XS = eye * omega * xi * dt / (4 * dx);
	foo2YS = eye * omega * xi * dt / (4 * dy);
	foo3X = xi * dt * idx2 / (Float)4;
	foo3XS = (Float)1 + xi * dt * idx2 / (Float)2;
	foo3Y = xi * dt * idy2 / (Float)4;
	foo3YS = (Float)1 + xi * dt * idy2 / (Float)2;
	foo3Z = xi * dt * idz2 / (Float)4;
	foo3ZS = (Float)1 + xi * dt * idz2 / (Float)2;
	foo4 = xi * dt / (Float)6;
	foo5X = (Float)1 - xi * dt * idx2 / (Float)2;
	foo5Y = (Float)1 - xi * dt * idy2 / (Float)2;
	foo5Z = (Float)1 - xi * dt * idz2 / (Float)2;

    // Helper variables for large part that need to be computed only once
	foo1Xl = eye * omega * xi * dt / (4 * dxl);
	foo1XSl = -xi * dt * idxl2 / (Float)4;
	foo1Yl = eye * omega * xi * dt / (4 * dyl);
	foo1YSl = -xi * dt * idyl2 / (Float)4;
	foo1ZSl = -xi * dt * idzl2 / (Float)4;
	foo2XSl = eye * omega * xi * dt / (4 * dxl);
	foo2YSl = eye * omega * xi * dt / (4 * dyl);
	foo3Xl = xi * dt * idxl2 / (Float)4;
	foo3XSl = (Float)1 + xi * dt * idxl2 / (Float)2;
	foo3Yl = xi * dt * idyl2 / (Float)4;
	foo3YSl = (Float)1 + xi * dt * idyl2 / (Float)2;
	foo3Zl = xi * dt * idzl2 / (Float)4;
	foo3ZSl = (Float)1 + xi * dt * idzl2 / (Float)2;
	foo4l = xi * dt / (Float)6;
	foo5Xl = (Float)1 - xi * dt * idxl2 / (Float)2;
	foo5Yl = (Float)1 - xi * dt * idyl2 / (Float)2;
	foo5Zl = (Float)1 - xi * dt * idzl2 / (Float)2;
#ifdef INIFILE
readdouble(inifile);
#endif
#ifdef INTERP
readinterp(interpfile);
#endif
	// Initial state
	filepath = path + "psi_ini.dat";
	fileini = fopen(filepath.c_str(), "w");
	t = 0.0;
	itime = 0;
	ktime = 0;
	// Initial all x grid
	for (i = 0; i <= Nx; i++)
	{
		if (Nxl<=i && i<=Nxl+Nxf){
			xgrid[i] = xl + (i-Nxl) * dx;
			}
		else if (i<Nxl){
			xgrid[i] = xl + (i-Nxl)*dxl;
			}
		else{
			xgrid[i] = xr + (i-Nxl-Nxf) * dxl;
			}

	}
	// Initial all y grid
	for (j = 0; j <= Ny; j++)
	{
		if (Nyl<=j && j<=Nyl+Nyf){
			ygrid[j] = yl + (j-Nyl) * dy;
			}
		else if (j<Nyl){
			ygrid[j] = yl + (j-Nyl)*dyl;
			}
		else{
			ygrid[j] = yr + (j-Nyl-Nyf) * dyl;
			}
	}
	// Initial all z grid
	for (k = 0; k <= Nz; k++)
	{
		if (Nzl<=k && k<=Nzl+Nzf){
			zgrid[k] = zl + (k-Nzl) * dz;
			}
		else if (k<Nzl){
			zgrid[k] = zl + (k-Nzl)*dzl;
			}
		else{
			zgrid[k] = zr + (k-Nzl-Nzf) * dzl;
			}
	}

	for (i = 0; i <= Nx; i++)
	{
		for (j = 0; j <= Ny; j++)
			for (k = 0; k <= Nz; k++)
		{
#ifdef GRAV
			Float x2 = SQ(xgrid[i]);
			Float y2 = SQ(ygrid[j]);
			Float z2 = SQ(zgrid[k]);
			Float r = sqrt(x2 + y2 + z2);
			#ifdef INTERP
			Float rho = interpini(r);
			#else
			Float rho = init(i, j, k);
			#endif
#ifndef INIFILE
			phi(i, j, k) = DMiniphi(i, j, k); //(Float)(-G * N / (r > (.25 * dx) ? r : .5 * dx));
			psi(i, j, k) = complex<Float>(sqrt(rho), 0);//complex initialization
#endif
			//phi(i, j, k) = DMiniphi(i, j, k); //(Float)(-G * N / (r > (.25 * dx) ? r : .5 * dx));
			phiBary(i,j,k) = 0.0;// BaryU(i, j, k);
			
#else
			psi(i, j, k) = complex<Float>(fermi(mu, i, j, k), 0);
#endif
			
		}		
	}
	norm_ini = get_normsimp();
	Float renorm = sqrt(N / norm_ini);
	for (i = 0; i <= Nx; i++)
	{
		for (j = 0; j <= Ny; j++)
			for (k = 0; k <= Nz; k++)
			{
				psi(i, j, k) *= renorm;
			}
	}	
	// 		fprintf(fileini, "\n");	// For Gnuplot		
	// }
	// fclose(fileini);	
	norm_ini = get_normsimp();
	printf("Initial norm is P=%11.4lg\n", norm_ini);
	fflush(stdout);
printf("Setting up boundary conditions...\n");
fflush(stdout);

	// Boundary conditions psi=0
	for (i = 0; i <= Nx; i++)
		for (j = 0; j <= Ny; j++)
	{
		psi(i, j, 0) = 0;
		psi(i, j, Nz) = 0;
	}
	for (j = 0; j <= Ny; j++)
		for (k = 0; k <= Nz; k++)
	{
		psi(0, j, k) = 0;
		psi(Nx, j, k) = 0;
	}
	for (k = 0; k <= Nz; k++)
		for (i = 0; i <= Nx; i++)
	{
		psi(i, 0, k) = 0;
		psi(i, Ny, k) = 0;
	}

printf("Setting up initial density...\n");
fflush(stdout);

	get_density();

printf("Setting up initial potential...\n");
fflush(stdout);

#ifdef GRAV
	get_phi();
	// filepath = path + "psi_phi.dat";
	// file_current = fopen(filepath.c_str(), "w");

	// for (i = 0; i <= Nx; i++)
	// {
	// 	for (j = 0; j <= Ny; j++)
	// 		for (k = 0; k <= Nz; k++)
	// 		{
	// 			fprintf(file_current, "%e %e %e %e %e %e\n", xgrid[i],ygrid[j],zgrid[k], real(psi(i, j, k)), imag(psi(i, j, k)), phi(i, j, k));
	// 		}	
	// 		fprintf(file_current, "\n");	// For Gnuplot		
	// }
	// fclose(file_current);
#else
	get_Vtr();
#endif
	for (i = 0; i <= Nx; i++)
	{
		for (j = 0; j <= Ny; j++)
			for (k = 0; k <= Nz; k++)
			{
				fprintf(fileini, "%e %e %e %e %e %e\n", xgrid[i], ygrid[j], zgrid[k], real(psi(i, j, k)), imag(psi(i, j, k)), phi(i, j, k));
			}	
			fprintf(fileini, "\n");	// For Gnuplot		
	}
	fclose(fileini);
printf("Initiate the iteration...\n");
fflush(stdout);

	// get initial U 
	filepath = path + "ergini.dat";
	fileerg = fopen(filepath.c_str(), "w");
	get_U(mu);

	norm =  get_normsimp();
	printf("N=%6d, t=%11.4lg, E=%11.4lg, P=%11.4lg\n", 0, 0.0,
			energy(mu, fileerg), norm);
	fflush(stdout);

	movie(itime);	// Output data for contour plots

	filepath = path + "psi21.dat";
	file21 = fopen(filepath.c_str(), "w");
	filepath = path + "psi22.dat";
	file22 = fopen(filepath.c_str(), "w");
	filepath = path + "psi23.dat";
	file23 = fopen(filepath.c_str(), "w");
	filepath = path + "psi24.dat";
	file24 = fopen(filepath.c_str(), "w");

	// Time loop

printf("Starting the time loop...\n");
fflush(stdout);

	for (itime = 1, ktime = 0; itime - ktime <= time_n; itime++)
	{
		bool bLoop;
		t += real(dt);
		//if ((itime - ktime) == despin_n) omega = 0.0;

		// Find psi'=Rx(psi_old)
		calc_rhs_x(foo1X, foo3X, foo4, foo5X, foo1Xl, foo3Xl, foo4l, foo5Xl);

		for (bLoop = true, nrmc = 0; bLoop;)// Solve Lx(psi'')=psi' and iterate
		{									// due to nonlinear term in U
			solve_x(foo1XS, foo2XS, foo3XS, foo4, foo1XSl, foo2XSl, foo3XSl, foo4l);
			get_density();
#ifdef GRAV
			get_phi();
#endif
			get_U(mu);
			Float norm1 = get_normsimp();
			//cout << "Normx Relerr: " << fabs((norm1 -  norm)/norm) << " " << norm << endl;
			if (fabs(norm1 - norm) < fabs(tolGPE * norm)) bLoop = false;
			else loopdetect(nrma, norm1, "X", nrmc);
			norm = norm1;
		}
#ifdef SHOW_LOOPS
	fprintf(stderr, "H_x required %d interations.\n", nrmc);
	fflush(stderr);
#endif
		// Find psi'''=Ry(psi'')
		calc_rhs_y(foo1Y, foo3Y, foo4, foo5Y, foo1Yl, foo3Yl, foo4l, foo5Yl);
		for (bLoop = true, nrmc = 0; bLoop;)
		{	// Solve Ly(psi'''')=psi'''
			solve_y(foo1YS, foo2YS, foo3YS, foo4, foo1YSl, foo2YSl, foo3YSl, foo4l);
			get_density();
#ifdef GRAV
			get_phi();
#endif
			get_U(mu);
			Float norm1 = get_normsimp();
			//cout <<"Normy Relerr: " <<  fabs((norm1-norm)/norm) << " " << norm << endl;
			if (fabs(norm1 - norm) < fabs(tolGPE * norm)) bLoop = false;
			else loopdetect(nrma, norm1, "Y", nrmc);
			norm = norm1;
		}
#ifdef SHOW_LOOPS
	fprintf(stderr, "H_y required %d interations.\n", nrmc);
	fflush(stderr);
#endif
		calc_rhs_z(foo3Z, foo4, foo5Z, foo3Zl, foo4l, foo5Zl);		// Find psi'''''=Rz(psi'''')
		for (bLoop = true, nrmc = 0; bLoop;)
		{
			solve_z(foo1ZS, foo3ZS, foo4, foo1ZSl, foo3ZSl, foo4l);	// Solve Lz(psi_new)=psi'''''
			get_density();
#ifdef GRAV
			get_phi();
#endif
			get_U(mu);
			Float norm1 = get_normsimp();
			//cout << "Normz Relerr: " << fabs((norm1-norm)/norm) << " " << norm << endl;
			if (fabs(norm1 - norm) < fabs(tolGPE * norm)) bLoop = false;
			else loopdetect(nrma, norm1, "Z", nrmc);
			norm = norm1;
		}
#ifdef SHOW_LOOPS
	fprintf(stderr, "H_z required %d interations.\n", nrmc);
	fflush(stderr);
#endif
if (imagt){
			Float renorm = sqrt(N / norm);
			for (i = 0; i <= Nx; i++)
				for (j = 0; j <= Ny; j++)
					for (k = 0; k <= Nz; k++)
				psi(i, j, k) *= renorm;
			// norm = get_normsimp();
			// printf("Checking at N=%6d, t=%11.4lg, P=%11.4lg\n", itime, t, norm);
			// fflush(stdout);
			}
		E0 = energy(mu, fileerg);
		printf("N=%6d, t=%11.4lg, E=%11.4lg, Erel=%11.4lg, P=%11.4lg\n", itime, t, E0, fabs((E0-E1)/E0), norm);
		fflush(stdout);

		if (itime > 0 && itime % nstep1 == 0)
		{
			filepath = path + "psi_phi.dat";
			file_current = fopen(filepath.c_str(), "w");

			for (i = 0; i <= Nx; i++)
			{
				for (j = 0; j <= Ny; j++)
					for (k = 0; k <= Nz; k++)
					{
						fprintf(file_current, "%e %e %e %e %e %e\n", xgrid[i],ygrid[j],zgrid[k], real(psi(i, j, k)), imag(psi(i, j, k)), phi(i, j, k));
					}	
					fprintf(file_current, "\n");	// For Gnuplot		
			}
			fclose(file_current);
		}
		if (itime > nstep0 && itime % nstep2 == 0)
		{
			movie(itime);			// Output data for contour plots
		}

		{
			Float x2, y2, alpha1, elz;

			get_cond(x2, y2, alpha1, omega, elz);

			fprintf(file21, "%lg %lg\n", t, alpha1);
			fprintf(file22, "%lg %lg\n", t, elz);
			fprintf(file23, "%lg %lg\n", t, E0);
			fprintf(file24, "%lg %lg\n", t, norm);

			fflush(file21);
			fflush(file22);
			fflush(file23);
			fflush(file24);
		}

		if (real(dt) < 1e-15)
		{
			//cout << "Energy Relerr: " << fabs((E0 - E1)/E0) << " " << E0 << endl;
			if (fabs(E0 - E1) < tolREL * fabs(E0) &&
				fabs(2 * E1 - E0 - E2) < tolREL * fabs(E0))
			{
				filepath = path + "psi_phi_ground.dat";
				file_current = fopen(filepath.c_str(), "w");

				for (i = 0; i <= Nx; i++)
				{
					for (j = 0; j <= Ny; j++)
						for (k = 0; k <= Nz; k++)
							fprintf(file_current, "%e %e %e %e %e %e\n", xgrid[i],ygrid[j],zgrid[k], real(psi(i, j, k)), imag(psi(i, j, k)), phi(i, j, k));
					fprintf(file_current, "\n");  // For Gnuplot
				}
				fclose(file_current);
				dt = tau;
				omega = omega0;
				gamma = gamma0;
				imagt = false;
#ifndef GRAV
				mu = 0;
#endif
				xi = (eye + gamma) / (gamma * gamma + 1);

				// Recompute helper variables
				foo1X = eye * omega * xi * dt / (4 * dx);
				foo1XS = -xi * dt * idx2 / (Float)4;
				foo1Y = eye * omega * xi * dt / (4 * dy);
				foo1YS = -xi * dt * idy2 / (Float)4;
				foo1ZS = -xi * dt * idz2 / (Float)4;
				foo2XS = eye * omega * xi * dt / (4 * dx);
				foo2YS = eye * omega * xi * dt / (4 * dy);
				foo3X = xi * dt * idx2 / (Float)4;
				foo3XS = (Float)1 + xi * dt * idx2 / (Float)2;
				foo3Y = xi * dt * idy2 / (Float)4;
				foo3YS = (Float)1 + xi * dt * idy2 / (Float)2;
				foo3Z = xi * dt * idz2 / (Float)4;
				foo3ZS = (Float)1 + xi * dt * idz2 / (Float)2;
				foo4 = xi * dt / (Float)6;
				foo5X = (Float)1 - xi * dt * idx2 / (Float)2;
				foo5Y = (Float)1 - xi * dt * idy2 / (Float)2;
				foo5Z = (Float)1 - xi * dt * idz2 / (Float)2;

				// Helper variables for large part
				foo1Xl = eye * omega * xi * dt / (4 * dxl);
				foo1XSl = -xi * dt * idxl2 / (Float)4;
				foo1Yl = eye * omega * xi * dt / (4 * dyl);
				foo1YSl = -xi * dt * idyl2 / (Float)4;
				foo1ZSl = -xi * dt * idzl2 / (Float)4;
				foo2XSl = eye * omega * xi * dt / (4 * dxl);
				foo2YSl = eye * omega * xi * dt / (4 * dyl);
				foo3Xl = xi * dt * idxl2 / (Float)4;
				foo3XSl = (Float)1 + xi * dt * idxl2 / (Float)2;
				foo3Yl = xi * dt * idyl2 / (Float)4;
				foo3YSl = (Float)1 + xi * dt * idyl2 / (Float)2;
				foo3Zl = xi * dt * idzl2 / (Float)4;
				foo3ZSl = (Float)1 + xi * dt * idzl2 / (Float)2;
				foo4l = xi * dt / (Float)6;
				foo5Xl = (Float)1 - xi * dt * idxl2 / (Float)2;
				foo5Yl = (Float)1 - xi * dt * idyl2 / (Float)2;
				foo5Zl = (Float)1 - xi * dt * idzl2 / (Float)2;
			}
			else
			{
				E2 = E1;
				E1 = E0;
				ktime++;
			}

#ifdef GRAV
			get_phi();
#endif
			get_U(mu);
		}

	}	// Close time loop

	fclose(file21);
	fclose(file22);
	fclose(file23);
	fclose(file24);
	fclose(fileerg);

#ifdef USECL
	finalizeCL();
#endif

	return true;
}

//*********************************************************************
// Initial TF density profile
//*********************************************************************
Float init(int i, int j, int k)		
{
	Float F, x, y, z, r, rh1, xrh1, rh2, xrh2;

	x = xgrid[i];
	y = ygrid[j];
	z = zgrid[k];
	r = sqrt((1 + ex) * x * x + (1 + ey) * y * y + (1 + ez) * z * z);
	rh1 = 3e5;
	xrh1 = r/0.707/rh1;
	rh2 = 23;
	xrh2 = r/0.707/rh2;
	
	F = 0.0;
	F = 1./(1.+3*SQ(xrh1))/1e10 + 1./pow(1.+3*SQ(xrh2),8);
	F *= (Float)rho0;


	return F;
}

//*********************************************************************
// Initial interpolated density profile
//*********************************************************************
Float interpini(Float r)		
{
	vector<Float> newx{r};
	vector<Float> res = interp1(r_interp, psi_interp, newx);

	return res[0];
}

//*********************************************************************
// Read a file as the initial state
//*********************************************************************
void readdouble(string file){
    Float *f_x = new Float[Nn];
    Float *f_y = new Float[Nn];
    Float *f_z = new Float[Nn];
	Float *f_psi_real = new Float[Nn];
    Float *f_psi_imag = new Float[Nn];
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
            ifs >> f_x[i] >> f_y[i] >> f_z[i] >> f_psi_real[i] >> f_psi_imag[i] >> f_phi[i];
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
			psi(i, j, k) = {f_psi_real[ijk(i, j, k)], f_psi_imag[ijk(i, j, k)]};
        }
    }

}

//*********************************************************************
// Initial guess of DM gravitational potential from a homogeneous 
// density profile
//*********************************************************************
Float DMiniphi(int i, int j, int k)		
{
	Float F, x, y, z, r;

	x = xgrid[i];
	y = ygrid[j];
	z = zgrid[k];
	r = sqrt((1 + ex) * x * x + (1 + ey) * y * y + (1 + ez) * z * z);
	F = 0.0;
	if (r<=R) F = (Float)(-G * N *(3*SQ(R)-SQ(r))/ (2*R*R*R));
	else F = (Float)(-G * N / r);
	return F;
}

//*********************************************************************
// Constant Baryonic gravitational potential
// Bulge Core component only
//*********************************************************************
Float BaryU(int i, int j, int k)		
{
	Float F, x, y, z, r, Mc1, Mc2, rc1, rc2;
        
        Mc1 = 0.0;//1.14715;
        Mc2 = 61181.5;//6.11815;
        rc1 = 2.7;
        rc2 = 0.42;

	x = xgrid[i];
	y = ygrid[j];
	z = zgrid[k];
 	r = sqrt((1 + ex) * x * x + (1 + ey) * y * y + (1 + ez) * z * z);
	F = 0.0;
 	F = (Float)((-G * Mc1 / sqrt( SQ(rc1) + SQ(r) ) ) + (-G * Mc2 / sqrt( SQ(rc2) + SQ(r) ) ));
 	return F;
 }

//Float BaryU(int i, int j, int k)
//{
//	Float F;
	//const Float d = (Float)0.25 * (SQ(xr - xl) + SQ(yr - yl) + SQ(zr - zl));
//	F = (Float)0.5 *SQ(omg)*((1 + ex) * SQ(xl + i * dx) +
//								(1 + ey) * SQ(yl + j * dy) +
//								(1 + ez) * SQ(zl + k * dz));
//	return F;
	
//}

//**********************************************************************
// Energy in lab frame, i.e. not calculating the rotational energy
//**********************************************************************
Float energy(Float mu, FILE *fileerg)
{
	Float EK, EU, EI, EG, mdpsisq, E;
	int i, j, k;
	const static Float dV = dx * dy * dz;

	EK = EU = EI = EG = 0;
	for (k = Nzl+1; k < Nzl+Nzf; k++)
		for (j = Nyl+1; j < Nyl+Nyf; j++)
			for (i = Nxl+1; i < Nxl+Nxf; i++)
	{
		mdpsisq = SQ(real(psi(i, j, k))) + SQ(imag(psi(i, j, k))); // density rho
		EI += SQ(mdpsisq); // rho^2
		EG += 0.5 * (phi(i, j, k)) * mdpsisq; 
		EU += (
#ifdef GRAV
			0.5 *		// See Wang, PRD64, 124009
#endif
			phiU(i, j, k) - mu) * mdpsisq;
		EK += SQ(real(psi(i + 1, j, k) - psi(i - 1, j, k))) * idx2; // kinetic energy in  
		EK += SQ(imag(psi(i + 1, j, k) - psi(i - 1, j, k))) * idx2;
		EK += SQ(real(psi(i, j + 1, k) - psi(i, j - 1, k))) * idy2;
		EK += SQ(imag(psi(i, j + 1, k) - psi(i, j - 1, k))) * idy2;
		EK += SQ(real(psi(i, j, k + 1) - psi(i, j, k - 1))) * idz2;
		EK += SQ(imag(psi(i, j, k + 1) - psi(i, j, k - 1))) * idz2;
	}
	EI *= (Float)0.5 * c * dV;
	EU *= dV;
	EG *= dV;
	EK *= (Float)0.125 * dV;  // 1/2*((\psi_i+1)-(\psi_i-1))/2dx)^2=1/8*...
	E = EK+EI+EU;
	fprintf(fileerg, "%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n", t, EK, EU, EG, EI, E);
	fflush(fileerg);
	return EK + EU + EI;
}

//**********************************************************************
//The function get_U computes the potential term (without rotation) in 
// the Gross-Pitaevskii equation: H = kinetic + U - Omega*L
//**********************************************************************
void get_U(Float mu)	// Find U
{
	int i, j, k;

	for (k = 0; k <= Nz; k++)
		for (j = 0; j <= Ny; j++)
			for (i = 0; i <= Nx; i++)
		U(i, j, k) = c * density(i, j, k) + phiU(i, j, k) - mu
#ifdef ISOTHERMAL
		+ kbT * log(density(i, j, k)+1e-20)
#endif
		;
}

//*********************************************************************
// Trapping potential
//*********************************************************************
void get_Vtr()	
{
	int i, j, k;
	//const Float d = (Float)0.25 * (SQ(xr - xl) + SQ(yr - yl) + SQ(zr - zl));
	const Float d = 0;
	for (k = 0; k <= Nz; k++)
		for (j = 0; j <= Ny; j++)
			for (i = 0; i <= Nx; i++)
	{
		//phi(i, j, k) = (Float)0.5 * ((1 + ex) * SQ(xl + i * dx) +
		//						(1 + ey) * SQ(yl + j * dy) +
		//						(1 + ez) * SQ(zl + k * dz) - d) / SQ(SQ(R));
		phi(i, j, k) = (Float)0.5 *SQ(omg)* ((1 + ex) * SQ(xgrid[i]) +
								(1 + ey) * SQ(ygrid[j]) +
								(1 + ez) * SQ(zgrid[k]));
	}
	}

//*********************************************************************
// Fermi-Thomas initial state, not used when GRAV is defined
//*********************************************************************
Float fermi(Float mu, int i, int j, int k)	
{
	Float x, y, z, r2, R2;
	// const Float norm = 15 * N * sqrt(2 * mu) * SQ(mu) / pi;
	const Float norm = sqrt(N)*pow(omg/pi, 3.0/4.0);
	x = xgrid[i];
	y = ygrid[j];
	z = zgrid[k];
	r2 = (1 + ex) * SQ(x) + (1 + ey) * SQ(y) + (1 + ez) * SQ(z);
	// R2 = SQ(R);
	//if (r2 < R2) return (Float)(1 - r2) * norm;//(Float)sqrt((0.5 * (R2 - r2)) * norm);
	//else return 0.0;
	return (Float)exp(-omg*r2/2.) * norm;
}

//**********************************************************************
// Calculate d_n, which is used to find psi_n+1 by Thomas algorirhm
// By ADI, we propogate along x for each y. So d/dy=0 and we fix j for each iteration
// Notedd - i*omega(xd/dy-yd/dx)/3 can be approximated by finite difference
//**********************************************************************
void calc_rhs_x(complex<Float> foo1,	// Find psi_n= Rx(psi)
			complex<Float> foo3, complex<Float> foo4, complex<Float> foo5,
			complex<Float> foo1l,
			complex<Float> foo3l, complex<Float> foo4l, complex<Float> foo5l)
{
	Float y;
	complex<Float> foo2;
	int i, j, k;

	for (i = 1; i < Nx; i++)
		for (j = 1; j < Ny; j++)
	{
		y = ygrid[j];
		for (k = 1; k < Nz; k++)
		{
		if (Nxl<i && i<Nxl+Nxf && Nyl<j && j<Nyl+Nyf && Nzl<k && k<Nzl+Nzf)		
		{
			foo2 = foo1 * y;
			psi_n(i, j, k) = (foo3 - foo2) * psi(i, j - 1, k) +
							 (foo5 - foo4 * U(i, j, k)) * psi(i, j, k) +
							 (foo3 + foo2) * psi(i, j + 1, k);
		}
		else
		{
			foo2 = foo1l * y;
			psi_n(i, j, k) = (foo3l - foo2) * psi(i - 1, j, k) +
							 (foo5l - foo4l * U(i, j, k)) * psi(i, j, k) +
							 (foo3l + foo2) * psi(i + 1, j, k);
		}
		}
	}
}

//**********************************************************************
// Calculate d_n, which is used to find psi_n+1 by Thomas algorirhm
// By ADI, we propogate along y for each x. So d/dx=0 and we fix i for each iteration
// Notedd - i*omega(xd/dy-yd/dx)/3 can be approximated by finite difference
//**********************************************************************
void calc_rhs_y(complex<Float> foo1,	// Find psi_n= Ry(psi)
			complex<Float> foo3, complex<Float> foo4, complex<Float> foo5,
			complex<Float> foo1l,
			complex<Float> foo3l, complex<Float> foo4l, complex<Float> foo5l)
{
	Float x;
	complex<Float> foo2;
	int i, j, k;

	for (i = 1; i < Nx; i++)
	{
		x = xgrid[i];
		for (j = 1; j < Ny; j++)
			for (k = 1; k < Nz; k++)
			{
			if (Nxl<i && i<Nxl+Nxf && Nyl<j && j<Nyl+Nyf && Nzl<k && k<Nzl+Nzf)
		{
			foo2 = foo1 * x;
			psi_n(i, j, k) = (foo3 + foo2) * psi(i, j - 1, k) +
							 (foo5 - foo4 * U(i, j, k)) * psi(i, j, k) +
							 (foo3 - foo2) * psi(i, j + 1, k);
		}
			else{
			foo2 = foo1l * x;
			psi_n(i, j, k) = (foo3l + foo2) * psi(i - 1, j, k) +
							 (foo5l - foo4l * U(i, j, k)) * psi(i, j, k) +
							 (foo3l - foo2) * psi(i + 1, j, k);
		}
			}
	}
}


//**********************************************************************
// Calculate d_n, which is used to find psi_n+1 by Thomas algorirhm
// By ADI, we propogate along z, so d/dx = d/dy = 0
//**********************************************************************
void calc_rhs_z(complex<Float> foo3,	// Find psi_n=Ry(psi)
				complex<Float> foo4, complex<Float> foo5,
				complex<Float> foo3l,
				complex<Float> foo4l, complex<Float> foo5l)
{
	int i, j, k;

	for (i = 1; i < Nx; i++)
		for (j = 1; j < Ny; j++)
			for (k = 1; k < Nz; k++)
			{
			if (Nxl<i && i<Nxl+Nxf && Nyl<j && j<Nyl+Nyf && Nzl<k && k<Nzl+Nzf)
		{
			psi_n(i, j, k) = foo3 * psi(i, j, k - 1) +
						 (foo5 - foo4 * U(i, j, k)) * psi(i, j, k) +
						 foo3 * psi(i, j, k + 1);
		}
			else{
			psi_n(i, j, k) = foo3l * psi(i, j, k - 1) +
						 (foo5l - foo4l * U(i, j, k)) * psi(i, j, k) +
						 foo3l * psi(i, j, k + 1);
		}
			}
}

//**********************************************************************
void solve_x(complex<Float> foo1,	// Solve Lx(psi)=psi_n for psi
			 complex<Float> foo2, complex<Float> foo3, complex<Float> foo4,
			 complex<Float> foo1l,
			 complex<Float> foo2l, complex<Float> foo3l, complex<Float> foo4l)
{
	Float y;
	complex<Float> l[Nx], d[Nx], u[Nx], r[Nx];
	int i, j, k;

	for (k = 1; k < Nz; k++)
		for (j = 1; j < Ny; j++)
	{
		y = ygrid[j];
		for (i = 1; i < Nx; i++)
		{
			if (Nxl<i && i<Nxl+Nxf && Nyl<j && j<Nyl+Nyf && Nzl<k && k<Nzl+Nzf){
				l[i] = foo1 + foo2 * y;
				d[i] = foo3 + foo4 * U(i, j, k);
				u[i] = foo1 - foo2 * y;
				r[i] = psi_n(i, j, k);
			}
			else{
				l[i] = foo1l + foo2l * y;
				d[i] = foo3l + foo4l * U(i, j, k);
				u[i] = foo1l - foo2l * y;
				r[i] = psi_n(i, j, k);
			}
		}
		thomas(l, d, u, r, Nx - 1);
		for (i = 1; i < Nx; i++) psi(i, j, k) = r[i];
	}
}

//**********************************************************************
void solve_y(complex<Float> foo1,	// Solve Ly(psi)=psi_n for psi
			 complex<Float> foo2, complex<Float> foo3, complex<Float> foo4,
			 complex<Float> foo1l,
			 complex<Float> foo2l, complex<Float> foo3l, complex<Float> foo4l)
{
	Float x;
	complex<Float> l[Ny], d[Ny], u[Ny], r[Ny];
	int i, j, k;

	for (i = 1; i < Nx; i++)
		for (k = 1; k < Nz; k++)
	{
		x = xgrid[i];
		for (j = 1; j < Ny; j++)
		{
			if (Nxl<i && i<Nxl+Nxf && Nyl<j && j<Nyl+Nyf && Nzl<k && k<Nzl+Nzf){
				l[j] = foo1 - foo2 * x;
				d[j] = foo3 + foo4 * U(i, j, k);
				u[j] = foo1 + foo2 * x;
				r[j] = psi_n(i, j, k);
			}
			else{
				l[j] = foo1l - foo2l * x;
				d[j] = foo3l + foo4l * U(i, j, k);
				u[j] = foo1l + foo2l * x;
				r[j] = psi_n(i, j, k);
			}
		}
		thomas(l, d, u, r, Ny - 1);
		for (j = 1; j < Ny; j++) psi(i, j, k) = r[j];
	}
}

//**********************************************************************
void solve_z(complex<Float> foo1,	// Solve Ly(psi)=psi_n for psi
			 complex<Float> foo3, complex<Float> foo4,
			 complex<Float> foo1l,	// Solve Ly(psi)=psi_n for psi
			 complex<Float> foo3l, complex<Float> foo4l)
{
	complex<Float> l[Nz], d[Nz], u[Nz], r[Nz];
	int i, j, k;

	for (i = 1; i < Nx; i++)
		for (j = 1; j < Ny; j++)
	{
		for (k = 1; k < Nz; k++)
		{
			if (Nxl<i && i<Nxl+Nxf && Nyl<j && j<Nyl+Nyf && Nzl<k && k<Nzl+Nzf){
				l[k] = u[k] = foo1;
				d[k] = foo3 + foo4 * U(i, j, k);
				r[k] = psi_n(i, j, k);
			}
			else{
				l[k] = u[k] = foo1l;
				d[k] = foo3l + foo4l * U(i, j, k);
				r[k] = psi_n(i, j, k);
			}
		}
		thomas(l, d, u, r, Nz - 1);
		for (k = 1; k < Nz; k++) psi(i, j, k) = r[k];
	}
}

//**********************************************************************
void get_cond(Float &x2, Float &y2, Float &alpha1,	// Find the mean values
			 Float omega, Float &elz)				// <x^2>, <y^2>, <elz>
{													// and the distortion
													// alpha of the condensate
	complex<Float> px, py, l;
	Float x, y;
	int i, j, k;

	x2 = 0.0;
	y2 = 0.0;
	elz = 0.0;
	l = complex<Float>(0, 0);

	for (i = Nxl; i <= Nxl+Nxf; i++)
	{
		x = xgrid[i];
		for (j = Nyl; j <= Nyl+Nyf; j++)
		{
			y = ygrid[j];
			for (k = Nzl; k <= Nzl+Nzf; k++)
			{
				x2 += x * x * (SQ(real(psi(i, j, k))) + SQ(imag(psi(i, j, k))));
				y2 += y * y * (SQ(real(psi(i, j, k))) + SQ(imag(psi(i, j, k))));
				if(i > 0 && i < Nx)
					px = (psi(i + 1, j, k) - psi(i - 1, j, k)) / (2 * dx);
				else if (i == 0) px = (psi(1, j, k) - psi(0, j, k)) / dx;
				else if (i == Nx) px = (psi(Nx, j, k) - psi(Nx - 1, j, k)) / dx;
				if (j > 0 && j < Ny)
					py = (psi(i, j + 1, k) - psi(i, j - 1, k)) / (2 * dy);
				else if (j == 0) py = (psi(i, 1, k) - psi(i, 0, k)) / dy;
				else if (j == Ny) py = (psi(i, Ny, k) - psi(i, Ny - 1, k)) / dy;
				l += conj(psi(i, j, k)) * eye * (x * py - y * px);
			}
		}
	}
	elz = sqrt(SQ(real(l)) + SQ(imag(l))) * dx * dy * dz;
	alpha1 = omega * ((x2 - y2) / (x2 + y2));
}

//**********************************************************************
void thomas(complex<Float> *l,	// Tridiagonal matrix solver with Thomas
			complex<Float> *d,	// algorithm. The arrays l,d,u and r contain
			complex<Float> *u,	// respectively the subdiagonal, the diagonal,
			complex<Float> *r,	// the superdiagonal and the right hand side.
			int m)				// m is the size of the system. At output the
{								// solution is put in r, while the original
								// right hand side is destroyed.
	int j;
	complex<Float> a;
    
	for (j = 2; j <= m; j++)
	{
		a = -l[j] / d[j - 1];
		d[j] = d[j] + a * u[j - 1];
		r[j] = r[j] + a * r[j - 1];
	}
	r[m] = r[m] / d[m];
	for (j = m - 1; j >= 1; j--)
		r[j] = (r[j] - u[j] * r[j + 1]) / d[j];
}

//***************************************************************************
Float get_normsimp()	// Find the norm using Simpson's rule
{
	static Float normz[Nx + 1][Ny + 1], normy[Nx + 1], norm;
	static Float normzl[Nx + 1][Ny + 1], normyl[Nx + 1], norml;
	int i, j, k;

	norm = 0.0;
	for (i = 0; i <= Nx; i++)
		for (j = 0; j <= Ny; j++)
	{
		normz[i][j] = SQ(real(psi(i, j, Nzl))) + SQ(imag(psi(i, j, Nzl))) + // Boundary index for inner box is Nzl~Nzl+Nzf, but we set right
					SQ(real(psi(i, j, Nzl+Nzf-1))) + SQ(imag(psi(i, j, Nzl+Nzf-1))); // boundary to be Nzl+Nzf-1 in order to make total number even
		normzl[i][j] = SQ(real(psi(i, j, 0))) + SQ(real(psi(i, j, 2*Nzl+Nzf-1)))+
					SQ(imag(psi(i, j, 0))) + SQ(imag(psi(i, j, 2*Nzl+Nzf-1))); // Boundary index for outer box is 0~2*Nzl+Nzf
		for (k = 1; k <= 2*Nzl+Nzf-3; k += 2) // third index from 1 to 2*Nzl+Nzf-2
			{
			if (Nzl+1<=k && k<=Nzl+Nzf-3) // third index from Nzl+1 to Nzl+Nzf-2
			normz[i][j] += 4 * (SQ(real(psi(i, j, k))) +
								SQ(imag(psi(i, j, k)))) +
						   2 * (SQ(real(psi(i, j, k + 1))) +
								SQ(imag(psi(i, j, k + 1))));
			else
			normzl[i][j] += (SQ(real(psi(i, j, k))) +
								SQ(imag(psi(i, j, k)))) +
						    (SQ(real(psi(i, j, k + 1))) +
								SQ(imag(psi(i, j, k + 1))));
			}
	}
	for (i = 0; i <=Ny; i++)
	{
		normy[i] = normz[i][Nyl] + normz[i][Nyl+Nyf-1];
		normyl[i] = normzl[i][0] + normzl[i][2*Nyl+Nyf-1];
		for (j = 1; j <= 2*Nyl+Nyf-3; j += 2)
		{
		if (Nyl+1<=j && j<=Nyl+Nyf-3){
			normy[i] += 4 * normz[i][j] + 2 * normz[i][j + 1];
		}
		else{
			normyl[i] += normzl[i][j] + normzl[i][j + 1];
		}
		}
		
	}

	norm = normy[Nxl] + normy[Nxl+Nxf-1];
	norml = normyl[0] + normyl[2*Nxl+Nxf-1];
	for (i = 1; i <= 2*Nxl+Nxf - 3; i += 2)
	{
	if (Nxl+1<=i && i<=Nxl+Nxf-3){
		norm += 4 * normy[i] + 2 * normy[i + 1];
	}
	else{
		norml += normyl[i] + normyl[i + 1];
	}
	}
	norm *= dx * dy * dz / 27;
	norml *= dxl * dyl * dzl;
	return norm+norml;
}
    
//**********************************************************************
void movieZ(int itime)	// Outputs files
{
	Float x, y, z, phpsi, depsi;
	char ch[17], cs[17], cu[17];
	int i, j, k;
	FILE *file8, *file10, *file11;
	string filepath;

	sprintf(ch, "densZ%07d.dat", itime);
	filepath = path + ch;
	file8 = fopen(filepath.c_str(), "w");
	sprintf(cs, "phasZ%07d.dat", itime);
	filepath = path + cs;
	file10 = fopen(filepath.c_str(), "w");
	sprintf(cu, "gravZ%07d.dat", itime);
	filepath = path + cu;
	file11 = fopen(filepath.c_str(), "w");

	// Output files for contour plots

	for (k = 0; k <= Nz; k += 2)
		for (j = 0; j <= Ny; j += 2)
	{
		for (i = 0; i <= Nx; i += 2)
		{
			x = xgrid[i];
			y = ygrid[j];
			z = zgrid[k];
			depsi = SQ(real(psi(i, j, k))) + SQ(imag(psi(i, j, k)));
			phpsi = atan2(imag(psi(i, j, k)) + 1e-15f, real(psi(i, j, k)));
			fprintf(file8, "%lf %lf %lf %lg\n", x, y, z, depsi);
			fprintf(file11, "%lf %lf %lf %lg\n", x, y, z, phiU(i, j, k));
			fprintf(file10, "%lf %lf %lf %lg\n", x, y, z, phpsi);
		}
		fprintf(file8, "\n");	// Leave blank line for Gnuplot before next j
		fprintf(file11, "\n");
		fprintf(file10, "\n");
	}

	fclose(file8);
	fclose(file10);
	fclose(file11);
}

//**********************************************************************
void movieX(int itime)	// Outputs files
{
	Float x, y, z, phpsi, depsi;
	char ch[17], cs[17], cu[17];
	int i, j, k;
	FILE *file8, *file10, *file11;
	string filepath;

	sprintf(ch, "densX%07d.dat", itime);
	filepath = path + ch;
	file8 = fopen(filepath.c_str(), "w");
	sprintf(cs, "phasX%07d.dat", itime);
	filepath = path + cs;
	file10 = fopen(filepath.c_str(), "w");
	sprintf(cu, "gravX%07d.dat", itime);
	filepath = path + cu;
	file11 = fopen(filepath.c_str(), "w");

	// Output files for contour plots

	for (i = 0; i <= Nx; i += 2)
		for (j = 0; j <= Ny; j += 2)
	{
		for (k = 0; k <= Nz; k += 2)
		{
			x = xgrid[i];
			y = ygrid[j];
			z = zgrid[k];
			depsi = SQ(real(psi(i, j, k))) + SQ(imag(psi(i, j, k)));
			phpsi = atan2(imag(psi(i, j, k)) + 1e-15f, real(psi(i, j, k)));
			fprintf(file8, "%lf %lf %lf %lg\n", x, y, z, depsi);
			fprintf(file11, "%lf %lf %lf %lg\n", x, y, z, phi(i, j, k));
			fprintf(file10, "%lf %lf %lf %lg\n", x, y, z, phpsi);
		}
		fprintf(file8, "\n");	// Leave blank line for Gnuplot before next j
		fprintf(file11, "\n");
		fprintf(file10, "\n");
	}

	fclose(file8);
	fclose(file10);
	fclose(file11);
}

//**********************************************************************
void movieY(int itime)	// Outputs files
{
	Float x, y, z, phpsi, depsi;
	char ch[17], cs[17], cu[17];
	int i, j, k;
	FILE *file8, *file10, *file11;
	string filepath;

	sprintf(ch, "densY%07d.dat", itime);
	filepath = path + ch;
	file8 = fopen(filepath.c_str(), "w");
	sprintf(cs, "phasY%07d.dat", itime);
	filepath = path + cs;
	file10 = fopen(filepath.c_str(), "w");
	sprintf(cu, "gravY%07d.dat", itime);
	filepath = path + cu;
	file11 = fopen(filepath.c_str(), "w");

	// Output files for contour plots

	for (j = 0; j <= Ny; j += 2)
		for (i = 0; i <= Nx; i += 2)
	{
		for (k = 0; k <= Nz; k += 2)
		{
			x = xgrid[i];
			y = ygrid[j];
			z = zgrid[k];
			depsi = SQ(real(psi(i, j, k))) + SQ(imag(psi(i, j, k)));
			phpsi = atan2(imag(psi(i, j, k)) + 1e-15f, real(psi(i, j, k)));
			fprintf(file8, "%lf %lf %lf %lg\n", x, y, z, depsi);
			fprintf(file11, "%lf %lf %lf %lg\n", x, y, z, phi(i, j, k));
			fprintf(file10, "%lf %lf %lf %lg\n", x, y, z, phpsi);
		}
		fprintf(file8, "\n");	// Leave blank line for Gnuplot before next j
		fprintf(file11, "\n");
		fprintf(file10, "\n");
	}

	fclose(file8);
	fclose(file10);
	fclose(file11);
}

//**********************************************************************
void movie(int itime)
{
	// movieZ(itime);
	movieX(itime);
	// movieY(itime);
}

//**********************************************************************
#ifdef USECL
inline void get_phi_kernel(bool bDir, Float G4pi, Float idx2, Float idy2,
						   Float idz2, Float h2)
{
	int a = 0;

	clSetKernelArg(kernelGet_res, a++, sizeof(cl_mem), &cldensity);
	clSetKernelArg(kernelGet_res, a++, sizeof(cl_mem), bDir ? &clphi : &clres);
	clSetKernelArg(kernelGet_res, a++, sizeof(cl_mem), bDir ? &clres : &clphi);
	clSetKernelArg(kernelGet_res, a++, sizeof(Float), &G4pi);
	clSetKernelArg(kernelGet_res, a++, sizeof(Float), &idx2);
	clSetKernelArg(kernelGet_res, a++, sizeof(Float), &idy2);
	clSetKernelArg(kernelGet_res, a++, sizeof(Float), &idz2);
	clSetKernelArg(kernelGet_res, a++, sizeof(Float), &h2);

	clEnqueueWriteBuffer(command_queue, cldensity, CL_TRUE, 0, sizeof(dens), dens,
						 0, NULL, NULL);
	clEnqueueWriteBuffer(command_queue, bDir ? clphi : clres, CL_TRUE, 0,
			bDir ? sizeof(phi) : sizeof(res), bDir ? phi : res, 0, NULL, NULL);
	clEnqueueNDRangeKernel(command_queue, kernelGet_res, 1, NULL, &globalThreads,
						   &localThreads, 0, NULL, NULL);
	clEnqueueReadBuffer(command_queue, bDir ? clres : clphi, CL_TRUE, 0,
			bDir ? sizeof(res) : sizeof(phi), bDir ? res : phi, 0, NULL, NULL);
}
#endif

void get_phi()	// Grav. potential via Poisson's Eq.
{
	Float rtot1, rtot2;
	const Float h2 = (Float)0.5 / (idx2 + idy2 + idz2);
	const Float h2l = (Float)0.5 / (idxl2 + idyl2 + idzl2);
	int i, j, k;
	const Float G4pi = 4 * pi * G;
	int nrmp = 0;

	rtot1 = -1;
	rtot2 = 1;
	for (nrmp = 0; fabs(rtot2 - rtot1) > tolPSN * fabs(rtot1); nrmp++)
	{
		if (nrmp > NRMP)
		{
			fprintf(stderr, "Poisson not converging after %d iterations "
							"(last: %lg -> %lg)\n", nrmp, rtot2, rtot1);
			fflush(stderr);
			return;
		}
		rtot2 = rtot1;
		rtot1 = 0;
#ifndef USECL
		for (k = 1; k < Nz; k++)
			for (j = 1; j < Ny; j++)
				for (i = 1; i < Nx; i++)
				{
				if (Nxl<i && i<Nxl+Nxf && Nyl<j && j<Nyl+Nyf && Nzl<k && k<Nzl+Nzf)
				{
					res(i, j, k) = ((phi(i - 1, j, k) + phi(i + 1, j, k)) * idx2 +
									(phi(i, j - 1, k) + phi(i, j + 1, k)) * idy2 +
									(phi(i, j, k - 1) + phi(i, j, k + 1)) * idz2 -
									G4pi * density(i, j, k)) * h2;
				}
				else{
					res(i, j, k) = ((phi(i - 1, j, k) + phi(i + 1, j, k)) * idxl2 +
									(phi(i, j - 1, k) + phi(i, j + 1, k)) * idyl2 +
									(phi(i, j, k - 1) + phi(i, j, k + 1)) * idzl2 -
									G4pi * density(i, j, k)) * h2l;
				}

				}
#else
		get_phi_kernel(true, G4pi, idx2, idy2, idz2, h2);
		get_phi_kernel(false, G4pi, idx2, idy2, idz2, h2);
#endif
		for (k = 1; k < Nz; k++)
			for (j = 1; j < Ny; j++)
				for (i = 1; i < Nx; i++)
				{				
					if (Nxl<i && i<Nxl+Nxf && Nyl<j && j<Nyl+Nyf && Nzl<k && k<Nzl+Nzf)
		{
#ifndef USECL
			phi(i, j, k) = ((res(i - 1, j, k) + res(i + 1, j, k)) * idx2 +
							(res(i, j - 1, k) + res(i, j + 1, k)) * idy2 +
							(res(i, j, k - 1) + res(i, j, k + 1)) * idz2 -
							G4pi * density(i, j, k)) * h2;
#endif
		}
					else{
			phi(i, j, k) = ((res(i - 1, j, k) + res(i + 1, j, k)) * idxl2 +
							(res(i, j - 1, k) + res(i, j + 1, k)) * idyl2 +
							(res(i, j, k - 1) + res(i, j, k + 1)) * idzl2 -
							G4pi * density(i, j, k)) * h2l;
		}
				rtot1 += phi(i, j, k);
				}
		//cout << rtot2 << " " << rtot1 << " " << fabs(rtot1-rtot2) <<endl;
	}

#ifdef SHOW_LOOPS
	fprintf(stderr, "Poisson required %d interations.\n", nrmp);
	fflush(stderr);
#endif



	// Add static baryonic potential
	for (i = 0; i <= Nx; i++)
		for (j = 0; j <= Ny; j++)
			for (k = 0; k <= Nz; k++)
	{
		phiU(i, j, k) = phi(i, j, k) 
		#ifdef BARY
		+ phiBary(i, j, k)
		#endif
		#ifdef SHIFTP
		+ shiftphi
		#endif
		;
	}


}

//**********************************************************************
void get_density()
{
	int i, j, k;

	for (i = 0; i <= Nx; i++)
		for (j = 0; j <= Ny; j++)
			for (k = 0; k <= Nz; k++)
				density(i, j, k) = SQ(real(psi(i, j, k))) + SQ(imag(psi(i, j, k)));
}

void readinterp(string file){
    Float *f_r = new Float[Nf];
    Float *f_psi = new Float[Nf];
    ifstream ifs(file, ios::in); // opening the file
    if (!ifs.is_open())
    {
        cout << "open interp file fail!" << endl;
    } 
    else
    {
        cout << "open interp file successful!" << endl;
        for (int i = 0; i < Nf; i++)
        {
            ifs >> f_r[i] >> f_psi[i];
        }
        ifs.close();
        cout << "Finished reading interp file! Number of entries: " << Nf << endl;
    }
    for (int i = 0; i <= Nf; i++)
	{
        r_interp[i] = f_r[i];
        psi_interp[i] = f_psi[i];
    }
}

template<typename T>
int nearestNeighbourIndex(std::vector<T>& x, T& value)
{
    T dist = std::numeric_limits<T>::max();
    T newDist = dist;
    size_t idx = 0;

    for (size_t i = 0; i < x.size(); ++i) {
        newDist = std::abs(value - x[i]);
        if (newDist <= dist) {
            dist = newDist;
            idx = i;
        }
    }

    return idx;
}


template<typename T>
std::vector<T> interp1(std::vector<T>& x, std::vector<T>& y, std::vector<T>& x_new)
{
    std::vector<T> y_new;
    T dx, dy, m, b;
    size_t x_max_idx = x.size() - 1;
    size_t x_new_size = x_new.size();

    y_new.reserve(x_new_size);

    for (size_t i = 0; i < x_new_size; ++i)
    {
        size_t idx = nearestNeighbourIndex(x, x_new[i]);

        if (x[idx] > x_new[i])
        {
            dx = idx > 0 ? (x[idx] - x[idx - 1]) : (x[idx + 1] - x[idx]);
            dy = idx > 0 ? (y[idx] - y[idx - 1]) : (y[idx + 1] - y[idx]);
        }
        else
        {
            dx = idx < x_max_idx ? (x[idx + 1] - x[idx]) : (x[idx] - x[idx - 1]);
            dy = idx < x_max_idx ? (y[idx + 1] - y[idx]) : (y[idx] - y[idx - 1]);
        }

        m = dy / dx;
        b = y[idx] - x[idx] * m;

        y_new.push_back(x_new[i] * m + b);
    }

    return y_new;
}
