#include "stdafx.h"
#include <string>

#include "parameters3.h"
// #include "phi-interpolator.h"
#include "mkpath.h"
#ifdef USECL
#include <CL/opencl.h>
#include <SDKCommon.hpp>
#include <SDKFile.hpp>
#endif

#pragma warning(disable:4996)

using namespace std;
using std::string;
const int Nn = (Nx + 1) * (Ny + 1) * (Nz + 1); // Number of total grid points

#define NRMN 100
#define NRMC 1000
#define NRMP 10000
#define NTOL 1e-6

complex<Float> eye, dt;
Float t, dx, dy, dz, idx2, idy2, idz2;
string path = prefix;

// Here flatten all the 3D data into 1D by ijk(i,j,k) function
complex<Float> *psi = new complex<Float>[Nn];  
complex<Float> *psi_n = new complex<Float>[Nn];
Float *dens = new Float[Nn];
Float *phi = new Float[Nn];
Float *phiU = new Float[Nn];
Float *UU = new Float[Nn];
Float *res = new Float[Nn];

// Here map all the 1D data back to 3D by spicifying the counting method
#define ijk(i,j,k) ((((i) * (Ny + 1)) + (j)) * (Nz + 1) + (k)) // Start from 0, increase by z then y then x. So (1,0,0) is the (Nz+1)*(Ny+1)-th point
#define psi(i,j,k) (psi[ijk(i,j,k)])
#define psi_n(i,j,k) (psi_n[ijk(i,j,k)])
#define density(i,j,k) (dens[ijk(i,j,k)])
#define phi(i,j,k) (phi[ijk(i,j,k)])
#define phiU(i,j,k) (phiU[ijk(i,j,k)])
#define U(i,j,k) (UU[ijk(i,j,k)])
#define res(i,j,k) (res[ijk(i,j,k)])


// Forward declarations
Float init(int, int, int);
void get_phi();
void get_density();
Float get_normsimp();

int _tmain(int argc, _TCHAR* argv[])
{
	complex<Float> xi;
	Float omega, gamma, mu;
	Float norm, norm0;
	Float E0, E1 = 0, E2 = 0;
	int i, j, k, itime, ktime;
    FILE *filephi, *filepsi;
	Float nrma[NRMN];
	int nrmc;
	double norm_ini;
	int mkdirretval;
    //mkdirretval=light::mkpath("foo2/bar",0755);
    //mkdirretval=light::mkpath("./lsl/foo2/bar");
    mkdirretval=light::mkpath(path.c_str());
    std::cout << mkdirretval << '\n';
    string filepath;

#ifdef USECL
	initializeCL();
#endif

	// Fixed parameters

	eye = complex<Float>(0, 1);
	dx = (xr - xl) / Nx;
	dy = (yr - yl) / Ny;
	dz = (zr - zl) / Nz;

	idx2 = 1 / dx / dx;
	idy2 = 1 / dy / dy;
	idz2 = 1 / dz / dz;

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

	dt = complex<Float>(0, -tau);
	omega = 0;
	gamma = 0;
#ifdef GRAV
	mu = 0;
#else
	mu = 0;//0.5 / SQ(R);
#endif



	// Initial state
    filepath = path + "psi_ini.dat";
	filepsi = fopen(filepath.c_str(), "w");
	t = 0.0;
	itime = 0;
	ktime = 0;
	for (i = 0; i <= Nx; i++)
	{
		for (j = 0; j <= Ny; j++)
			for (k = 0; k <= Nz; k++)
		{
#ifdef GRAV
			Float rho = init(i, j, k);
			Float x2 = SQ(xl + i * dx);
			Float y2 = SQ(yl + j * dy);
			Float z2 = SQ(zl + k * dz);
			Float r = sqrt(x2 + y2 + z2);

			psi(i, j, k) = sqrt(rho);
			phi(i, j, k) = (Float)(-G * N / (r > (.25 * dx) ? r : .5 * dx));
#else
			psi(i, j, k) = fermi(mu, i, j, k);
#endif
        }
        fprintf(filepsi, "\n");	// For Gnuplot
    }
    fclose(filepsi);
    norm_ini = get_normsimp();
	printf("Initial norm is P=%11.4lg\n", norm_ini);
	fflush(stdout);
    printf("Setting up initial density...\n");
    fflush(stdout);

	get_density();

printf("Setting up initial potential...\n");
fflush(stdout);

#ifdef GRAV
	get_phi();
#else
	get_Vtr();
#endif
// Save trap information
filepath = path + "phi_ini.dat";
filephi = fopen(filepath.c_str(), "w");
for (i = 0; i <= Nx; i++)
{
    for (j = 0; j <= Ny; j++)
        for (k = 0; k <= Nz; k++)
    {
            fprintf(filephi, "%lg %lg %lg %lg\n", xl + i * dx, yl + j * dy,
                                        zl + k * dz, phi(i, j, k));
    }
    fprintf(filephi, "\n");	// For Gnuplot
}
fclose(filephi);
}



Float init(int i, int j, int k)		
{
	Float F, x, y, z, r;

	x = xl + i * dx;
	y = yl + j * dy;
	z = zl + k * dz;
	r = sqrt((1 + ex) * x * x + (1 + ey) * y * y + (1 + ez) * z * z);
	F = 0.0;
	if (r > R) F = 0.0;
	else F = (Float)N  / (4 * pi * R * R * R / 3);
	return F;
}

void get_phi()	// Grav. potential via Poisson's Eq.
{
	Float rtot1, rtot2;
	const Float h2 = (Float)0.5 / (idx2 + idy2 + idz2);
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
			res(i, j, k) = ((phi(i - 1, j, k) + phi(i + 1, j, k)) * idx2 +
							(phi(i, j - 1, k) + phi(i, j + 1, k)) * idy2 +
							(phi(i, j, k - 1) + phi(i, j, k + 1)) * idz2 -
							G4pi * density(i, j, k)) * h2;
		}
#else
		get_phi_kernel(true, G4pi, idx2, idy2, idz2, h2);
		get_phi_kernel(false, G4pi, idx2, idy2, idz2, h2);
#endif
		for (k = 1; k < Nz; k++)
			for (j = 1; j < Ny; j++)
				for (i = 1; i < Nx; i++)
		{
#ifndef USECL
			phi(i, j, k) = ((res(i - 1, j, k) + res(i + 1, j, k)) * idx2 +
							(res(i, j - 1, k) + res(i, j + 1, k)) * idy2 +
							(res(i, j, k - 1) + res(i, j, k + 1)) * idz2 -
							G4pi * density(i, j, k)) * h2;
#endif
			rtot1 += phi(i, j, k);
		}
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
		double x = xl + i * dx;
		double y = yl + j * dy;
		double z = zl + k * dz;
		phiU(i, j, k) = phi(i, j, k) 
		#ifdef BARY
		+ Ub(x, y, z)
		#endif
		;
	}


}

void get_density()
{
	int i, j, k;

	for (i = 0; i <= Nx; i++)
		for (j = 0; j <= Ny; j++)
			for (k = 0; k <= Nz; k++)
				density(i, j, k) = SQ(real(psi(i, j, k))) + SQ(imag(psi(i, j, k)));
}

Float get_normsimp()	// Find the norm using Simpson's rule
{
	static Float normz[Nx + 1][Ny + 1], normy[Nx + 1], norm;
	int i, j, k;

	norm = 0.0;
	for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++)
	{
		normz[i][j] = SQ(real(psi(i, j, 0))) + SQ(imag(psi(i, j, 0))) +
					SQ(real(psi(i, j, Nz - 1))) + SQ(imag(psi(i, j, Nz - 1)));
		for (k = 1; k <= Nz - 3; k += 2)
			normz[i][j] += 4 * (SQ(real(psi(i, j, k))) +
								SQ(imag(psi(i, j, k)))) +
						   2 * (SQ(real(psi(i, j, k + 1))) +
								SQ(imag(psi(i, j, k + 1))));
	}
	for (i = 0; i < Nx; i++)
	{
		normy[i] = normz[i][0] + normz[i][Ny - 1];
		for (j = 1; j <= Ny - 3; j += 2)
			normy[i] += 4 * normz[i][j] + 2 * normz[i][j + 1];
	}
	norm = normy[0] + normy[Nx - 1];
	for (i = 1; i <= Nx - 3; i += 2)
		norm += 4 * normy[i] + 2 * normy[i + 1];
	return norm * dx * dy * dz / 27;
}