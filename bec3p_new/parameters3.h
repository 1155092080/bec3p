// BEC3P.cpp : Self-gravitating Bose-Einstein condensate
//
// This version dated 2013/08/03.

#pragma once

// Set this to enable OPENCL (experimental;
//                            GPU needed, negligible gain in performance)
//#define USECL

#ifdef USECL
#define Float float
#else
#define Float double
#endif
#ifndef KERNEL
const Float pi = 4 * atan((Float)1);
#endif

#define SQ(x) ((x)*(x))

// User-configurable parameters begin here

// Self-gravitating or trapped condensate?
#define GRAV
#define BARY
#define INIFILE
// Real-time or imaginary-time evolution?
//#define REALTIME
// Use an isothermal wall potential?
#define ISOWALL
// Grid size
#define Nx 300//120
#define Ny 300//120
#define Nz 300//120

std::string prefix = "./no_bary/omega0/";
std:: string inifile = "./no_bary/omega0/psi_phi.dat";


const Float xl = -25.0f, yl = -25.0f, zl = -25.0f;
const Float xr = 25.0f, yr = 25.0f, zr = 25.0f;

#ifndef KERNEL
// Simulation parameters
const Float tau = 0.01;//10;					// Time step (units of [T])
const int time_n = 100;//0000;				// Number of iterations to run
const Float G = 0.0001;//0.0667;				// Newton's constant (may be scaled)
const Float N = 100000.0;//2.0;					// Particle number (may be scaled)
const Float a = 0.001;//0.5 * G * SQ(R/pi);		// Scattering length (TF default)
const Float c = 4 * pi * a;				// BEC interaction coupling strength
const Float R = pi*sqrt(a/G);//50.0;					// Size of initial condensate (in [L])
const Float rho0 = pi*N/(4*R*R*R);			// Initial density (in [L^-3])
const Float ex = 0.0;					// Softening parameters
const Float ey = 0.0;
const Float ez = 0.0;                 // In this way, the trap frequency along z is 3 of that of x and y direction
const Float r0 = 12.0;					// Isowall radius (in [L])
const Float isowallpot = 100.0;       // Isothermal wall potential
const Float omega0 = 0.0;//0.0001;				// Initial angular velocity (in rad/[T])
const Float gamma0 = 0.0;				// Softening parameter
//const int despin_n = 1;				// When to stop spinning the condensate
// const Float aho = 1.0;                // harmonic length (in [L])
const Float omg = 1.587;                  // harmonic trap (in rad/[T])
// Iteration tolerances
const Float tolGPE = 1e-6;				// GPE nonlinear term iteration
const Float tolPSN = 1e-4;				// Poisson relaxation method iteration
const Float tolREL = 1e-8;				// Imaginary time system relaxation

// Output control
const int nstep0 = 1;	// number of steps of initial transient without output
const int nstep1 = 10;	// every how many steps energy is output
const int nstep2 = 100;	// every how many steps contour plot is output
#endif
