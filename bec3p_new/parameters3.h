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
#define ISOTHERMAL
// Grid size
#define Nx 228//120
#define Ny 228//120
#define Nz 228//120

std::string prefix = "./small_test/fine15MilkyRot0.000/";
std:: string inifile = "./ini_vort/fine15MilkyRot0.112/psi_ini.dat";


const Float xl = -15.0f, yl = -15.0f, zl = -15.0f;
const Float xr = 15.0f, yr = 15.0f, zr = 15.0f;

#ifndef KERNEL
// Simulation parameters
const Float tau = 0.005;//10;					// Time step (units of [T])
const int time_n = 100;//0000;				// Number of iterations to run
const Float G = 0.0001;//0.0667;				// Newton's constant (may be scaled)
const Float N = 100000.0;//2.0;					// Particle number (may be scaled)
const Float a = 0.001;//0.5 * G * SQ(R/pi);		// Scattering length (TF default)
const Float c = 4 * pi * a;				// BEC interaction coupling strength
const Float R = pi*sqrt(a/G);//50.0;					// Size of initial condensate (in [L])
const Float rho0 = pi*N/(4*R*R*R);			// Initial density (in [L^-3])
const Float Chi = 50;                   // Dimensionless parameter for contraction
const Float kbT = rho0/Chi*c;           // Temperature    
const Float ex = 0.0;					// Softening parameters
const Float ey = 0.0;
const Float ez = 0.0;                 // In this way, the trap frequency along z is 3 of that of x and y direction
const Float omega0 = 0.00;//0.0001;				// Initial angular velocity (in rad/[T])
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
