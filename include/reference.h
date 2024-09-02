#pragma once
#include "vec3.h"
#include <complex>
#include <cmath>


// declare analytic function
// this function gives pressure fluctuations produced by a monopole sound source
void monopole(
	double Ap, 
	vec3 x, 
	vec3 y, 
	double M0, 
	double C0, 
	double rho0, 
	double omega, 
	double t, 
	double& p, 
	vec3& u);

void dipole(
	double Ap,
	vec3 x,
	vec3 y,
	double M0,
	double C0,
	double rho0,
	double omega,
	double t,
	double& p,
	vec3& u);

void monopole2D(
	double Ap,
	vec3 x,
	vec3 y,
	double M0,
	double C0,
	double rho0,
	double omega,
	double t,
	double& p,
	vec3& u);

std::complex<double> Hkel1(
	double omega, 
	double C0, 
	double x1, 
	double x2, 
	double beta2);

void dipole2D(
	double Ap,
	vec3 x,
	vec3 y,
	double M0,
	double C0,
	double rho0,
	double omega,
	double t,
	double& p,
	vec3& u);

