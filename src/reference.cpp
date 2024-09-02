#include <vector>
#include "reference.h"


using namespace std;
typedef std::complex<double> cplx;
// 
void monopole(double Ap, vec3 x, vec3 y, double M0, double C0, double rho0, double omega, double t, double& p, vec3& u)
{
	double pi = 2 * acos(0.0);

	// square of beta
	double betaSq = 1.0 - M0 * M0;

	// amplitude distance
	double rp = sqrt(pow(x[0] - y[0], 2) + betaSq * (pow(x[1] - y[1], 2) + pow(x[2] - y[2], 2)));
	// phase distance
	double rq = (rp - M0 * (x[0] - y[0])) / betaSq;

	//! derivatives of the two radius
	vec3 rpd, rqd;
	rpd[0] = (x[0] - y[0]) / rp;
	rpd[1] = betaSq * (x[1] - y[1]) / rp;
	rpd[2] = betaSq * (x[2] - y[2]) / rp;
	rqd[0] = (rpd[0] - M0) / betaSq;
	rqd[1] = (rpd[1] - 0.0) / betaSq;
	rqd[2] = (rpd[2] - 0.0) / betaSq;

	//! calculating the results
	cplx cj(0., 1.);
	cplx ftemp = Ap * exp(cj * omega * t - cj * omega * rq) / (4.0 * rp * pi);   //速度势，声速为1

	//! the complex velocities
	vector<cplx> utemp;
	utemp.resize(3);
	utemp[0] = -(rpd[0] / rp + cj * omega * rqd[0]) * ftemp;//对x1求导
	utemp[1] = -(rpd[1] / rp + cj * omega * rqd[1]) * ftemp;//对x2求导
	utemp[2] = -(rpd[2] / rp + cj * omega * rqd[2]) * ftemp;
	u[0] = utemp[0].real();
	u[1] = utemp[1].real();
	u[2] = utemp[2].real();

	//! the pressure, and its gradient
	cplx ptemp = -(cj * omega * ftemp + M0 * utemp[0]);//自由流密度为1
	p = ptemp.real();
}


void dipole(double Ap, vec3 x, vec3 y, double M0, double C0, double rho0, double omega, double t, double& p, vec3& u) {

	double dh = 5E-5;

	vec3 yp = y + vec3(0, dh, 0);
	vec3 ym = y - vec3(0, dh, 0);

	double pp, pm;
	vec3   up, um;
	monopole( Ap, x, yp, M0, C0, rho0, omega, t,  pp, up);
	monopole( Ap, x, ym, M0, C0, rho0, omega, t,  pm, um);

	p = (pp - pm) / dh / 2.0;//中心差分代替求导；
	u = (up - um) / dh / 2.0;
}

void monopole2D(double Ap, vec3 x, vec3 y, double M0, double C0, double rho0, double omega, double t, double& p, vec3& u)
{
	// square of beta
	double betaSqr = 1.0 - M0 * M0;

	cplx cj(0., 1.);

	cplx phif = 0.25 * Ap * cj / sqrt(betaSqr) * exp(cj * omega * t + cj * M0 * omega / C0 * (x[0]-y[0]) / betaSqr);

	cplx hk1 = Hkel1(omega, C0, x[0], x[1], betaSqr);

	cplx phit = phif * hk1 * cj * omega;

	double h = 1E-7;
	cplx hk2 = Hkel1(omega, C0, x[0] + h, x[1], betaSqr);
	cplx hk3 = Hkel1(omega, C0, x[0] - h, x[1], betaSqr);
	cplx hkx = (hk2 - hk3) / 2.0 / h;

	cplx hk4 = Hkel1(omega, C0, x[0] , x[1]+h, betaSqr);
	cplx hk5 = Hkel1(omega, C0, x[0] , x[1]-h, betaSqr);
	cplx hky = (hk4 - hk5) / 2.0 / h;

	cplx phix =  ((cj * M0 * omega / C0 / betaSqr) * hk1 + hkx) * phif;
	cplx phiy =  phif * hky;
	u[0] = phix.real();
	u[1] = phiy.real();
	u[2] = 0.0;

	//! the pressure, 
	cplx pref = -rho0 * (phit + M0 * C0 * phix);
	p = pref.real();


}



std::complex<double> Hkel1(double omega, double C0, double x1, double x2,double beta2)
{
	cplx cj(0., 1.);
	//汉克尔函数参数
	double hk = omega / C0 / beta2 * sqrt(pow(x1, 2) + beta2 * pow(x2, 2));

	double H1 = cyl_bessel_j(0, hk);
	double H2 = cyl_neumann(0, hk);
	complex<double> hkel = H1 - cj * H2;
	return hkel;

}

void dipole2D(double Ap, vec3 x, vec3 y, double M0, double C0, double rho0, double omega, double t, double& p, vec3& u) {

	double dh = 5E-5;

	vec3 xp = x + vec3(0, dh, 0);
	vec3 xm = x - vec3(0, dh, 0);

	double pp, pm;
	vec3   up, um;
	monopole2D(Ap, xp, y, M0, C0, rho0, omega, t, pp, up);
	monopole2D(Ap, xm, y, M0, C0, rho0, omega, t, pm, um);

	p = (pp - pm) / dh / 2.0;//中心差分代替求导；
	u = (up - um) / dh / 2.0;
}


