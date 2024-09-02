#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "reference.h"
#include "Config.h"
#include "vec3.h"
/*****************************************************************/
//
//	TODO: Isolate file write/
// 
//
/*****************************************************************/

using namespace std;

int main() {

	string configFile = "config.inp";  
	Config CF(configFile); 
	/*****************************************************************/
	// define freestream condition
	double MxInf = 0.0;
	double rhoInf = 0.0;
	double cInf = 0.0;
	double pInf = 0.0;
	MxInf = CF.Read("FreeStreamMachNumber", MxInf);  
	rhoInf = CF.Read("FreestreamDensity", rhoInf);
	cInf = CF.Read("FreestreamSoundSpeed", cInf);
	pInf = CF.Read("FreestreamPressure", pInf);
	cout << "Read freestream condition" << endl;
	cout << "Mach number in x+ direction " << MxInf << endl;
	cout << "Speed of sound " << cInf << endl;
	cout << "Density " << rhoInf << endl;
	cout << "Pressure " << pInf << endl;
	/*****************************************************************/
	
	/**************** define observers *******************************/
	// set observer points of interest
	vector<vec3> observers;
	int	obsNum = 30;
	double	obsR = 30.;
	const double pi = 2 * acos(0.0);
	for (int i = 0; i < obsNum; i++)
	{
		vec3 point;
		double theta = i * (pi / (obsNum - 1));
		point[0] = obsR * cos(theta);
		point[1] = obsR * sin(theta);
		point[2] = 0.0;
		observers.push_back(point);
	}

	
	string observerFileName;
	observerFileName = CF.Read("observerFileName", observerFileName);
	cout << "Write observer file: " << observerFileName << endl;
	ofstream observerFile(observerFileName.c_str(), ios::out);
	// time-step and cell number
	observerFile << obsNum << " " << endl;
	for (int io = 0; io < obsNum; io++) {
		observerFile << observers[io][0] << " " << observers[io][1] << " " << observers[io][2] << endl;
	}
	observerFile.close();
	/*****************************************************************/


	/**************** read flow variables ****************************/
	/*************** write integral surface **************************/
	// define integral surface using a sphere
	// TODO:
	//	this integral surface should be read from a file
	int nTau = 5000;
	double tauMax = 50;
	double tauMin = 0.;
	double dTau = (tauMax - tauMin) / nTau;
	vector<double> tau(nTau, 0.0);         
	for (int i = 0; i < nTau; i++) {
		tau[i] = double(i) * dTau + tauMin;
	}

	double R = 1.0;
	size_t nTheta = 20;
	size_t nPhi = 20;
	double dTheta = pi / nTheta;
	vector<double> theta(nTheta, 0.0);
	for (int it = 0; it < nTheta; it++)
		theta[it] = (it + 0.5) * (pi / nTheta);

	double dPhi = 2 * pi / nPhi;
	vector<double> phi(nPhi, 0.0);
	for (int ip = 0; ip < nPhi; ip++)
		phi[ip] = (ip + 0.5) * dPhi;

	size_t nPanels = nTheta * nPhi;
	vector<vector<vec3>> x;
	vector<vector<vec3>> n;
	vector<double> s;

	s.resize(nPanels, 0.0);
	for (size_t it = 0; it < nTau; it++)
	{
		x.push_back(vector<vec3>(nPanels, vec3(0., 0., 0.)));
		n.push_back(vector<vec3>(nPanels, vec3(0., 0., 0.)));
	}

	for (size_t it = 0; it < 1; it++)
	{
		for (size_t ab = 0; ab < nTheta; ab++)
		{
			for (size_t cd = 0; cd < nPhi; cd++)
			{
				size_t ip = ab * nPhi + cd;
				s[ip] = R * R * sin(theta[ab]) * dTheta * dPhi;
				x[it][ip][0] = R * sin(theta[ab]) * cos(phi[cd]);
				x[it][ip][1] = R * sin(theta[ab]) * sin(phi[cd]);
				x[it][ip][2] = R * cos(theta[ab]);
				n[it][ip][0] = R * sin(theta[ab]) * cos(phi[cd]);
				n[it][ip][1] = R * sin(theta[ab]) * sin(phi[cd]);
				n[it][ip][2] = R * cos(theta[ab]);
			}
		}
	}
	string intSurfaceFileName; 
	intSurfaceFileName = CF.Read("intSurfaceFileName", intSurfaceFileName);
	cout << "Write integral surface file: " << intSurfaceFileName << endl;
	
	// TDODO: deine output path and file name
	ofstream fwhSurfacefile(intSurfaceFileName.c_str(), ios::out);
	// time-step and cell number
	
	fwhSurfacefile << 1 << " " << nPanels << " " << endl;
	for (int ip = 0; ip < nPanels; ip++) {
		fwhSurfacefile << s[ip] << " ";
	}
	for (int it = 0; it < 1; it++) {
		for (int ip = 0; ip < nPanels; ip++) {
			fwhSurfacefile << x[it][ip][0] << " " << x[it][ip][1] << " " << x[it][ip][2] << " ";
			fwhSurfacefile << n[it][ip][0] << " " << n[it][ip][1] << " " << n[it][ip][2] << " ";
			fwhSurfacefile << endl;
		}
	}
	fwhSurfacefile.close();
	/********************************************************************/


	/*************** write 2D integral surface **************************/
	// define integral surface using a circle
	double Rr = 1.0;
	size_t nTheta1 = 360;
	double dTheta1 = 2 * pi / nTheta1;
	vector<double> theta1(nTheta1, 0.0);
	for (int it = 0; it < nTheta1; it++)
		theta1[it] = (it + 0.5) * dTheta1;

	vector<vec3> x1;
	vector<vec3> n1;
	double darea1 = Rr * dTheta1;
	x1.resize(nTheta1, vec3(0., 0., 0.));
	n1.resize(nTheta1, vec3(0., 0., 0.));

	for (size_t it = 0; it < nTheta1; it++)
	{
		x1[it][0] = Rr * cos(theta1[it]) ;
		x1[it][1] = Rr * sin(theta1[it]) ;
		x1[it][2] = 0.0 ;
		n1[it][0] = Rr * cos(theta1[it]) ;
		n1[it][1] = Rr * sin(theta1[it]) ;
		n1[it][2] = 0.0;
	}

	string SurfaceFileName2D;
	SurfaceFileName2D = CF.Read("SurfaceFileName2D", SurfaceFileName2D);
	cout << "Write 2D integral surface file: " << SurfaceFileName2D << endl;
	// TDODO: deine output path and file name
	ofstream fwh2DSurfacefile(SurfaceFileName2D.c_str(), ios::out);

	
	fwh2DSurfacefile << 1 << " " << nTheta1 << " " << endl;
	
	fwh2DSurfacefile << darea1 << " " << endl;

		for (int ip = 0; ip < nTheta1; ip++) {
			fwh2DSurfacefile << x1[ip][0] << " " << x1[ip][1] << " " << x1[ip][2] << " ";
			fwh2DSurfacefile << n1[ip][0] << " " << n1[ip][1] << " " << n1[ip][2] << " ";
			fwh2DSurfacefile << endl;
		}
	
	fwh2DSurfacefile.close();
	/********************************************************************/


	/*************** write 3D CFD Flow variables **************************/
	
	double freq = 0.5;
	double omega = 2.0 * pi * freq;
	vec3 y(0., 0., 0.);
	// allocate memory space before read data

	string flowDataFileName;
	flowDataFileName = CF.Read("flowDataFileName", flowDataFileName);
	cout << "Write flow-solutions file: " << flowDataFileName << endl;
	ofstream fwhFlowfile(flowDataFileName.c_str());
	fwhFlowfile << nTau << " " << nPanels << " " << endl;
	// Note:
	// default precison is 6 figures, which may lead to  truncation errors
	/********************************************************************/
	fwhFlowfile.precision(12);
	
	vec3 u;
	double p, rho;
	cout << "pInf" << pInf << endl;
	for (int it = 0; it < nTau; it++)
	{
		double& t = tau[it];
		printf("\r");
		cout << "3D Current progress: " << floor(100 * (it + 1) / nTau) << " % ";
		fwhFlowfile << t << " ";
		for (int ip = 0; ip < nPanels; ip++)
		{
			//monopole(1.0 / 340., x[0][ip], y, MxInf, cInf, rhoInf, omega, t, p, u);  //monopole
			dipole(1.0 / 340., x[0][ip], y, MxInf, cInf, rhoInf, omega, t, p, u);  //dipole
			rho = p / cInf / cInf + rhoInf;
			p += pInf;
			fwhFlowfile << rho << " " << p << " " << u[0] +MxInf*cInf<< " " << u[1] << " " << u[2] << " ";
			//fwhFlowfile << 1.0 << " " << 1.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " ";
		}
		fwhFlowfile << endl;
	}
	fwhFlowfile.close();
	printf("\n");
	/********************************************************************/
	

	/*************** write 2D CFD Flow variables **************************/
	string flowDataFileName2D;
	flowDataFileName2D = CF.Read("flowDataFileName2D", flowDataFileName2D);
	cout << "Write 2D flow-solutions file: " << flowDataFileName2D << endl;
	ofstream fwhFlowfile2D(flowDataFileName2D.c_str(), ios::out);
	fwhFlowfile2D << nTau << " " << nTheta1 << " " << endl;
	// Note:
	// default precison is 6 figures, which may lead to  truncation errors
	/********************************************************************/
	fwhFlowfile2D.precision(12);

	
	vec3 u1;
	double p1, rho1;
	for (int it = 0; it < nTau; it++)
	{
		double& t1 = tau[it];
		printf("\r");
		cout << "2D Current progress: " << floor(100 * (it + 1) / nTau) << " % ";
		fwhFlowfile2D << t1 << " ";
		for (int ip = 0; ip < nTheta1; ip++)//360
		{
			//monopole2D(1.0 / 340., x1[ip], y, MxInf, cInf, rhoInf, omega, t1, p1, u1);  
			dipole2D(1.0 / 340., x1[ip], y, MxInf, cInf, rhoInf, omega, t1, p1, u1);  
			rho1 = p1 / cInf / cInf + rhoInf;
			p1 += pInf;
			fwhFlowfile2D << rho1 << " " << p1 << " " << u1[0] + MxInf * cInf << " " << u1[1] << " " << u1[2] << " ";
			//fwhFlowfile2D << 1.0 << " " << 1.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " ";
		}
		fwhFlowfile2D << endl;
	}
	fwhFlowfile2D.close();
	printf("\n");
	/********************************************************************/

	double timeDelayMin;
	timeDelayMin = 16.0;
	
	/*************** write panalytic solutions **************************/
	vector<double> obsTime(tau);
	for (int it = 0; it < obsTime.size(); it++)
		obsTime[it] += timeDelayMin; // add time shift

	vector<vector<double>> pAnalytic;
	pAnalytic.resize(nTau, vector<double>(obsNum, 0.));//

	string refSolutionFileName = "pAnalytic.dat";
	cout << "Write analytic solution file: " << refSolutionFileName << endl;

	ofstream pRefFile(refSolutionFileName.c_str(), ios::out);
	for (int it = 0; it < obsTime.size(); it++)//
	{
		double& t = obsTime[it];
		double& t2 = tau[it];
		pRefFile << t << " ";
		for (int io = 0; io < obsNum; io++)//
		{
			//****3D*****/
			//monopole(1.0 / 340., observers[io], y, MxInf, cInf, rhoInf, omega, t, p, u);//A
			dipole(1.0 / 340., observers[io], y, MxInf, cInf, rhoInf, omega, t, p, u);//A
			
			//****2D*****/
			//monopole2D(1.0 / 340., observers[io], y, MxInf, cInf, rhoInf, omega, t2, p, u);//A
			//dipole2D(1.0 / 340., observers[io], y, MxInf, cInf, rhoInf, omega, t, p, u);

			pRefFile << p << " ";
		}
		pRefFile << endl;
	}
	pRefFile.close();
	return 0;
}