/*
edited by: markus
     date: 2015-12-12
*/

#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;

//------------------------------------------------------------------------
void vector(double* k, const double mu, const double r1, const double r2, const double r3, const double r4){
	k[0] = r3;
	k[1] = r4;
	k[2] = r1+2.0*r4-(1.0-r2)*(r1+r2)/sqrt(pow(pow(r1+mu,2.0)+pow(r2,2.0),2))-mu*(r1-1.0+mu)/sqrt(pow(pow(r1-1.0+mu,2.0)+pow(r2,2.0),3));
	k[3] = r2+2.0*r3-(1.0-mu)*r3/sqrt(pow(pow(r1+mu,2.0)+pow(r2,2.0),3))-mu*r2/sqrt(pow(pow(r1-1.0+mu,2.0)+pow(r2,2.0),3));
}
//------------------------------------------------------------------------
void RKstep(const double* r, const double mu, const double dt, double* k1, double* k2, double* k3, double* k4, double* k5, double* k6, double* k7){
	double r1=r[0], r2=r[1], r3=r[2], r4=r[3];

	const double a21=1.0/5.0;
	const double a31=3.0/40.0, a32=9.0/40.0;
	const double a41=44.0/45.0, a42=-56.0/15.0, a43=32.0/9.0;
	const double a51=19372.0/6561.0, a52=-25360.0/2187.0, a53=64448.0/6561.0, a54=-212.0/729.0;
	const double a61=9017.0/3168.0, a62=-355.0/33.0, a63=46732.0/5247.0, a64=49.0/176.0, a65=-5103.0/18656.0;
	const double a71=35.0/384.0, a72=0, a73=500.0/1113.0, a74=125.0/192.0, a75=-2187.0/6784.0, a76=11.0/84.0;

	vector(k1, mu, r1, r2, r3, r4);
	vector(k2, mu, r1+dt*(a21*k1[0]), r2+dt*(a21*k1[1]), r3+dt*(a21*k1[2]), r4+dt*(a21*k1[3]));
	vector(k3, mu, r1+dt*(a31*k1[0]+a32*k2[0]), r2+dt*(a31*k1[1]+a32*k2[1]), r3+dt*(a31*k1[2]+a32*k2[2]), r4+dt*(a31*k1[3]+a32*k2[3]));
	vector(k4, mu, r1+dt*(a41*k1[0]+a42*k2[0]+a43*k3[0]), r2+dt*(a41*k1[1]+a42*k2[1]+a43*k3[1]), r3+dt*(a41*k1[2]+a42*k2[2]+a43*k3[2]), r4+dt*(a41*k1[3]+a42*k2[3]+a43*k3[3]));
	vector(k5, mu, r1+dt*(a51*k1[0]+a52*k2[0]+a53*k3[0]+a54*k4[0]), r2+dt*(a51*k1[1]+a52*k2[1]+a53*k3[1]+a54*k4[1]), r3+dt*(a51*k1[2]+a52*k2[2]+a53*k3[2]+a54*k4[2]), r4+dt*(a51*k1[3]+a52*k2[3]+a53*k3[3]+a54*k4[3]));
	vector(k6, mu, r1+dt*(a61*k1[0]+a62*k2[0]+a63*k3[0]+a64*k4[0]+a65*k5[0]), r2+dt*(a61*k1[1]+a62*k2[1]+a63*k3[1]+a64*k4[1]+a65*k5[1]), r3+dt*(a61*k1[2]+a62*k2[2]+a63*k3[2]+a64*k4[2]+a65*k5[2]), r4+dt*(a61*k1[3]+a62*k2[3]+a63*k3[3]+a64*k4[3]+a65*k5[3]));
	vector(k7, mu, r1+dt*(a71*k1[0]+a72*k2[0]+a73*k3[0]+a74*k4[0]+a75*k5[0]+a76*k6[0]), r2+dt*(a71*k1[1]+a72*k2[1]+a73*k3[1]+a74*k4[1]+a75*k5[1]+a76*k6[1]), r3+dt*(a71*k1[2]+a72*k2[2]+a73*k3[2]+a74*k4[2]+a75*k5[2]+a76*k6[2]), r4+dt*(a71*k1[3]+a72*k2[3]+a73*k3[3]+a74*k4[3]+a75*k5[3]+a76*k6[3]));
}
//------------------------------------------------------------------------
void RK4(double* r, const double mu, double dt, double* k1, double* k2, double* k3, double* k4, double* k5, double* k6, double* k7){
	const double b14=5179.0/57300.0, b24=0, b34=7571.0/16695.0, b44=393.0/640.0, b54=-92097.0/339200.0, b64=187.0/2100.0, b74=1.0/40.0;
	RKstep(r, mu, dt, k1, k2, k3, k4 , k5, k6, k7);
	for(int i=0; i<4; i++) r[i]+=dt*(b14*k1[i]+b24*k2[i]+b34*k3[i]+b44*k4[i]+b54*k5[i]+b64*k6[i]+b74*k7[i]);
	
}
//------------------------------------------------------------------------
void RK5(double* s, const double mu, double dt, double* k1, double* k2, double* k3, double* k4, double* k5, double* k6, double* k7){
	const double b15=35.0/384.0, b25=0, b35=500.0/1113.0, b45=125.0/192.0, b55=-2187.0/6784.0, b65=11.0/84.0;
	RKstep(s, mu, dt, k1, k2, k3, k4 , k5, k6, k7);
	for(int i=0; i<4; i++) s[i]+=dt*(b15*k1[i]+b25*k2[i]+b35*k3[i]+b45*k4[i]+b55*k5[i]+b65*k6[i]);
}
//------------------------------------------------------------------------
void ssc(const double* r, const double* s, double& maxi){
	double norm[4];
	for(int i=0; i<4; i++) norm[i]=abs(r[i]-s[i]);
	double max1=max(norm[0],norm[1]);
	double max2=max(norm[2],norm[3]);
	maxi=max(max1,max2);
}
//------------------------------------------------------------------------
int main(){

	ofstream out("solution");
	const int dim = 4;
	const double mu = 0.012277471;
	double dt = 1e-3;
	double T = 17.065216560157;
	double t = 0;
	double k1[dim],k2[dim],k3[dim],k4[dim],k5[dim],k6[dim],k7[dim];
	double r[dim]={0.994, 0, 0, -2.00158510637908};
	double s[dim]={0.994, 0, 0, -2.00158510637908};
	double u[dim],v[dim];
	double maxi;
	double tol=1e-3;

	out << t << "\t" << r[0] << "\t" << r[1] << "\t" << s[0] << "\t" << s[1] << endl;			// printing initial values

	while(t<=T){

		for(int i=0; i<4; i++){u[i]=r[i];v[i]=s[i];}

		RK4(r, mu, dt, k1, k2, k3, k4, k5, k6, k7);

		RK5(s, mu, dt, k1, k2, k3, k4, k5, k6, k7);

		ssc(r, s, maxi);

		dt*=pow(tol/maxi, 0.2);

		for(int i=0; i<4; i++){r[i]=u[i];s[i]=v[i];}

		RK4(r, mu, dt, k1, k2, k3, k4, k5, k6, k7);

		RK5(s, mu, dt, k1, k2, k3, k4, k5, k6, k7);

		t+=dt;


	out << t << "\t" << r[0] << "\t" << r[1] << "\t" << s[0] << "\t" << s[1] << endl;	
	}

	out.close();
	return 0;
}
