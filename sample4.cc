// table 1 (4 fields) of 1901.03714

#include<iostream>
#include<cmath>
#include"simplebounce.h"
using namespace std;

class model4 : public genericModel{
  public:
	double c0;
	double c1;
	double c2;
	double c3;
	double c4;
	model4(){
		c0 = 0.534808;
		c1 = 0.77023;
		c2 = 0.838912;
		c3 = 0.00517238;
		c4 = 0.258889;
		nphi = 4;
		dvdphi = new double[nphi];
	}
	~model4(){
		delete[] dvdphi;
	}
	// potential of scalar field(s)
	double vpot(const double* phi) const {
		double r1 = (
			c0*(phi[0]-1.)*(phi[0]-1.)
			+ c1*(phi[1]-1.)*(phi[1]-1.)
			+ c2*(phi[2]-1.)*(phi[2]-1.)
			+ c3*(phi[3]-1.)*(phi[3]-1.)
		);
		double r2 = (
			phi[0]*phi[0]
			+ phi[1]*phi[1]
			+ phi[2]*phi[2]
			+ phi[3]*phi[3]
		);
		return (r1-c4)*r2;
	}

	// derivative of potential of scalar field(s)
	void calcDvdphi(const double* phi) const {
		double r1 = (
			c0*(phi[0]-1.)*(phi[0]-1.)
			+ c1*(phi[1]-1.)*(phi[1]-1.)
			+ c2*(phi[2]-1.)*(phi[2]-1.)
			+ c3*(phi[3]-1.)*(phi[3]-1.)
		);
		double r2 = (
			phi[0]*phi[0]
			+ phi[1]*phi[1]
			+ phi[2]*phi[2]
			+ phi[3]*phi[3]
		);
		dvdphi[0] = 2.*c0*(phi[0]-1.)*r2 + 2.*phi[0]*(r1-c4);
		dvdphi[1] = 2.*c1*(phi[1]-1.)*r2 + 2.*phi[1]*(r1-c4);
		dvdphi[2] = 2.*c2*(phi[2]-1.)*r2 + 2.*phi[2]*(r1-c4);
		dvdphi[3] = 2.*c3*(phi[3]-1.)*r2 + 2.*phi[3]*(r1-c4);
	}
};






int main() {

	int n = 100; // number of grid
	double rmax = 1.; // number of space dimension 
	int dim = 3; // phi(rmax) = phi(False vacuum)

	bounce c(n, rmax, dim);
	model4 Model;
	c.setModel(&Model);

	double phiTV[4] = {1.,1.,1.,1.}; // a point at which V<0
	double phiFV[4] = {0.,0.,0.,0.}; // false vacuum
	c.setVacuum(phiTV, phiFV);

	// calcualte the bounce solution
	c.solve();

	// show the results
	c.printBounce();

	// Euclidean action
	cerr << c.action() << endl;

	return 0;
}