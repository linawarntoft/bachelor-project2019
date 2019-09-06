// Author: Florent Renaud, Lund University (2018)

#ifndef MYGALAXY_HEAD
#define MYGALAXY_HEAD

#include <cmath>
#include "units.h"

//// User definition of the galactic potential (modeB only)

class MYGALAXY { 
  public:
    MYGALAXY(UNITS* uni) { // constructor: initialization function (called once, to compute constants)
    	u = uni;
    	init = true;
    	phi = potential(0,0,0,0);
    	init = false;
    }

    double potential(double x, double y, double z, double t) { // user defined potential
    	phi = 0.0; // reset galactic potential

    	/////////////////////////////////////////////////////////////////////
    	///// Put here quantities computed once per position and time    ////
    	/////////////////////////////////////////////////////////////////////
    	double x2 = x*x;
    	double y2 = y*y;
    	double z2 = z*z;
    	double r2d2 = x2+y2;
    	double r2 = r2d2+z2;
    	double r = sqrt(r2);
    	/////////////////////////////////////////////////////////////////////

    	/////////////////////////////////////////////////////////////////////
    	///// Add here components to the potential                       ////
    	///// (all functions used (even if not used at some positions    ////
    	///// or times) must be called here, to be properly initialized) ////
    	/////////////////////////////////////////////////////////////////////
    	pointmass(r);
    	//plummer(r2);
    	// hernquist(r);
    	// nfw(r);
    	// nfw2(r);
    	// isothermal(r2);
    	// nakedmndisk(r2d2, z2);
    	// bulgediskhalo(r2, r2d2, z2);
    	// spiralarms(x, y, z, r, t);
    	// nfwcosmo(r, t);
    	/////////////////////////////////////////////////////////////////////

      return phi;
    }

  private:
		UNITS* u;
		bool init;
        double phi;

/////////////////////////////////////////////////////////////////////
//// Potential components:
/////////////////////////////////////////////////////////////////////

    void pointmass(double r) {
    	static double mass;
    	if(init) {
    		mass = 1e10 / u->msun;
    	}
    	phi += -mass / r;
    }

//////////////////////////////////////////////////////////////////////////////////////////////////

    void plummer(double r2) {
    	static double mass, r02;
    	if(init) {
    		mass = 1e10 / u->msun;
    		r02 = 1000.0 / u->pc;
    		r02 = r02*r02;
    	}
    	phi += -mass / sqrt(r02 + r2);
    }

//////////////////////////////////////////////////////////////////////////////////////////////////

    void hernquist(double r) {
			static double mass, r0;
    	if(init) {
    		mass = 1e10 / u->msun;
    		r0 = 1000. / u->pc;
    	}
    	phi += -mass / sqrt(r0 + r);
    }

//////////////////////////////////////////////////////////////////////////////////////////////////

    void nfw(double r) {
      static double mass, r0;
    	if(init) {
    		mass = 1e10 / u->msun;
    		r0 = 1000. / u->pc;
    	}
    	phi += -mass / r * log(1.0 + r/r0);
    }

//////////////////////////////////////////////////////////////////////////////////////////////////

    void nfw2(double r) {
    	static double mvir, c, rs;
    	if(init) {
    		mvir = 23805e12 / u->msun; // virial mass
    		c = 15.8163; // concentration				
				double h0 = 71.0; // km/s/Mpc
    		h0 = h0 * 1e-6 / u->kms;
        rs = pow(mvir / ( (4./3.)*M_PI*200.0 * 3.*h0*h0/(8.*M_PI) ), 1./3.) / c;
    	}
    	phi += -mvir * log(1.0 + r/rs) / (log(1.0 + c) - c/ (1.0 + c));
    } 

//////////////////////////////////////////////////////////////////////////////////////////////////

    void isothermal(double r2) {
    	static double v02;
    	if(init) {
    		v02 = 220.0 / u->kms;
    		v02 = 0.5 * v02*v02;
    	}
    	phi += v02 * log(r2);
    }
     
//////////////////////////////////////////////////////////////////////////////////////////////////

    void nakedmndisk(double r2d2, double z2) {
    	static double mass, a, b2;
    	if(init) {
    		mass = 1e11 / u->kms;
    		a = 5000. / u->pc;
    		b2 = 300.0 / u->pc;
    		b2 = b2*b2;
    	}
    	phi += -mass /sqrt(r2d2+pow((a+sqrt(z2+b2)),2));
    } 

//////////////////////////////////////////////////////////////////////////////////////////////////

    void bulgediskhalo(double r2, double r2d2, double z2) {
    	static double mb[2], rb2[2], md[3], rd[3], bd2, vh2, rh2;
    	if(init) {
    		// bulge parameters
        mb[0] = 3.0e9 / u->msun;
        rb2[0] = 2700.0 / u->pc;
        rb2[1] = rb2[0]*rb2[0];
        mb[1] = 1.6e10 / u->msun;
        rb2[1] = 420.0 / u->pc;
        rb2[1] = rb2[1]*rb2[1];

    		// disk parameters
        md[0] = 8.9e10 / u->msun;
        md[1] = -6.9e10 / u->msun;
        md[2] = 2.8e10 / u->msun;
        rd[0] = 5000. / u->pc;
        rd[1] = 15800. / u->pc;
        rd[2] = 33000. / u->pc;
        bd2 = 300.0 / u->pc;
        bd2 = bd2*bd2;

  			// halo parameters
  			vh2 = 225.0 / u->kms;
    		vh2 = 0.5 * vh2*vh2;
    		rh2 = 8400. / u->pc;
    		rh2 = rh2*rh2;
    	}

    	phi += -mb[0] / sqrt(r2 + rb2[0]);
    	phi += -mb[1] / sqrt(r2 + rb2[1]);
    	for(int i=0; i<3; ++i) {
    		phi += -md[i] / sqrt(r2d2+pow(rd[i]+sqrt(z2+bd2),2));
    	 }
    	phi += vh2 * log(r2 + rh2);
    } 

//////////////////////////////////////////////////////////////////////////////////////////////////

    void spiralarms(double x, double y, double z, double r, double t) { // Spiral arm parameters for Cox & Gomez (2002), A&A
    	static double n, alpha, h, rs, r0, a, omega;
    	if(init) {
        n = 2;
        alpha = 15.5*M_PI/180.0;
        h = 300. / u->pc;
        rs = 2600.0 / u->pc;
        r0 = 5600.0 / u->pc;
        a = 0.0336 / u->msun * pow(u->pc, 3); // Msun / pc^3
        omega = 20.0 * 0.001023 * u->myr; // km/s/kpc -> Myr^-1
    	}
    	
      double theta = -omega*t;
      double xp = cos(theta)*x - sin(theta)*y;
      double yp = sin(theta)*x + cos(theta)*y;
      double kh = h*n/(r*sin(alpha));
    	double beta = kh*(1.+0.4*kh);
    	double d = (1. + kh + 0.3*kh*kh) / (1.+0.3*kh);

      phi += (-4.*M_PI*h*a*exp((r0-r)/rs) * cos(2.*(atan2(yp,xp) - log(r/r0)/tan(alpha)))) / (pow(cosh(kh*z/beta), beta)*kh*d);
    } 

//////////////////////////////////////////////////////////////////////////////////////////////////

    void nfwcosmo(double r, double t) { // Buist & Helmi
    	static double t0, omega_m, p1, h0, p2, mgm, rs;
    	if(init) {
    		t0 = 1169.16 / u->myr; // age of the universe at t = 0
				omega_m = 0.31; // Omega_matter
				p1 = pow((1.0-omega_m)/omega_m, 1.0/3.0); // (1 - Omega_m / Omega_m )^(1/3)
    		h0 = 68.0 * 1.023 * 1.0e-6 * u->myr; // Hubble constant km/s/Mpc
				p2 = 2.0/(3.0 * h0 * sqrt(1.0-omega_m)); // 2 / 3  * 1/ ( H0 * sqrt(1-Omega_m)  )
        mgm = 1.5e11 / u->msun;
        rs = 16000. / u->pc;
    	}
    	double redshift = p1 * pow(sinh((t+t0)/p2),-2.0/3.0) -1.0; // normal evolution	
			// double redshift = p1 * pow(sinh((12563.2/u->myr - t +t0)/p2),-2.0/3.0) -1.0;  // reverse evolution
			// double redshift = 5.; // static potential
    	
    	double growth = exp(-0.1*redshift); // normal growth
			// double growth = exp(-0.1*pow(5.,0.8)*pow(redshift,0.2)); // late growth

    	phi += -mgm * growth*growth / r * log(1.0 + r/(rs*growth));
    } 

//////////////////////////////////////////////////////////////////////////////////////////////////

};

#endif
