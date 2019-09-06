#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>
#include <time.h>

/// Functions related to the galaxy: tidal accelerations and orbital motion of the cluster

#include "mygalaxy.h"
#include "tides.h"



TIDES::~TIDES() {}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

TIDESA::TIDESA(UNITS* uni) { // constructor: read input file and compute corresponding lambda1's. Set initial tensor and lambda1. Open dump file
  long unsigned int i;
  double unit, offset, unit2;
  double dummy;

  u = uni;
  inow = 0;
  std::ifstream in;
  in.open("tt.dat");
  in >> n >> unit >> offset;
  unit2 = unit*unit;

  ttime = new double[n];
  ttall = new double *[n];
  for(i=0; i<n; ++i) {
    ttall[i] = new double[7];
  }
  
  for(i=0; i<n; ++i) {
    in >> ttime[i] >> ttall[i][0] >> ttall[i][1] >> ttall[i][2] >> dummy >> ttall[i][3] >> ttall[i][4] >> dummy >> dummy >> ttall[i][5];

    // convert units  
    ttime[i] = ttime[i] * unit + offset;
    for(int j=0; j<6; ++j) {
      ttall[i][j] = ttall[i][j] / unit2;
    }
    maineigenval(ttall[i]); // store max eigenvalue in ttall[i][6];
  }
  in.close();
  update(0.0);

  dumpfile.open("tides_tl1.out", std::ios::out | std::ios::trunc);
}

////////////////////////////////////////////////////

TIDESA::~TIDESA() { // destructor: close dump file and free memory
  dumpfile.close();
  for(int i=0; i<n; ++i) {
    delete[] ttall[i];
  }
  delete[] ttall;
  delete[] ttime;
}

////////////////////////////////////////////////////

void TIDESA::update(double t) { // update tensor and lambda1 by quadratic interpolation on the input table
  double dt21, dt32, dt31, dtdtdt;
  double tt1, tt2, tt3, tta, ttb;
  int j;

  if(t == ttime[0])  { // at first snapshot
    for(j=0; j<7 ; ++j) {
      tt[j] = ttall[0][j];
    }
    lambda1 = tt[6];
    return;
  }

  if((t < ttime[0]) || (t >= ttime[n-1]))  { // before first snapshot or after last snapshot: no tides
    for(j=0; j<6 ; ++j) {
      tt[j] = 0.;
    }
    lambda1 = 0.;
    return;
  }

  while(t > ttime[inow]) { // update the index of the next tensor
    ++inow;
  }

  dt21 = ttime[inow] - ttime[inow-1];  
  if(inow == n-1) {  // between the second to last and the last snapshots
    dt32 = dt21;
    dt31 = 2.*dt21;
  } else { // normal
    dt32 = ttime[inow+1] - ttime[inow];
    dt31 = ttime[inow+1] - ttime[inow-1];
  }

  dtdtdt = dt21*dt32*dt31;
  for(j=0; j<7 ; ++j) { // interpolate tensor and lambd1
    tt1 = ttall[inow-1][j];
    tt2 = ttall[inow][j];
    tta = (tt2-tt1)/dt21;
    if(inow == n-1) {
      tt3 = 0.0;
    } else {
      tt3 = ttall[inow+1][j]; 
    }
    ttb = ((tt3-tt1)/dt21-(tt2-tt1)*dt31)/dtdtdt;
    tt[j] = tt1 + tta * (t-ttime[inow-1]) + ttb * (t-ttime[inow-1])*(t-ttime[inow]);
  }
  lambda1 = tt[6];
}

////////////////////////////////////////////////////

void TIDESA::acceleration(double px, double py, double pz, double &pax, double &pay, double &paz) const { // tidal acceleration from *current* tensor at position x, y, z
  pax = tt[0]*px + tt[1]*py + tt[2]*pz; 
  pay = tt[1]*px + tt[3]*py + tt[4]*pz;
  paz = tt[2]*px + tt[4]*py + tt[5]*pz;
}

////////////////////////////////////////////////////

void TIDESA::dump(double t) {
  dumpfile << t*u->myr << " " << lambda1*u->myrm2 << "\n";
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

TIDESB::TIDESB(UNITS* uni, double initx, double inity, double initz, double initvx, double initvy, double initvz) { // constructor: set initial position, velovity, acceleration and lambda1 on cluster centre
  u = uni;
  mygal = new MYGALAXY(u); // initialize mygalaxy

  x = initx;
  y = inity;
  z = initz;
  vx = initvx;
  vy = initvy;
  vz = initvz;
  dx = 4e-5*sqrt(x*x+y*y+z*z);
  time = 0.;

  galacticacceleration(true, time, x, y, z, ax, ay, az, lambda1); // get acceleration and lambda1 (but do not change position and velocity)

  dumpfile.open("tides_txyzl1.out", std::ios::out | std::ios::trunc);
}

////////////////////////////////////////////////////

TIDESB::~TIDESB() { // destructor: close file and free memory
  dumpfile.close();
  delete mygal;
}

////////////////////////////////////////////////////

void TIDESB::update(double t) { // update position, acceleration and lambda1 of cluster centre
  double l1;

  // Leap-frog
  x = x + dt * vx;
  y = y + dt * vy;
  z = z + dt * vz;
  galacticacceleration(true, t, x, y, z, ax, ay, az, l1);
  vx = vx + dt * ax;
  vy = vy + dt * ay;
  vz = vz + dt * az;

  lambda1 = l1;
  dx = 4e-5 * sqrt(x*x+y*y+z*z); // update spacestep
  time = t; // store the time at which the cluster centre has been updated
}

////////////////////////////////////////////////////

void TIDESB::acceleration(double px, double py, double pz, double &pax, double &pay, double &paz) const { // get tidal acceleration as differential acceleration (= at position - at cluster centre), using *current* acceleration on cluster centre
  double tmp;
  galacticacceleration(false, time, x+px, y+py, y+pz, pax, pay, paz, tmp); // no need to update lambda1
  pax = pax - ax;
  pay = pay - ay;
  paz = paz - az;
}

///////////////////////////////////////////////////

void TIDESB::galacticacceleration(bool getlambda1, double t, double px, double py, double pz, double &pax, double &pay, double &paz, double &l1) const {
  double tt[7];
  double dx2sq, dx12;
  double phi[25], phic;
  
  // galactic potential at the points needed for acceleration only
  for(int i=0 ; i< 12 ; ++i) {
    phi[i] = mygal->potential(px+xstencil[i]*dx, py+ystencil[i]*dx, pz+zstencil[i]*dx, t);
  }

  // galactic acceleration
  dx12 = 12.*dx;
  pax = (phi[11]-8.*phi[ 5]+8*phi[ 2]-phi[ 8])/dx12;
  pay = (phi[10]-8.*phi[ 4]+8*phi[ 1]-phi[ 7])/dx12;
  paz = (phi[ 9]-8.*phi[ 3]+8*phi[ 0]-phi[ 6])/dx12;
  

  if(getlambda1) {
    for(int i=12 ; i< 24 ; ++i) { // galactic potential at the additional points needed for tidal tensor
      phi[i] = mygal->potential(px+xstencil[i]*dx, py+ystencil[i]*dx, pz+zstencil[i]*dx, t);
    }
    phic = mygal->potential(px, py, pz, t);
  
    dx2sq = 4.*dx*dx;
    tt[0] = (2*phic-phi[11]-phi[ 8]) / dx2sq; // Txx
    tt[1] = (phi[21]-phi[23]-phi[20]+phi[22]) / dx2sq; // Txy
    tt[2] = (phi[17]-phi[13]-phi[19]+phi[15]) / dx2sq; // Txz
    tt[3] = (2*phic-phi[10]-phi[ 7]) / dx2sq; // Tyy
    tt[4] = (phi[16]-phi[18]-phi[12]+phi[14]) / dx2sq; // Tyz
    tt[5] = (2*phic-phi[ 9]-phi[ 6]) / dx2sq; // Tzz
    maineigenval(tt); // store max eigenvalue in tt[6];
    l1 = tt[6];
  }
}

////////////////////////////////////////////////////

void TIDESB::dump(double t) {
  dumpfile << t*u->myr << " " << x*u->pc << " " << y*u->pc << " " << z*u->pc << " " << lambda1*u->myrm2 << "\n";
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void maineigenval(double *tt) { // Diagonalisation algorithm from Kopp (2006, paragraph 3)
  double c2,c1,c0,p,q,r,s,d,phi,cphi,sphi;
  c2 = -tt[0]-tt[3]-tt[5];
  c1 = tt[0]*tt[3] + tt[0]*tt[5] + tt[3]*tt[5] - tt[1]*tt[1] - tt[2]*tt[2] - tt[4]*tt[4];
  c0 = tt[0]*tt[4]*tt[4] + tt[3]*tt[2]*tt[2] + tt[5]*tt[1]*tt[1] - tt[0]*tt[3]*tt[5] - 2.*tt[2]*tt[1]*tt[4]; 
  p = c2*c2 - 3.*c1;
  q = -13.5*c0 - c2*c2 + 4.5*c2*c1; 
  if (p<0) {
    r=0;
  } else {
    r = sqrt(p)/3.0;
  }
  s = -c2/3.0;
  d = p*p*p - q*q;
  if (q == 0) {
    phi = M_PI/6.;
  } else {
    if (d<0) {
      phi = 0;
    } else {
      phi = atan(sqrt(d)/q) / 3.0;
    }
  }
  cphi = cos(phi);
  sphi = sin(phi);
  tt[6] = r*fmax(fmax(2*cphi, -cphi - sqrt(3)*sphi),-cphi + sqrt(3)*sphi)+s;  // maximum eigenvalue
  // get others to compute density for DF.
}


