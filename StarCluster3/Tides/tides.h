#ifndef TIDES_HEAD
#define TIDES_HEAD

#include <fstream>
#include "mygalaxy.h"
#include "units.h"

void maineigenval(double*);

////////////////////////////////////////////////////////////////////////////////////////////////////////
class TIDES {  // General tides class (for either Mode A or B)
  public:
    double dt; // timestep
    double lambda1;
    UNITS *u;

    virtual ~TIDES() =0;
    virtual void update(double) =0;  // update the tensor (mode A) or update the position of the cluster (mode B)
    virtual void acceleration(double, double, double, double&, double&, double&) const =0; // returns the tidal acceleration at current time, and given position
    virtual void dump(double) =0; // output tidal information in a separate file

  protected:
    std::ofstream dumpfile;
    
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
class TIDESA : public TIDES { // inherit from TIDES
  public:
    TIDESA(UNITS*);
    ~TIDESA();
    void update(double);
    void acceleration(double, double, double, double&, double&, double&) const;
    void dump(double);

  protected:
    unsigned long int n; // number of tensors;
    double *ttime; // timestamps of stored tensors
    double **ttall; // table of all tensors and lambda1's
    unsigned long int inow; // index of the next tensor (saved here for speed up)
    double tt[7]; // current tensor (and lambda1 during interpolation)
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
class TIDESB : public TIDES { // inherit from TIDES
  public:
    TIDESB(UNITS*, double, double, double, double, double, double); //Changed from TIDES to TIDESB 
    ~TIDESB(); //Changed from TIDES to TIDESB
    void update(double);
    void acceleration(double, double, double, double&, double&, double&) const;
    void dump(double);

  protected:
    MYGALAXY* mygal; // mygalaxy model
    double time; // time used in the evaluation of the potential (this is updated in update)
    double dx; // spacestep
    double x, y, z, vx, vy, vz, ax, ay, az; // position, velocity, acceleration of the cluster centre within the galaxy

    // 2nd order stencil used in finite differences to compute the galactic acceleration and tensor (Warning: different order than in Nbody6tt)
    //                      0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
    double xstencil[24] = { 0, 0,-1, 0, 0, 1, 0, 0,-2, 0, 0, 2, 0,-1, 0, 1, 0,-1, 0, 1,-1, 1,-1, 1};
    double ystencil[24] = { 0,-1, 0, 0, 1, 0, 0,-2, 0, 0, 2, 0,-1, 0, 1, 0,-1, 0, 1, 0,-1,-1, 1, 1};
    double zstencil[24] = {-1, 0, 0, 1, 0, 0,-2, 0, 0, 2, 0, 0,-1,-1,-1,-1, 1, 1, 1, 1, 0, 0, 0, 0};
    void galacticacceleration(bool, double, double, double, double, double&, double&, double&, double&) const;
};

#endif