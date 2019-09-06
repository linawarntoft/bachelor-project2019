#ifndef UNITS_H
#define UNITS_H

#include <cmath>

class UNITS { // convert to code = val / units.x        convert to physical = val * units.x
  public:
    double pc, msun, kms, myr, myrm2;  // units (G=1, 1 pc, 1 Msun). Code units * unit.X = pc, km/s, Msun, Myr, Myr^-2

    UNITS(double scale_l=1.0, double scale_m=1.0) {
			pc = scale_l;
			msun = scale_m;
			kms = sqrt(4.3010795e-3 * msun / pc);
			myr = pc/(kms*1.023);
			myrm2 = 1./(myr*myr);
		}

};

#endif