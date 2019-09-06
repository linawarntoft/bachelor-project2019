//Bachelor project by Lina Warntoft 2018
/* This is the first version of the code for simulating star cluster formation and evolution. The differential equations, the argument-evaluating functions and the Runge-Kutta 4 function will be declared and defined in this header file.
 
 The differential equations are:
 - dN/dt = -N*xi/t_rh
 - dr/dt = r*mu/t_rh
 - dE/dt = |E|*zeta/t_rh */

#ifndef STAR_CLUSTER_1
#define STAR_CLUSTER_1

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>
#include <time.h>


//------------------------------- Global constants ---------------------------------------

const double KAPPA = 0.25;
const double ZETA = 0.111;
const double GAMMA = 0.11;
const double G = 1;
const double XI1 = 0.0142;
const int N1 = 38252;
const double R1 = 0.145;
const double X = 0.75;
const double Z = 1.61;
//const double MG = 1.5*pow(10.0,5.0);
//const double RG = 7706.9;

//----------------------------------------------------------------------------------------

//------------------------------- Differential equations ---------------------------------

//Differential equation for N
double diff_N(double var_trh, double var_xi, double var_N){
  return -var_N*var_xi/var_trh;
}

//Differential equation for r_h
double diff_r(double var_trh, double var_mu, double var_r){
  return var_r*var_mu/var_trh;
}

//Differential equation for E
double diff_E(double var_trh, double var_zeta, double var_E){
  return fabs(var_E)*var_zeta/var_trh;
}


//-----------------------------------------------------------------------------------------

//------------------------------- Evaluate arguments --------------------------------------

//Function to evaluate mu. Arguments: Xi
double eval_mu(const double var_xi){
  return (ZETA - 2.0*var_xi);
}

//Function to evaluate t_rh. Argumens: N, r, m
double eval_trh(const double var_N, const double var_r, const double var_m){
  double rh = 4.0*var_r*KAPPA;
  return 0.138*(sqrt((double)(var_N*pow(rh, 3.0)/var_m))/log((double)(GAMMA*var_N)));
}

//Function to evaluate escape rate R. Arguments: N, m, r, rate constant
/*double eval_R(const double var_N, const double var_r, const double var_m, const double var_alpha){
	double rj = RG*pow(G*var_N*var_m/(3*MG),1.0/3.0);
	return var_r/rj;
}*/

//Function to evaluate escape rate R. Arguments: N, r, m and alpha
double eval_R(const double var_N, const double var_r, const double var_m, const double var_alpha){
  return var_alpha*var_r/pow((double)(var_N*var_m), 1.0/3.0);
}

//Function to evaluate xi. Argumens: r, N, m
double eval_xi(const double var_N, const double var_r, const double var_m, const double var_alpha){
  double R = eval_R(var_N, var_r, var_m, var_alpha);
  double ln_1 = log((double)(GAMMA*N1));
  double ln_2 = log((double)(GAMMA*var_N));
  double P = pow(R/R1, Z)*pow((double)(var_N*ln_1/(N1*ln_2)), 1.0-X);
  return XI1*(1.0-P) + 3.0/5.0*ZETA*P;
} 


//-----------------------------------------------------------------------------------------

//------------------------------- 4th order Runge-Kutta -----------------------------------


//Runge-Kutta function to evaluate differential equations at one timestep.
double rungekutta4(const double h, const double time, const double arg, const double y0, double (*diff)(double, double, double )){

  double k1 = (*diff)(time, arg, y0);
  double k2 = (*diff)(time, arg, y0+k1*h/2.0);
  double k3 = (*diff)(time, arg, y0+k2*h/2.0);
  double k4 = (*diff)(time, arg, y0+k3*h);
  
  double y1 = y0 + h/6.0*(k1 + k2*2.0 + k3*2.0 + k4);
  
  return y1;
}

//-----------------------------------------------------------------------------------------


#endif  // STAR_CLUSTER_1