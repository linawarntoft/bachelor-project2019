//Bachelor project by Lina Warntoft 2018
/* This is the third version of the code for simulating star cluster formation and evolution.
 The differential equations, the argument-evaluating functions and the Runge-Kutta 4 function
 will be declared and defined in this header file.
 
 The differential equations are:
 - dN/dt = -N*xi_e/t_rh
 - dr_h/dt = r_h*mu/t_rh
 - dE/dt = -E*epsilon/t_rh 
 - dkappa/dt = kappa*lambda/t_rh
 - dr_c/dt = r_c*delta/t_rh*/
 
 
#ifndef STAR_CLUSTER_3
#define STAR_CLUSTER_3

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>
#include <time.h>

//------------------------------- Global constants ---------------------------------------

const double KAPPA0_B = 0.2;
const double KAPPA1_UB = 0.295;
const double KAPPA1_B = 0.265;
const double ZETA = 0.1;
const double GAMMA = 0.11;
const double G = 1.0;
const double XI_E1 = 0.0142;
const int N1 = 15000;
const int N2 = 12;
const int N3 = 15000;
const double DELTA1 = -0.09;
const double DELTA2 = -0.002;
const double RHO_C0 = 0.055;
const double ALPHA = 2.2;
const double R_CH0_UB = 0.1;
const double R_CH0_B = 0.22;
const double R_VJ1 = 0.145;
const double f = 0.3;
const double X = 0.75;
const double Z = 1.61;

//----------------------------------------------------------------------------------------

//------------------------------- Differential equations ---------------------------------

//Differential equation for N
double diff_N(double var_trh, double var_xi, double var_N){
  return -var_N*var_xi/var_trh;
}

//Differential equation for r_h
double diff_rh(double var_trh, double var_mu, double var_rh){
  return var_rh*var_mu/var_trh;
}

//Differential equation for E
double diff_E(double var_trh, double var_eps, double var_E){
  return -var_E*var_eps/var_trh;
}

//Differential equation for kappa
double diff_kappa(double var_trh, double var_lambda, double var_kappa){
	return var_kappa*var_lambda/var_trh;
}

//Differential equation for r_c
double diff_rc(double var_trh, double var_delta, double var_rc){
	return var_rc*var_delta/var_trh;
}

//-----------------------------------------------------------------------------------------

//------------------------------- Evaluate arguments --------------------------------------

//Function to evaluate t_rh
double eval_trh(const double var_m, const double var_N, const double var_rh){
  return 0.138*(sqrt((double)var_N*pow(var_rh, 3.0)/(G*var_m))/log(GAMMA*(double)var_N));
}

//Function to evaluate t_rc
double eval_trc(const double var_m, const double var_N, const double var_rc,const double var_rhoc, const bool balanced){
	if (balanced){
		return 0;
	}
	else {
		double sigma_c2 = 8.0*M_PI*G*RHO_C0*pow(var_rc, (double)(2.0-ALPHA))/3.0;
		return pow(sigma_c2, 3.0/2.0)/(15.4*pow(G, 2.0)*var_m*var_rhoc*log(GAMMA*(double)var_N));
	}
}

//Function to evaluate R_hJ
double eval_Rhj(const double var_beta, const double var_m, const double var_N, const double var_rh){
  return var_beta*var_rh*pow((double)var_N*var_m, -1.0/3.0);
}

//Function to evaluate R_ch
double eval_Rch(const double var_N){
	return pow(((double)N2/(double)var_N + (double)N2/(double)N3), 2.0/3.0);
}

//Alternative function to evaluate R_ch
double eval_Rch(const double var_rc, const double var_rh){
	return var_rc/var_rh;
}

//Function to evaluate R_vj
double eval_Rvj(const double var_beta, const double var_m, const double var_N, const double var_rv){
	return var_beta*var_rv*pow((double)var_N*var_m, -1.0/3.0);
}

//Function to evaluate r_v
double eval_rv(const double var_kappa, const double var_rh){
	return var_rh/(4.0*var_kappa);
}

//Funnction to evaluate kappa(R_ch)
double eval_kappa_Rch(const double var_Rch){
	return (KAPPA1_B + (KAPPA0_B - KAPPA1_B)*erf(var_Rch/R_CH0_B));
}

//Function to evaluate Kappa
double eval_Kappa(const double var_kappa, const double var_kappa0, const double var_Rch, const double var_rh, const double var_rv, const bool balanced){
	if (balanced){
		return var_Rch*2.0*(KAPPA0_B - KAPPA1_B)*exp(-pow(var_Rch/R_CH0_B, 2.0))/(var_kappa*sqrt(M_PI)*R_CH0_B);
	}
	else{
		return var_Rch*2.0*(var_kappa0 - KAPPA1_UB)*exp(-pow(var_Rch/R_CH0_UB, 2.0))/(var_kappa*sqrt(M_PI)*R_CH0_UB);
	}
}

//Function to evaluate mu
double eval_mu(const double var_delta, const double var_kappa, const double var_Kappa, const double var_N, const double var_Rhj, const double var_xi_e, const bool balanced){
	if (balanced){
		return (ZETA + (2.0/3.0*var_Kappa/(1.0 + (double)var_N/(double)N3) - 2.0)*var_xi_e);
	}
	else{
		return ((var_Rhj/var_kappa - 2.0)*var_xi_e + var_Kappa*var_delta)/(1.0 + var_Kappa);
	}	
}

//Fuction to evaluate epsilon
double eval_epsilon(const double var_lambda, const double var_mu, double const var_xie, const bool balanced){
	if (balanced){
		return ZETA;
	}
	else{
		return (-var_lambda + var_mu + 2*var_xie);
	}
}

//Function to evaluate delta
double eval_delta(const double var_mu, const double var_N, const double var_xi_e, const double var_trc, const double var_trh, const bool balanced){
	if (balanced){
		return (2.0/3.0*var_xi_e/(1.0 + (double)var_N/(double)N3) + var_mu);
	}
	else{
		return (DELTA1 + DELTA2*var_trh/var_trc);
	}
}

//Function to evaluate xi_e
double eval_xie(const double var_beta, const double var_m, const double var_N, const double var_Rch, const double var_Rch_min, const double var_rh, const double var_Rvj, const bool balanced){
	double ln_1 = log(GAMMA*(double)N1);
	double ln_2 = log(GAMMA*(double)var_N);
	double P = pow(var_Rvj/R_VJ1, Z)*pow(((double)var_N/(double)N1*ln_1/ln_2), (1.0-X));
	
	if (balanced){
		return XI_E1*(1.0-P) + 3.0/5.0*ZETA*P;
	}
	else{
		double F = var_Rch_min/var_Rch;
		return F*XI_E1*(1.0 - P) + (f + (1.0 - f)*F)*3.0/5.0*ZETA*P;
	}
}

//Function to evaluate lambda
double eval_lambda(const double var_delta,  const double var_kappa, const double var_Kappa, const double var_mu,  const double var_Rch, const bool balanced){
	if (balanced){
		double kappa_Rch = eval_kappa_Rch(var_Rch);
		return (var_Kappa*(var_delta - var_mu) + (kappa_Rch - var_kappa)/kappa_Rch);
	}
	else{
		return var_Kappa*(var_delta - var_mu);
	}
}

//Function to evaluate rho_c
double eval_rho_c(const double var_m, const double var_N, const double var_Rch, const double var_rc, const double var_rh, const bool balanced){
	if (balanced){
		double rho_h = 3.0*(double)var_N*var_m/(8.0*M_PI*pow(var_rh, 3.0));
		return rho_h/pow(var_Rch, 2.0);
	}
	else{
		return RHO_C0*pow(var_rc, (double)(-ALPHA));
	}
}

//Function to evaluate kappa_0
double eval_kappa0(const double var_rh, double const var_rv, const bool balanced){
	if (balanced){
		return 0;
	}
	else {
		return var_rh/(4.0*var_rv);
	}
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


#endif  // STAR_CLUSTER_3
