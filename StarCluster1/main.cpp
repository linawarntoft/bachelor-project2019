//Bachelor project by Lina Warntoft 2018
/* In this main part of the project, the variables will be asked for and then fed into the called functions declared and defined in starcluster1.h */

#include "starcluster1.h"

using namespace std;

int main(){
  
  //For measuring run time
  clock_t start_time = clock();
  
  //Variables where initial value is fed from input file
  int N;
  double r;
  double R0;
  
  
  ifstream input_values;
  input_values.open("InputFiles/input1.txt");

  if (input_values.fail()){
    cout << "The input file did not open correctly." << endl;
    exit(1);
  }

  string dummy;
  bool found = false;
  while(found == false){
    getline(input_values, dummy);
    if (dummy.find("R_0:") != string::npos){
      found = true;
    }
  }
  
  input_values >> N >> r >> R0;
  cout <<  N << " " << r << " " << R0 << endl;
  input_values.close();
  
  //Evaluate m
  double m = 1.0/N;

  
  //Evaluate energy
  double E = -G*pow(N*m,2.0)/(4.0*r);
  
  //Variables evaluated by equations
  double trh;
  double mu;
  double xi;
  double R;
  
  //Point where half-mass crossing time t_cr~t_rh and balanced evolution is no longer a valid assumption:
  int N_end = 200;
  
  //To determine escape rate alpha (different depending on initial conditions)
  double alpha = R0*pow(N*m, 1.0/3.0)/r;
  cout <<"Constant rate: " << alpha << endl;
  R = R0;
  
  //Parameters for while-loop
  double t_step;
  long int count = 0;
  double tot_time = 0;
  double end_time = pow(10,25); //something big
  
  //Restults being streamed to a file
  ofstream output_results;
  output_results.open("Results/test.txt");
  
  if (output_results.fail()){
    cout << "The output file did not open correctly." << endl;
    exit(1);
  }
  
  output_results << "Result from one simulation by StarCluster1. The numbers displayed below are ordered: Number of stars (N), cluster radius (r), total energy of the cluster (E), ratio for radius to Jacobi radius (R), the adaptive time step and the time. The initial values were set to: N = " << N << ", r = " << r << ", R_0 = " << R << "." << endl << endl;
  
  //Fractional change in N and r from start to after the core collaps
  tot_time = eval_trh(N, r, m)*13.513; //EMACSS 20 not 9.2, turns out about 13.51

  output_results << N << " " << r << " " << E << " " << R << " " << t_step << " " << tot_time << endl;
  
  
  /*Evaluating cluster formation using the differential equations and 4th order Runge-Kutta

  - dN/dt = -N*xi/t_rh
  - dr/dt = r*mu/t_rh
  - dE/dt = E*zeta/t_rh 
 
  Ends before N<200 and reasonable time. */
  
  while ((N >= N_end)&&(tot_time < end_time)){
    trh = eval_trh(N, r, m);
    xi = eval_xi(N, r, m, alpha);
    R = eval_R(N, r, m, alpha); //EMACSS erased rateconst
    mu = eval_mu(xi);
    
    t_step = 0.1*eval_trh(N, r, m); //Adaptive step size
    tot_time += t_step;
    
    N = rungekutta4(t_step, trh, xi, N, diff_N);
    r = rungekutta4(t_step, trh, mu, r, diff_r);
    E = rungekutta4(t_step, trh, ZETA, E, diff_E);
    
    output_results << N << " " << r << " " << E << " " << R << " " << t_step << " " << tot_time << endl;
    count += 1;
  }
  
  output_results.close();
  
  cout << "Iterations: " << count << endl << "Run time: " << (float)(clock()-start_time)/CLOCKS_PER_SEC << " sec" << endl;
  
  return 0;
}