//Bachelor project by Lina Warntoft 2018
/* In this main part of the project, the variables will be asked for and then fed into the called functions declared and defined in starcluster2.h */

#include "starcluster2.h"

using namespace std;

int main(){
	
	//For measuring run time
	clock_t start_time = clock();
	
	//Initial values from input file
	int N;
	double rh, rc, R_hj, kappa;
	
	ifstream input_values;
	input_values.open("InputFiles/test_in.txt");

	if (input_values.fail()){
		cout << "Error: The input file did not open correctly." << endl;
		exit(1);
	}

	string dummy;
	bool found = false;
	while(found == false){
		getline(input_values, dummy);
		if (dummy.find("kappa:") != string::npos){
		found = true;
		}
	}

	input_values >> N >> rh >> rc >> R_hj >> kappa;
	cout << "N_0: " <<  N << ", r_h0: " << rh << " r_c0: " << rc << ", R_hj0: " << R_hj << ", kappa_0: " << kappa << endl;
	input_values.close();
	
	
	//Evaluate m
	double m = 1.0/(double)N;
	
	//Evaluate energy
	double E = -G*kappa*pow((double)N*m, 2.0)/rh;
	
	//Variables evaluated by equations
	double trh;
	double trc;
	double rv;
	double R_ch;
	double R_ch_min;
	double R_vj;
	double xi_e;
	double mu;
	double delta;
	double Kappa;
	double epsilon;
	double lambda;
	double rho_c;
	double kappa0;
	
	//Point where half-mass crossing time t_cr~t_rh and balanced evolution is no longer a valid assumption:
	int N_end = 200;
	
	//To determine escape rate alpha (different depending on initial conditions)
	double alpha = R_hj*pow((double)N*m, 1.0/3.0)/rh;
	cout <<"Constant rate: " << alpha << endl;
	
	//Parameters for while-loop
	double t_step;
	long int count = 0;
	bool phase = false;
	double tot_time = 0;
	double end_time = pow(10,25); //something big
	double R_ch_N = eval_Rch(N);
	double R_ch_r = eval_Rch(rc, rh);
	
	
	//Restults being streamed to a file
	ofstream output_results;
	output_results.open("Results/test_out.txt");
  
	if (output_results.fail()){
		cout << "Error: The output file did not open correctly." << endl;
		exit(1);
	}
	
	output_results << __TIMESTAMP__ << " Result from one simulation by StarCluster2. The numbers displayed below are ordered: 1. Unbalanced (0) or balanced phase (1), 2. Number of stars (N), 3. cluster half-mass radius (r_h), 4. core radius (r_c), 5. total energy of the cluster (E)," 
					<< "6. form factor depending on density profile (kappa), 7. escape rate (xi_e), 8. ratio for virial radius to Jacobi radius (R_vJ), 9. ratio for radius to Jacobi radius (R_hj), 10. the adaptive time step,"
					<< "11. core relaxation time (t_rc) and 12. the time. The initial values were set to: N_0 = " << N << ", r_h0 = " << rh << ", R_hj0 = " << R_hj << ", kappa_0 = " << kappa << "." << endl << endl;
	output_results << phase << " " << N << " " << rh << " " << rc << " " << E << " " << kappa << " " << xi_e << " " << R_vj << " " << R_hj << " " << t_step << " " << trc << " " << tot_time << endl;
	

	/*Evaluating cluster formation using the differential equations and 4th order Runge-Kutta

	- dN/dt = -N*xi_e/t_rh
	- dr_h/dt = r_h*mu/t_rh
	- dE/dt = -E*epsilon/t_rh 
	- dkappa/dt = kappa*lambda/t_rh
	- dr_c/dt = r_c*delta/t_rh

	Ends before N<200 and reasonable time. */
	
	while ((N >= N_end)&&(tot_time < end_time)){
		//UNBALANCED PHASE
		if (!(R_ch_r <= R_ch_N)){ //Condition for UB or B
			R_hj = eval_Rhj(alpha, m, N, rh);
			R_ch = R_ch_r;
			R_ch_min = R_ch_N;
			rho_c = eval_rho_c_pre(rc);
			rv = eval_rv(kappa, rh);
			R_vj = eval_Rvj(alpha, m, N, rv);
			kappa0 = eval_kappa0(rh, rv);
			Kappa = eval_Kappa_pre(kappa, kappa0, R_ch, rh, rv);
			xi_e = eval_xie_pre(N, R_ch, R_ch_min, R_vj);
			
			trh = eval_trh(m, N, rh);
			trc = eval_trc(m, N, rc, rho_c);
			
			delta = eval_delta_pre(trc, trh);
			mu = eval_mu_pre(delta, kappa, Kappa, R_hj, xi_e);
			lambda = eval_lambda_pre(delta, Kappa, mu);
			epsilon = eval_epsilon_pre(lambda, mu, xi_e);
			
			phase = false;
			
			t_step = 1.0/(1.0/(100.0*trc) + 1.0/(0.1*trh));
			tot_time += t_step;
			
			//Maybe add if-statement to only print every 0.1*trh
			N = rungekutta4(t_step, trh, xi_e, N, diff_N);
			rh = rungekutta4(t_step, trh, mu, rh, diff_rh);
			E = rungekutta4(t_step, trh, epsilon, E, diff_E);
			kappa = rungekutta4(t_step, trh, lambda, kappa, diff_kappa);
			rc = rungekutta4(t_step, trh, delta, rc, diff_rc);
			
			output_results << phase << " " << N << " " << rh << " " << rc << " " << E << " " << kappa << " " << xi_e << " " << R_vj << " " << R_hj << " " << t_step << " " << trc << " " << tot_time << endl;
			count += 1;
			R_ch_N = eval_Rch(N);
			R_ch_r = eval_Rch(rc, rh);
		}
		//BALANCED PHASE
		else{
			R_hj = eval_Rhj(alpha, m, N, rh);
			R_ch = eval_Rch(rc, rh);
			rho_c = eval_rho_c_post(m, N, R_ch, rh);
			rv = eval_rv(kappa, rh);
			R_vj = eval_Rvj(alpha, m, N, rv);
			Kappa = eval_Kappa_post(kappa, R_ch);
			xi_e = eval_xie_post(alpha, m, N, rh, R_vj);
			epsilon = eval_epsilon_post();
			mu = eval_mu_post(N, Kappa, xi_e);
			delta = eval_delta_post(mu, N, xi_e);
			lambda = eval_lambda_post(delta, kappa, Kappa, mu, R_ch);
			
			trh = eval_trh(m, N, rh);
			
			phase = true;
			
			t_step = 0.1*trh;
			tot_time += t_step;
			
			N = rungekutta4(t_step, trh, xi_e, N, diff_N);
			rh = rungekutta4(t_step, trh, mu, rh, diff_rh);
			E = rungekutta4(t_step, trh, epsilon, E, diff_E);
			kappa = rungekutta4(t_step, trh, lambda, kappa, diff_kappa);
			rc = rungekutta4(t_step, trh, delta, rc, diff_rc);
			
			output_results << phase << " " << N << " " << rh << " " << rc << " " << E << " " << kappa << " " << xi_e << " " << R_vj << " " << R_hj << " " << t_step << " " << trc << " " << tot_time << endl;
			count += 1;
		}
	}
	
	output_results.close();
	
	cout << "Iterations: " << count << endl << "Run time: " << (double)(clock()-start_time)/CLOCKS_PER_SEC << " sec" << endl;
	
	return 0;
}