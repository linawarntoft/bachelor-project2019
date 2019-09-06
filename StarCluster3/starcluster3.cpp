//Bachelor project by Lina Warntoft 2018
/* In this main part of the project, the variables will be asked for and then fed into the called functions declared and defined in starcluster3.h */

#include "starcluster3.h"
#include "potentialEscaper.h"
#include "Tides/mygalaxy.h"
#include "Tides/tides.cpp"
#include "Tides/units.h"

using namespace std;

int main(){
	
	//For measuring run time
	//clock_t start_time = clock();
	
	//=================== Initial values from input file ===================
	
	int N;
	double rh, rc, kappa;
	double follow_esc; //true if we want to follow the escapers
	double N_strip; //the posibility of cluster swallowing escapers again
	double x_init = 0, y_init = 0, z_init = 0; //Only for ModeB
	double vx_init = 0, vy_init = 0, vz_init = 0; //Only for ModeB
	
	ifstream input_values;
	input_values.open("InputFiles/test.txt");

	if (input_values.fail()){
		cout << "Error: The input file did not open correctly." << endl;
		exit(1);
	}

	string dummy;
	bool found = false;
	while(found == false){
		getline(input_values, dummy);
		if (dummy.find("v_z:") != string::npos){
		found = true;
		}
	}

	input_values >> N >> rh >> rc >> kappa >> follow_esc >> N_strip;
	cout << "N_0: " <<  N << ", r_h0: " << rh << ", r_c0: " << rc << ", kappa_0: " << kappa << endl;
	if (!(input_values.eof())){
		input_values >> x_init >> y_init >> z_init >> vx_init >> vy_init >> vz_init;
	}
	input_values.close();
	
	
	//======================================================================
	//============================ Convert units ===========================
	UNITS u;
	x_init = x_init/u.pc;
	y_init = y_init/u.pc;
	z_init = z_init/u.pc;
	
	vx_init = vx_init/u.kms;
	vy_init = vy_init/u.kms;
	vz_init = vz_init/u.kms;
	
	rh = rh/u.pc;
	rc = rc/u.pc;
	
	//======================================================================
	//============================= Choose Mode ============================
	
	TIDES* tides;
	
	if (sqrt(x_init*x_init + y_init*y_init + z_init+z_init) != 0.0){
		tides = new TIDESB(&u, x_init, y_init, z_init, vx_init, vy_init, vz_init);
		cout << "Mode B was used." << endl;
	}
	else{
		tides = new TIDESA(&u);
		cout << "Mode A was used." << endl;
	}
	
	//======================================================================
	//============================== Variables =============================
	//Initial N
	int N_0 = N;
	
	//Evaluate m
	double m = 1.0/(double)N;
	
	//Evaluate energy
	double E = -G*kappa*pow((double)N*m, 2.0)/rh;
	
	//Variables evaluated by equations
	double trh, trc;
	double rv;
	double R_ch, R_ch_min, R_vj;
	double xi_e, mu, delta, epsilon, lambda;
	double Kappa, kappa0, rho_c;
	double R_hj, beta;
	
	//Point where half-mass crossing time t_cr~t_rh and balanced evolution is no longer a valid assumption:
	int N_end = 200;
	
	//Parameters for while-loop
	double t_step;
	long int count = 0;
	bool phase = false;
	double tot_time = 0;
	double end_time = pow(10,25); //something big
	double R_ch_N = eval_Rch(N);
	double R_ch_r = eval_Rch(rc, rh);
	
	double smaller_dt = 0.001;
	
	//Potential escapers
	int PE_count = 0; //count PE's
	double dE, E_pot, E_k;
	int N_pre;
	int escaped = 0; //counts escaped stars
	int PE_new = 0;
	
	PotEsc pot_escaper; //initialize PE-class object
	
	//======================================================================
	//============================= Output file ============================
	
	//Restults being streamed to a file
	ofstream output_results;
	output_results.open("Results/test.txt");
  
	if (output_results.fail()){
		cout << "Error: The output file did not open correctly." << endl;
		exit(1);
	}
	
	output_results << __TIMESTAMP__ << " Result from one simulation by StarCluster2. The numbers displayed below are ordered: 1. Unbalanced (0) or balanced phase (1), 2. Number of stars (N), 3. cluster half-mass radius (r_h), 4. core radius (r_c), 5. total energy of the cluster (E)," 
					<< "6. form factor depending on density profile (kappa), 7. escape rate (xi_e), 8. ratio for virial radius to Jacobi radius (R_vJ), 9. ratio for radius to Jacobi radius (R_hj), 10. number of potential escapers, 11 escaped stars, 12. the adaptive time step,"
					<< "11. core relaxation time (t_rc) and 12. the time. The initial values were set to: N_0 = " << N << ", r_h0 = " << rh << ", R_hj0 = " << R_hj << ", kappa_0 = " << kappa << "." << endl << endl;
	output_results << phase << " " << N << " " << rh << " " << rc << " " << E << " " << kappa << " " << xi_e << " " << R_vj << " " << R_hj << " " << PE_count << " " << escaped << " " << t_step << " " << trc << " " << tot_time*u.myr << endl;
	
	
	/*ofstream outfile; //Just to reset coords output file
	outfile.open("Results/coords.txt");
	if (outfile.fail()){
		cout << "Error: The output file did not open correctly in potentialEscapers.h." << endl;
		exit(1);
	}
	outfile << endl;
	outfile.close();*/
	//======================================================================
	//============================= While-loop =============================

	/*Evaluating cluster formation using the differential equations and 4th order Runge-Kutta

	- dN/dt = -N*xi_e/t_rh
	- dr_h/dt = r_h*mu/t_rh
	- dE/dt = -E*epsilon/t_rh 
	- dkappa/dt = kappa*lambda/t_rh
	- dr_c/dt = r_c*delta/t_rh

	Ends before N<200 and reasonable time. */
	
	while (((N_0-escaped) >= N_end)&&(tot_time < end_time)){ //LINA
		phase = true; //True for balanced phase
		R_ch = eval_Rch(rc, rh);
		
		//----------------- Check condition for UB and B ---------------
		if (!(R_ch_r <= R_ch_N)){ //Condition for UB or B
			phase = false; //False for unbalanced phase
			R_ch = R_ch_r;
			R_ch_min = R_ch_N;
		}
		//------------------ Update coupled variables ------------------
		tides->update(tot_time);
		tides->dump(tot_time);
		R_hj = pow(G*(double)N*m/tides->lambda1, 1.0/3.0); //LINA START
		beta = R_hj*pow((double)N*m, 1.0/3.0)/rh; //??? LINA
		rho_c = eval_rho_c(m, N, R_ch, rc, rh, phase); //LINA
		rv = eval_rv(kappa, rh);
		R_vj = eval_Rvj(beta, m, N, rv);
		kappa0 = eval_kappa0(rh, rv, phase);
		Kappa = eval_Kappa(kappa, kappa0, R_ch, rh, rv, phase);
		xi_e = eval_xie(beta, m, N, R_ch, R_ch_min, rh, R_vj, phase);
		
		trh = eval_trh(m, N, rh);
		trc = eval_trc(m, N, rc, rho_c, phase);
		
		if (phase){
			mu = eval_mu(delta, kappa, Kappa, N, R_hj, xi_e, phase);
			delta = eval_delta(mu, N, xi_e, trc, trh, phase);
		}
		else{
			delta = eval_delta(mu, N, xi_e, trc, trh, phase);
			mu = eval_mu(delta, kappa, Kappa, N, R_hj, xi_e, phase);
		}
		lambda = eval_lambda(delta, kappa, Kappa, mu, R_ch, phase);
		epsilon = eval_epsilon(lambda, mu, xi_e, phase);
		
		
		dE = E; //For updating PE's
		N_pre = N; //For updating PE's
		
		//--------------------- Adaptive step size ---------------------
		if (phase){
			t_step = 0.1*trh;
		}
		else {
			t_step = 1.0/(1.0/(100.0*trc) + 1.0/(0.1*trh));
		}
		t_step = t_step*smaller_dt;
		tot_time += t_step;
		tides->dt = t_step;
		
		//----------------------- Runge-Kutta 4 ------------------------
		
		N = rungekutta4(t_step, trh, xi_e, N, diff_N);
		rh = rungekutta4(t_step, trh, mu, rh, diff_rh);
		E = rungekutta4(t_step, trh, epsilon, E, diff_E);
		kappa = rungekutta4(t_step, trh, lambda, kappa, diff_kappa);
		rc = rungekutta4(t_step, trh, delta, rc, diff_rc);
		
		//---------------------- Generating PE's -----------------------
		if (PE_count != 0){
			pot_escaper.leap_frog(tot_time, t_step, m, N, rh, tides, follow_esc, N_strip);
			//pot_escaper.fstream_coords(u);
		}
		dE = E - dE;
		E_pot = -G*(double)N*m/(2*rh); //POINTMASS
		PE_new = N_pre - N;
		dE = dE/(double)PE_new;
		E_k = dE - E_pot;
		for (int i = 0; i < PE_new; i++){
			pot_escaper.append(tot_time, t_step, m, N, rh, E_k, tides, N_strip);
		}
		PE_count += PE_new;
		escaped = pot_escaper.count_loss();
		
		//----------------------- Output results -----------------------
		output_results << phase << " " << N << " " << rh << " " << rc << " " << E << " " << kappa << " " << xi_e << " " << R_vj << " " << R_hj << " " <<  PE_count << " " << escaped << " " << t_step << " " << trc << " " << tot_time*u.myr << endl;
		count++;
		if (!(phase)){
			R_ch_N = eval_Rch(N);
			R_ch_r = eval_Rch(rc, rh);
		}
		//==================================================================
	}
	
	output_results.close();
	delete tides;
	
	//======================================================================
	//========================== Output to terminal ========================
	cout << "PE's: " << PE_count << ", escaped: " << escaped << endl;
	cout << "Iterations: " << count << endl << "Run time: " << (double)(clock()-start_time)/CLOCKS_PER_SEC << " sec" << endl;
	//cout << "CPS: " << CLOCKS_PER_SEC << ", start time: " << start_time << ", clock: " << clock() << endl;
	//======================================================================
	
	return 0;
}
