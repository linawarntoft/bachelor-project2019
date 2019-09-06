//Bachelor project by Lina Warntoft 2018

#ifndef POTENTIAL_ESCAPER_H
#define POTENTIAL_ESCAPER_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <time.h>
#include "nodearr.h"
#include "starcluster3.h"
#include "Tides/tides.h"

using namespace std;


//|===================================================================================================================|
//|-------------------- POTENTIAL ESCAPER CLASS FOR GENERATING STARS ABOUT TO LEAVE A STAR CLUSTER -------------------|
//|===================================================================================================================|
class PotEsc{
public:
	PotEsc();
	PotEsc(const double time, const double dtime, const double m, const double N, const double rh, const double E_k, TIDES* tides, const double N_strip);
	~PotEsc();

	//If coordinates or velocities want to be printed to file
	void fstream_coords(UNITS units);
	void fstream_veloc(UNITS units);
	
	//Updates motion of PE's
	void leap_frog(const double time, const double dtime, const double m, const double N, const double rh, TIDES* tides, const double trace, const double N_strip);
	//Adds another PE
	void append(const double time, const double dtime, const double m, const double N, const double rh, const double E_k, TIDES* tides, const double N_strip);
	//Returns amount of PE's that have escaped
	int count_loss();
	
private:
	//Variables saved as linked lists
	double x_, y_, z_;
	double vx_, vy_, vz_;
	double ax_, ay_, az_;
	double dt_ = 0, t_ = 0;
	int index_ = 0;
	
	Node<double>* var_; 
	/*The indices of the variables in the list are:
	0:x, 1:y, 2:z, 3:vx, 4:vy, 5:vz, 6:ax, 7:ay, 8:az,
	9: escaped_stars, 10:index*/
	
	//Checks if PE has escaped
	void potential_escaper(const double N, const double m, const double rh, const double N_strip, TIDES* tides);
	//Generating random position and velocity
	void rand_position(const double E_k, const double m, const double rh);
};
//---------------------------------------------------------------------------------------------------------------------
//------------------------------------------ FUNCTION DEFINITIONS FOR POTESC ------------------------------------------

PotEsc::PotEsc(){
	var_ = nullptr;
}


PotEsc::PotEsc(const double time, const double dtime, const double m, const double N, const double rh, const double E_k, TIDES* tides, const double N_strip){
	index_ = 1;
	dt_ = dtime;
	t_ = time;
	
	//Position and velocity
	rand_position(E_k, m, rh);
	
	//Acceleration
	tides->acceleration(x_, y_, z_, ax_, ay_, az_);
	double eps2 = 0.5874*rh*rh;
	double d = sqrt(x_*x_ + y_*y_ + z_*z_ + eps2);
	double ac = -1.0*m/(d*d*d);
	ax_ = ac*x_ + ax_;
	ay_ = ac*y_ + ay_;
	az_ = ac*z_ + az_;
	
	//Initialize data in node
	double init_list[VARIABLES] = {x_, y_, z_, vx_, vy_, vz_, ax_, ay_, az_, 0.0, index_};
	var_ = new Node<double>(init_list, var_);
	potential_escaper(N, m, rh, N_strip, tides);
}


PotEsc::~PotEsc(){
	list_clear(var_);
}


void PotEsc::fstream_coords(UNITS units){
	ofstream outfile;
	outfile.open("Results/coords.txt", ios::app);
	
	if (outfile.fail()){
		cout << "Error: The output file did not open correctly in potentialEscapers.h." << endl;
		exit(1);
	}
	Node<double>* head = var_;
	while (var_ != nullptr){
		if (var_->get_data(10) == 1.0){
			outfile << var_->get_data(10) << " " << var_->get_data(9) << " " << var_->get_data(0) << " " << var_->get_data(1) << " " << var_->get_data(2) << " " << dt_ << " " << t_*units.myr << endl;
		}
		var_ = var_->get_link();
	}
	var_ = head;
	outfile.close();
}


void PotEsc::fstream_veloc(UNITS units){
	ofstream outfile;
	outfile.open("Results/velocities.txt", ios::app);
	
	if (outfile.fail()){
		cout << "Error: The output file did not open correctly in potentialEscapers.h." << endl;
		exit(1);
	}
	Node<double>* head = var_;
	while (var_ != nullptr){
		if (var_->get_data(10) == 1.0){
			outfile << var_->get_data(10) << " " << var_->get_data(9) << " " << var_->get_data(3) << " " << var_->get_data(4) << " " << var_->get_data(5) << " " << dt_ << " " << t_*units.myr << endl;
		}
		var_ = var_->get_link();
	}
	var_ = head;
	outfile.close();
}


void PotEsc::rand_position(const double E_k, const double m, const double rh){
	
	srand(time(NULL));
	double theta = rand()%361;
	srand(time(NULL) + 10);
	double phi = rand()%181;
	
	x_ = rh*sin(theta*M_PI/180.0)*cos(phi*M_PI/180.0);
	y_ = rh*sin(theta*M_PI/180.0)*sin(phi*M_PI/180.0);
	z_ = rh*cos(theta*M_PI/180.0);
	
	srand(time(NULL) + 20);
	theta = rand()%361;
	srand(time(NULL) + 30);
	phi = rand()%181;
	
	double v_tot = sqrt(2.0*abs(E_k)/m);
	vx_ = v_tot*sin(theta*M_PI/180.0)*cos(phi*M_PI/180.0);
	vy_ = v_tot*sin(theta*M_PI/180.0)*sin(phi*M_PI/180.0);
	vz_ = v_tot*cos(theta*M_PI/180.0);
	
}


void PotEsc::potential_escaper(const double N, const double m, const double rh, const double N_strip, TIDES* tides){
	//Returns true if the PE has escaped
	Node<double>* head = var_;
	double dist;
	double rJ = ((G*(double)N*m)/tides->lambda1, 1.0/3.0);
	
	while (var_ != nullptr){
		dist = sqrt(pow(var_->get_data(0), 2.0) + pow(var_->get_data(1), 2.0) + pow(var_->get_data(2), 2.0));
		if (dist > N_strip*rJ){
			var_->set_data(1.0, 9);
		}
		var_ = var_->get_link();
	}
	var_ = head;
}


void PotEsc::leap_frog(const double time, const double dtime, const double m, const double N, const double rh, TIDES* tides, const double trace, const double N_strip){
	Node<double>* head = var_;
	dt_ = dtime;
	t_ = time;
	double eps2;
	double d;
	double ac;
	double compare;
	
	if (trace == 0.0){
		compare = 1.0;
	}
	else{
		compare = 2.0;
	}
	
	while (var_ != nullptr){
		if (var_->get_data(9) != compare){ //If the PE has not left (bool has turned into a double in init-function)
			//Coordinates
			x_ = var_->get_data(0) + var_->get_data(3)*dt_;
			y_ = var_->get_data(1) + var_->get_data(4)*dt_;
			z_ = var_->get_data(2) + var_->get_data(5)*dt_;
			
			//Velocities
			vx_ = var_->get_data(3) + var_->get_data(6)*dt_;
			vy_ = var_->get_data(4) + var_->get_data(7)*dt_;
			vz_ = var_->get_data(5) + var_->get_data(8)*dt_;
			
			//Acceleration
			tides->acceleration(var_->get_data(0), var_->get_data(1), var_->get_data(2), ax_, ay_, az_);
			
			eps2 = 0.5874*rh*rh;
			d = sqrt(x_*x_ + y_*y_ + z_*z_ + eps2);
			ac = -1.0*m/(d*d*d);
			ax_ = ac*x_ + ax_;
			ay_ = ac*y_ + ay_;
			az_ = ac*z_ + az_;
			
			var_->set_data(x_,0);
			var_->set_data(y_,1);
			var_->set_data(z_,2);
			var_->set_data(vx_,3);
			var_->set_data(vy_,4);
			var_->set_data(vz_,5);
			var_->set_data(ax_,6);
			var_->set_data(ay_,7);
			var_->set_data(az_,8);
		}	
		var_ = var_->get_link();
	}
	var_ = head;
	potential_escaper(N, m, rh, N_strip, tides);
}


void PotEsc::append(const double time, const double dtime, const double m, const double N, const double rh, const double E_k, TIDES* tides, const double N_strip){
	index_++;
	dt_ = dtime;
	t_ = time;
	
	//Position and velocity
	rand_position(E_k, m, rh);
	
	//Acceleration
	tides->acceleration(x_, y_, z_, ax_, ay_, az_);
	double eps2 = 0.5874*rh*rh;
	double d = sqrt(x_*x_ + y_*y_ + z_*z_ + eps2);
	double ac = -1.0*m/(d*d*d);
	ax_ = ac*x_ + ax_;
	ay_ = ac*y_ + ay_;
	az_ = ac*z_ + az_;
	
	//Initialize data in node
	double init_list[VARIABLES] = {x_, y_, z_, vx_, vy_, vz_, ax_, ay_, az_, 0.0, index_};
	var_ = new Node<double>(init_list, var_);
	potential_escaper(N, m, rh, N_strip, tides);
}


int PotEsc::count_loss(){
	Node<double>* head = var_;
	int count = 0;
	while (var_ != nullptr){
		if (var_->get_data(9) == 1.0){ //Is true if star has escaped
			count++;
		}
		var_ = var_->get_link();
	}
	var_ = head; 
	
	return count;
}

//|===================================================================================================================|

#endif //POTENTIAL_ESCAPER_H