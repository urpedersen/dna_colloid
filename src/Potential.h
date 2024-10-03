/*
 * Potential.h
 *
 *  Created on: Oct 6, 2011
 *      Author: urp
 */

#ifndef POTENTIAL_H_
#define POTENTIAL_H_

#include <iostream>
using namespace std;

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>

#include "Particle.h"

class Potential {
public:

	vector<Particle*> *particles;

	// Potential parameters
	double cut_distance ;
	double r2_cutoff    ;
	double cutoff_shift_energy ;

	//double pair_eps;   // TODO put in a parameter array
	//double pair_sigma;
	vector<double> parameters ;

	int  num_extra_sticky_spots ;
	bool extra_sticky_spots     ;
	vector<double> parameters_extra_sticky_spots;

	int only_bind_to_lowest_arm_energy;

	// Variables for order parameter
	//double r2_cutoff_orientational_order_max            ;
	//double r2_cutoff_orientational_order_min            ;
	unsigned int type_of_orientational_order              ;
	bool		 orientational_order_threshold_is_enabled ;
	double		 orientational_order_threshold_max		  ;
	double		 orientational_order_threshold_min		  ;
	unsigned int num_orientational_order_vectors          ;
	vector<vector<double> > orientational_order_vectors   ;
	bool do_orientational_order_umbrella                  ;
	double umbrella_orientational_order_kappa             ;
	double umbrella_orientational_order_center            ;

	int pin_particles;

	Potential ( ) ;

	virtual ~Potential ( ) ;

	void reinitialize ( ) ;

	void set_parameters ( string ) ;

	double get_particle_size ( ) ;
	double get_arm_eps       ( ) ;
	double get_arm_length    ( ) ;

	void print_pair_energy      (string , int )                ;
	void print_arm_energy       (string , int , int)           ;
	void print_arm_energy_angle (string , int )                ;
	void print_arm_energy_distance_angle (string , int , int ) ;
	//void print_arm_energy_reverse(string,int,int);

	double get_particle_energy  ( Particle* ) ;
	double get_internal_energy     ( ) ;
	double get_total_energy ( ) ;

	double get_pair_energy(double) ;
	double get_pair_energy(Particle*,Particle*) ;

	double get_arm_energy( double , double , double ) ;
	double get_arm_energy( double[3] , double* , Arm* ,Arm* ) ;
	double get_arm_energy( Particle* , Particle* ) ;

	double get_pinning_energy ( Particle * p ) ;
	double get_pinning_energy ( ) ;
	void   setup_pinning     ( string ) ;

	double get_orientational_order_energy ( ) ;
	void   set_parameters_orientational_order   ( string )           ;
	double get_orientational_order              ( )           ;
	double get_orientational_order              ( Particle* ) ;
	double get_orientational_order_legendre4    ( Particle* ) ;
	double get_orientational_order_legendre4_max_bonds ( unsigned int , Particle* ) ;
	double get_orientational_order_legendre4_fix_bonds ( unsigned int , Particle* ) ;
	double get_orientational_order_legendre4_neighbours ( unsigned int , Particle* ) ;
	double get_orientational_order_hexdia ( Particle* ) ;
	double get_orientational_order_legendre4_tagged ( Particle* ) ;
	double get_orientational_order_qubaticJCP115   ( Particle* ) ;

	void get_p0p1_unit_vector(double[3] , double* , Particle* , Particle* ) ;

	//void reset_particle_indexes ( ) ;
	void reset_arm_allegiance_variables ( Particle * ) ;
	//void only_select_consistant_bonds ( Particle * ) ;

};

#endif /* POTENTIAL_H_ */
