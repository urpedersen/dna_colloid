/*
 * Sim.h
 *
 *  Created on: Oct 6, 2011
 *      Author: urp
 */

#ifndef BOX_H_
#define BOX_H_

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
#include "Potential.h"

class Sim {
public:

	int verbose;
	string variables;

	vector<double> bbox;

	vector<Particle*> particles;
	//int particle_counter;

	Potential potential ;

	double beta;	// inverse temperature


	double time;
	int neighbour_list_updates;

	int seed;		// For random number generator

	double skin_distance;

	// For particle insertion, grand canonical computations
	double chem_pot; 			// Chemical potential
	bool do_num_par_umbrella;
	double num_par_center;
	double num_par_kappa;
	//vector<int> num_par_histogram;

	// Variables for volume changes
	bool do_volume_moves    ;
	int type_of_volume_move ;
	double delta_L_max      ;
	double pressure         ;

	Sim ( string& ) ;

	virtual ~Sim ( ) ;

	void update                         ( ) ;
	void update_particles_indexes       ( ) ;
	void update_neighbour_list          ( ) ;
	void update_arm_allegiance_variable ( ) ;

	void add_particle ( ) ;
	void add_particle ( int ,
				        double , double , double ,
				        int    , int    , int ,
				        double , double , double , double ) ;

	void load_xyz     ( string&,string& ) ;

	void print        ( ) ;
	void print_energy ( ) ;

	string get_ener_string         ( ) ;
	string get_coordinates_xyz     ( ) ;
	string get_coordinates_network ( ) ;
	string get_coordinates_vmd     ( int ) ;
	string get_coordinates_vmd     ( int , double , double ) ;
	string get_coordinates_vmd     ( int , double , double , double , double ) ;

	double get_volume ( ) ;
	//double get_local_orientational_order ( ) ;
	double get_avg_pair_list_length ( ) ;

	int attempt_a_move_translation ( double* ) ;
	int attempt_a_move_rotation    ( double* ) ;

	int  attempt_a_move_num_particles   ( ) ;
	int  attempt_a_move_num_particles_OLD    ( ) ;
	void enable_umbrella_num_particles ( double , double ) ;

	int  attempt_a_move_volume ( ) ;
	void enable_volume_changes ( int  , double  , double  ) ;

	void wrap_coordinates ( ) ;

	void initialize_random_number ( ) ;
	double my_rand                ( ) ;
};

#endif /* BOX_H_ */
