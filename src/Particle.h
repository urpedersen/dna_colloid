/*
 * Particle.h
 *
 *  Created on: Oct 6, 2011
 *      Author: urp
 */

#ifndef PARTICLE_H_
#define PARTICLE_H_

#include "Arm.h"

class Particle {
public:

	// Coordinate variables
	unsigned int index;
	unsigned int type;
	vector<double> x;	// Position vector (3 numbers between 0-1)
	vector<double> q;	// Quaternion vector (4 numbers between 0-1)

	vector<double> rotation_matrix;

	vector<int> img;	// Image Position vector (3 integers)
	vector<double> *bbox; // Size of boundary box

	vector<Arm*> arms;	// Arms on particle

	// Pair list variables
	vector<Particle*> *particles;	// All particles
	vector<Particle*> pairs;		// Pair list
	vector<Particle*> particles_bounded;

	vector<double> dx_since_last_pair_list_update;
	double *cut_distance;			// Largest interaction distance
	double *skin_distance;			// Skin that particles can move in before updating

	int *neighbour_list_updates ;

	// Variables for pinning particles in space
	//bool enable_pinning;
	vector<double> pinning_x ;
	vector<int> pinning_img  ;
	double pinning_kappa     ;

	// Constructor
	Particle(unsigned int, int,
			vector<double> *,
			double,double,double,
			int,int,int,
			double,double,double,double,
			vector<Particle*> *,
			double*,double*,
			int*
	);

	virtual ~Particle();

	// Arm
	//void reset_arm_allegiance_variable();

	// Pair list
	void make_pair_list_this();
	void make_pair_list_all();
	bool pair_list_is_old();

	// Arms
	void set_arms(string);

	// Movers
	void translate_step(int,double);

	void rotate_step(vector<double>);
	void wrap_coordinates();

	// Get positions in the minimum image
	double get_min_img(int);
	double get_min_img_x();
	double get_min_img_y();
	double get_min_img_z();

	double get_min_img(int , Arm*);
	double get_min_img_x(Arm*);
	double get_min_img_y(Arm*);
	double get_min_img_z(Arm*);

	double get_dx_img(Particle *,int);

	// Return pair distances
	double r2_mi(Particle*);
	double r_mi(Particle*);

	double r2_pinning();

	// Arms
	void get_vector_in_global_frame(vector<double>*,vector<double>*);

	// Orientations
	void update_orientation();
	void normalize_q();
	void construct_rotation_matrix();
	void set_arm_orientations();

	// Print
	void print();

};

#endif /* PARTICLE_H_ */
