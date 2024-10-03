/*
 * Arm.h
 *
 *  Created on: Oct 9, 2011
 *      Author: urp
 */

#ifndef ARM_H_
#define ARM_H_

#include <iostream>
using namespace std;

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>

//#include "Particle.h"

class Arm {
public:

	vector<double> u_local;
	vector<double> u_global;
	double length;

	//Particle * host_particle;

	//int my_particle_index;

	int allegiance;
	vector<int>    possible_allegiances_connecting_particle ;
	vector<int>    possible_allegiances_connecting_arm ;
	vector<bool>   possible_allegiances_is_bound ;
	vector<double> possible_allegiances_energies ;

	//Arm(double,double,double, Particle * );
	Arm ( double , double , double ) ;

	virtual ~Arm ( ) ;

	void clear_allegiance ( ) ;
	void add_possible_allegiance ( double , int , int , bool ) ;

	bool arm_have_a_bond_in_possible_allegiance ( ) ;
	double get_energy_from_possible_allegiance ( ) ;

};

#endif /* ARM_H_ */
