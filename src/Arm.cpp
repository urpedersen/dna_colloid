/*
 * Arm.cpp
 *
 *  Created on: Oct 9, 2011
 *      Author: urp
 *
 */

#include "Arm.h"

//Arm::Arm(double x0,double x1,double x2, Particle * in_host_particle) {
Arm::Arm(double x0,double x1,double x2) {

	//host_particle = in_host_particle;

	length = sqrt ( x0*x0 + x1*x1 + x2*x2 ) ;

	u_local.clear();
	u_local.push_back(x0/length);
	u_local.push_back(x1/length);
	u_local.push_back(x2/length);

	u_global.clear();
	u_global.push_back(x0/length);
	u_global.push_back(x1/length);
	u_global.push_back(x2/length);

	//cout << "  " << u_local.at(0)*u_local.at(0) + u_local.at(1)*u_local.at(1) + u_local.at(2)*u_local.at(2) << endl;

	//holding_hand = null;

	// Have initially no allegiance
	allegiance = -1 ;
	possible_allegiances_connecting_particle.clear (  ) ;
	possible_allegiances_connecting_arm.clear (  ) ;
	possible_allegiances_energies.clear (  ) ;
	possible_allegiances_is_bound.clear (  ) ;

}

Arm::~Arm() {
	// Auto-generated destructor stub
}

void Arm::clear_allegiance ( ) {
	allegiance = -1 ;
	possible_allegiances_connecting_particle.clear (  ) ;
	possible_allegiances_connecting_arm.clear (  ) ;
	possible_allegiances_energies.clear (  ) ;
	possible_allegiances_is_bound.clear (  ) ;
}

void Arm::add_possible_allegiance ( double energy , int in_particle , int in_arm , bool in_is_bound ) {

	possible_allegiances_connecting_particle.push_back ( in_particle ) ;
	possible_allegiances_connecting_arm.push_back ( in_arm ) ;
	possible_allegiances_energies.push_back( energy ) ;
	possible_allegiances_is_bound.push_back( in_is_bound ) ;

	//cout << in_is_bound << endl ;

}

/**
 * Return true if arm have a connection in possible_allegiance array
 */
bool Arm::arm_have_a_bond_in_possible_allegiance(){

	bool out = false;

	for ( unsigned int i = 0 ; i < possible_allegiances_is_bound.size() ; i ++ ) {
		if( possible_allegiances_is_bound.at(i) ) out = true ;
	}

	return out;
}


/**
 * Return energy from arm from possible_allegiance arrays.
 */
double Arm::get_energy_from_possible_allegiance ( ) {

	double out = 0.0 ;

	for ( unsigned int i = 0 ; i < possible_allegiances_is_bound.size() ; i ++ ) {

		if( possible_allegiances_is_bound.at(i) ) out += possible_allegiances_energies.at(i) ;

	}

	return out ;

}

