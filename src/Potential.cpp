/**
 * Potential.cpp
 *
 *  Created on: Oct 6, 2011
 *      Author: urp
 *
 *      TODO Add more than one hand on each arm
 *      TODO Only one holding-hands per arm.
 *
 */
#include "Potential.h"

#include <iostream>
using namespace std;

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>

#include "helpers_io.h"

#include "Particle.h"

#include "omp.h"

/**
 * Constructor
 *
 */
Potential::Potential() {
	cut_distance = 0.0;
	cutoff_shift_energy = 0.0;
}

Potential::~Potential() {
	// TODO Auto-generated destructor stub
}

/*void Potential::add_particle(Particle *input){
	particles.push_back(input);
}*/

/**
 * Run this function whenever potential parameters are changed
 */
void Potential::reinitialize(){

	r2_cutoff = cut_distance*cut_distance;

	cutoff_shift_energy = 0.0;
	cutoff_shift_energy += get_pair_energy(r2_cutoff);

}

/**
 * Set the parameter vector with parameters for the potential.
 */
void Potential::set_parameters(string variables) {

	cout << "Note: For Lennard-Jones potential, put [ pair_eps=0 ; arm_eps=0 ; pair_12_eps=4 ; pair_6_eps=-4 ] or [ pair_eps=0 ; arm_eps=0 ; pair_12_eps=1 ; pair_6_eps=-2 ]." << endl;

	parameters.clear();

	parameters.push_back(get_variable(variables,"pair_eps"  , 1.0)) ;
	parameters.push_back(get_variable(variables,"pair_sigma", 1.0)) ;

	parameters.push_back(get_variable(variables, "arm_eps"       , 1.0 ) ) ;
	parameters.push_back(get_variable(variables, "arm_r_center"  , 1.0 ) ) ;
	parameters.push_back(get_variable(variables, "arm_r_width"   , 0.2 ) ) ;

	double max_angle = get_variable( variables, "max_angle" , 40.0 )    ;
	double arm_sin_width = sin ( max_angle * 3.14159265358979 / 180.0 ) ;
	cout << "Max angle: theta_max = " << max_angle << " degree = " << max_angle * 3.14159265358979 / 180.0 << " radians." << endl;
	parameters.push_back( arm_sin_width );

	parameters.push_back(get_variable(variables,"pair_12_eps"  , 0.0));
	parameters.push_back(get_variable(variables,"pair_6_eps"  , 0.0));
	parameters.push_back(get_variable(variables,"pair_4_eps"  , 0.0));

	cout << "Note: Pair potential: pair_eps*(pair_sigma/r)^18" << endl;

	cout << "potential.set_parameters: pair_eps = "   << parameters.at(0) << endl;
	cout << "potential.set_parameters: pair_sigma = " << parameters.at(1) << endl;

	cout << "potential.set_parameters: arm_eps = "        << parameters.at(2) << endl;
	cout << "potential.set_parameters: arm_r_center = "   << parameters.at(3) << endl;
	cout << "potential.set_parameters: arm_r_width = "    << parameters.at(4) << endl;
	cout << "potential.set_parameters: arm_sin_width = "  << parameters.at(5) << endl;

	cout << "potential.set_parameters: pair_12_eps = "   << parameters.at(6) << endl;
	cout << "potential.set_parameters: pair_6_eps = "    << parameters.at(7) << endl;
	cout << "potential.set_parameters: pair_4_eps = "    << parameters.at(8) << endl;

	cut_distance = get_variable( variables , "cut_distance"   , parameters.at(3) + parameters.at(4) ) ;

	// Parameters for only allowing one bond per arm TODO write documentation for parameters only_bind_to_lowest_arm_energy
	cout << "Note: only_bind_to_lowest_arm_energy = {0,1,2} " << endl ;
	only_bind_to_lowest_arm_energy = get_variable(variables,"only_bind_to_lowest_arm_energy",0);
	//if ( int_only_bind_to_lowest_arm_energy == 0 ) {
	//	only_bind_to_lowest_arm_energy = false;
	//} else {
	//	only_bind_to_lowest_arm_energy = true;
	//}
	if ( only_bind_to_lowest_arm_energy == 0 ) {
		cout << "Note: Bind to all arms. " << endl ;
	} else if(only_bind_to_lowest_arm_energy == 2 ) {
		cout << "Note: Arms only bind to the lowest energy arm (must be lowest energy of THIS arms)." << endl;
	} else if ( only_bind_to_lowest_arm_energy == 1 ) {
		cout << "Note: Arms only bind to the lowest energy arm (just lowest energy of BOTH arms). " << endl;
	} else {
		cout << "Error: Unknown value of only_bind_to_lowest_arm_energy=" << only_bind_to_lowest_arm_energy << ". Exit 32047." << endl;
		exit(32047);
	}

	// Parameters for adding extra sticky spots.
	num_extra_sticky_spots = get_variable( variables , "num_extra_sticky_spots" , 0 ) ;
	parameters_extra_sticky_spots.clear();
	if ( num_extra_sticky_spots > 1 ) {
		cout << "Add extra sticky spots: num_extra_sticky_spots = " << num_extra_sticky_spots << "." << endl;
		extra_sticky_spots = true ;

		for ( int i = 0 ; i < num_extra_sticky_spots ; i++ ) {
			stringstream name;

			// add parameters
			name << "arm_eps" << i;
			parameters_extra_sticky_spots.push_back(
					get_variable( variables , name.str() , parameters.at(2) )
			);
			cout << "parameters_extra_sticky_spots: " << name.str() << " = " << parameters_extra_sticky_spots.back() << endl;
			name.str("");

			name << "arm_r_center" << i;
			parameters_extra_sticky_spots.push_back(
					get_variable( variables , name.str() , parameters.at(3) )
			);
			cout << "parameters_extra_sticky_spots: " << name.str() << " = " << parameters_extra_sticky_spots.back() << endl;
			name.str("");

			name << "arm_r_width" << i;
			parameters_extra_sticky_spots.push_back(
					get_variable( variables , name.str() , parameters.at(4) )
			);
			cout << "parameters_extra_sticky_spots: " << name.str() << " = " << parameters_extra_sticky_spots.back() << endl;
			name.str("");

			name << "arm_sin_width" << i;
			parameters_extra_sticky_spots.push_back(
					get_variable( variables , name.str() , parameters.at(5) )
			);
			cout << "parameters_extra_sticky_spots: " << name.str() << " = " << parameters_extra_sticky_spots.back() << endl;
			name.str("");

		}
	} else {
		cout << "Do not add extra sticky spots." << endl;
		extra_sticky_spots = false ;
	}
	reinitialize();

	// Do not pin particles as defauld
	pin_particles = 0 ;

	// Set parameters for orientaional order-parameter
	set_parameters_orientational_order( variables ) ;

}

double Potential::get_particle_size(){
	return parameters.at(1);
}


double Potential::get_arm_eps(){
	return parameters.at(2);
}


double Potential::get_arm_length(){
	return parameters.at(3);
}

/**
 * Print pair energy to file.
 */
void Potential::print_pair_energy(string filename,int number_of_points){

	FILE * ofile = fopen(filename.c_str(),"w");

	double d=cut_distance/(double)number_of_points;
	while(d<cut_distance){
		fprintf ( ofile , "%f %f\n", d , get_pair_energy(d*d) ) ;
		//cout << "pair energy: " << d << " " << get_pair_energy(d*d);
		d+=cut_distance/(double)number_of_points;
	}

	fclose(ofile);
}



/**
 * Print arm energy to file
 */
void Potential::print_arm_energy(string filename,int number_of_points,int number_of_points_angle){
	FILE * ofile = fopen(filename.c_str(),"w");

	for ( double angle  = 0.0 ;
		         angle  < 3.14159265358979 ;
		         angle += 3.14159265358979/number_of_points  ) {

		fprintf ( ofile , "# cos = %f \n" , cos(angle) ) ;

		for ( double d = 0.0 ; d <= cut_distance ; d += cut_distance/(double)number_of_points ) {

			fprintf ( ofile , "%f %f\n", d , get_arm_energy(d*d,cos(angle),1.0) ) ;

		}

		fprintf ( ofile , "\n" ) ;
	}

	fclose(ofile);

}


void Potential::print_arm_energy_angle( string filename , int number_of_points ) {

	FILE * ofile = fopen(filename.c_str(),"w");

	double r_min = parameters.at(3) ;

	fprintf ( ofile , "# [ angle , angle energy at r = %f ]\n",r_min);

	for ( double angle  = 0.0 ;
		         angle  < 3.14159265358979 ;
		         angle += 3.14159265358979/number_of_points  ) {
		fprintf( ofile , "%f %f\n" , angle , get_arm_energy( r_min * r_min , cos(angle) , 1.0 ) );
	}

	fclose(ofile);

}

void Potential::print_arm_energy_distance_angle( string filename , int number_of_points_distance , int number_of_points_angle ) {

	FILE * ofile = fopen(filename.c_str(),"w");

	//double r_min = parameters.at(3) ;
	//fprintf ( ofile , "# [ angle , angle energy at r = %f ]\n",r_min);
	fprintf ( ofile , "\n" ) ;

	for ( double d = 0.0 ; d <= cut_distance ; d += cut_distance/(double)number_of_points_distance ) {
		for ( double angle  = 0.0 ;
					 angle  < 3.14159265358979 ;
					 angle += 3.14159265358979/number_of_points_angle  ) {
			fprintf( ofile , "%f %f %f\n" , d, angle/3.14159265358979*180.0 , get_arm_energy( d * d , cos(angle) , 1.0 ) );
		}
	}
	fclose(ofile);
}



/**
 * Return energy of particle.
 */
double Potential::get_particle_energy (Particle *p) {
	double out = 0.0;


	for (unsigned int i = 0 ; i < p->pairs.size() ; i++ ) {

		out += get_pair_energy ( p , p -> pairs.at(i) ) ;

		//out += get_arm_energy  ( p , p -> pairs.at(i) ) ;
		//cout << get_arm_energy  ( p , p -> pairs.at(i) )  << endl;
	}

	// Reset reset_arm_allegiance_variables for this particle, and its neighbours
	reset_arm_allegiance_variables ( p ) ;

	//if ( only_bind_to_lowest_arm_energy == 1 ) {
		//for(unsigned int i = 0 ; i < p->pairs.size() ; i++ ) {
			//cout << i << endl;
			//reset_arm_allegiance_variables ( p -> pairs.at(i) );
		//}
	//}
	//only_select_consistant_bonds( p );

	// Do special thing, of particle
	/*if ( only_bind_to_lowest_arm_energy == 1 ) {

		double lowest_energy = 0.0;

		for(unsigned int ia = 0 ; ia < p->arms.size() ; ia++ ){
			for ( unsigned int i = 0 ; i < p->arms.at(ia)->possible_allegiances_connecting_arm.size() ; i++ ) {
				if( p->arms.at(ia)->possible_allegiances_is_bound.at(i) ){
					double next_energy = p->arms.at(ia)->possible_allegiances_energies.at(i);
					if( next_energy < lowest_energy ){
						lowest_energy = next_energy ;
					}
				}
			}
		}

		out += lowest_energy ;

	} else {*/
	    // TODO add openMP in parallel
	    //omp_set_num_threads(p->arms.size());
        #pragma omp for reduction(+:out)
		for(unsigned int ia = 0 ; ia < p->arms.size() ; ia++ ){
			for ( unsigned int i = 0 ; i < p->arms.at(ia)->possible_allegiances_connecting_arm.size() ; i++ ) {
				if( p->arms.at(ia)->possible_allegiances_is_bound.at(i) ){
					out += p->arms.at(ia)->possible_allegiances_energies.at(i);
				}
			}
		}

	//}

	return out;
}

/**
 * Return total energy (without umbrella energy)
 */
double Potential::get_internal_energy ( ) {
	double out = 0.0;
	// Sum over all pairs
	for(unsigned int i = 0 ; i < particles->size() ; i++ ){
		out += 0.5*get_particle_energy(particles->at(i));
	}
	return out;
}

/**
 * Return total energy INCLUDING orientationa order umbrella energy
 */
double Potential::get_total_energy ( ) {

	return get_internal_energy()
	     + get_orientational_order_energy()
	     + get_pinning_energy() ;

}


/**
 * Return pair-energy
 *
 * Input: Squared pair distance, [pair distance]^2.
 *
 */
double Potential::get_pair_energy( double r2 ) {

	if( r2 <= r2_cutoff ) {

		double *pair_eps    = &parameters.at(0) ;
		double *pair_sigma  = &parameters.at(1) ;

		double *pair_12_eps = &parameters.at(6) ;
		double *pair_6_eps  = &parameters.at(7) ;
		double *pair_4_eps  = &parameters.at(8) ;

		double r2inv  = 1 / r2 ;
		       r2inv *= (*pair_sigma) * (*pair_sigma) ;

		double r4inv  = r2inv * r2inv ;

		double r6inv  = r2inv * r2inv * r2inv ;

		double r12inv = r6inv * r6inv         ;
		double r18inv = r6inv * r6inv * r6inv ;

		return (   (*pair_eps   ) * r18inv
				 + (*pair_12_eps) * r12inv
				 + (*pair_6_eps ) * r6inv
				 + (*pair_4_eps ) * r4inv   ) - cutoff_shift_energy ;

	} else {

		return 0.0 ;

	}

}

/**
 * Return pair energy
 */
double Potential::get_pair_energy( Particle *p0 , Particle *p1 ) {

	double energy = 0.0;

	double r2 = p0 -> r2_mi ( p1 ) ;
	energy += get_pair_energy ( r2 ) ;

	return energy ;

}



/**
 * Return arm energy:
 *
 *  u = eps * [ (r - r0)^2/r_w^2 + ( sin^2(theta0) + sin^2(theta1) ) / s_w^2 - 1]
 *        if (r - r0)^2/r_w^2 + sin^2(theta0)/s_w^2 < 1, and
 *  u = 0
 *        otherwise.
 *
 *  Note: input takes cos(theta)
 *     and function use that sin^2(theta) = 1 - cos^2(theta)
 *
 *
 * Input:
 * 	(double) Squared pair distance, [pair distance]^2.
 *  (double) cos(theta),
 *
 *      Note: Let cos carry a sign
 *      so that anti-parallel vectors have negative cos2.
 *      Code use that sin^2 = 1 - cos^2.
 *
 */
double Potential::get_arm_energy( double r2 , double cos0 , double cos1 ) {

	if( r2 <= r2_cutoff and cos0 > 0.0 and cos1 > 0.0 ) {

		if( extra_sticky_spots ) {

			double sum = 0.0 ;

			for ( int spot = 0 ; spot < num_extra_sticky_spots ; spot++ ) {
				double arm_eps       = parameters_extra_sticky_spots.at(4*spot+0) ;
				double arm_r_center  = parameters_extra_sticky_spots.at(4*spot+1) ;
				double arm_r_width   = parameters_extra_sticky_spots.at(4*spot+2) ;
				double arm_sin_width = parameters_extra_sticky_spots.at(4*spot+3) ;

				double sin0sin0 = ( 1 - cos0 * cos0 ) ;
				double sin1sin1 = ( 1 - cos1 * cos1 ) ;

				double d2 = 0.0;

				d2 += sin0sin0 + sin1sin1;
				d2 /= ( arm_sin_width * arm_sin_width ) ;

				if ( d2 < 1.0 ) {

					double delta_r = ( sqrt(r2) - arm_r_center ) / arm_r_width;
					d2 += delta_r*delta_r;

					if ( d2 < 1.0 ) {

						sum += arm_eps * ( d2 - 1.0 );

					} else {

						// sum += 0.0 ;

					}

				} else {

					// sum += 0.0 ;

				}
			}

			return sum ;

		}else{ // Only one sticky spot

			double arm_eps = parameters.at(2) ;
			double arm_r_center = parameters.at(3) ;
			double arm_r_width  = parameters.at(4) ;
			double arm_sin_width = parameters.at(5) ;

			double sin0sin0 = ( 1 - cos0 * cos0 ) ;
			double sin1sin1 = ( 1 - cos1 * cos1 ) ;

			double d2 = 0.0;

			d2 += sin0sin0 + sin1sin1;
			d2 /= ( arm_sin_width * arm_sin_width ) ;

			/*
			cout << "cos0 = " << cos0;
			cout << " cos1 = " << cos1;
			cout << endl;
			*/

			//cout << sin0sin0 << " " << sin1sin1 << " " <<  delta_r << " " << d2 << endl;

			if ( d2 < 1.0 ) {

				double delta_r = ( sqrt(r2) - arm_r_center ) / arm_r_width;
				d2 += delta_r*delta_r;

				if ( d2 < 1.0 ) {

					return arm_eps * ( d2 - 1.0 );

				} else {

					return 0.0;

				}

			} else {

				return 0.0 ;

			}
		}

	} else {

		return 0.0 ;

	}

}



/**
 * Return arm energy.
 */
double Potential::get_arm_energy( double* p0p1_unit_vector , double *r2 , Arm *a0 , Arm *a1) {

	if( *r2 <= r2_cutoff ) {

		/*
		vector<double> arm_unit_vector0;
		arm_unit_vector0.assign(3,0.0);
		p0 -> get_vector_in_global_frame ( &arm_unit_vector0 , &(a0->u_local) );

		vector<double> arm_unit_vector1;
		arm_unit_vector1.assign(3,0.0);
		p1 -> get_vector_in_global_frame ( &arm_unit_vector1 , &(a1->u_local) );
		*/

		double cos0 = 0.0;
		double cos1 = 0.0;
		for (int i = 0 ; i < 3 ; i++ ) cos0 += a0->u_global.at(i)          * p0p1_unit_vector[i] ;
		for (int i = 0 ; i < 3 ; i++ ) cos1 += a1->u_global.at(i) * (-1.0) * p0p1_unit_vector[i] ;

		/*cout << " p0 = " << p0->index;
		cout << " v01 = " << p0p1_unit_vector.at(0) << " " << p0p1_unit_vector.at(1) << " "  << p0p1_unit_vector.at(2) ;
		cout << " v_arm0 = " << arm_unit_vector0.at(0) << " " << arm_unit_vector0.at(1) << " "  << arm_unit_vector0.at(2) ;
		cout << " v_arm1 = " << arm_unit_vector1.at(0) << " " << arm_unit_vector1.at(1) << " "  << arm_unit_vector1.at(2) ;
		cout << " cos0 = " << cos0;
		cout << " cos1 = " << cos1;
		cout << endl;*/

		return get_arm_energy ( *r2 , cos0 , cos1 ) ;

	} else {

		return 0.0;

	}

}






/**
 * Return arm energy
 */
double Potential::get_arm_energy ( Particle *p0 , Particle *p1 ) {

	double p0p1_unit_vector[] = { p0->get_dx_img(p1,0) ,
		                          p0->get_dx_img(p1,1) ,
		                          p0->get_dx_img(p1,2) } ;
	double r2 = 0.0 ;
	for ( int i = 0 ; i < 3 ; i++ ) r2+= p0p1_unit_vector[i]*p0p1_unit_vector[i];
	double r = sqrt(r2);
	for ( int i = 0 ; i < 3 ; i++ ) p0p1_unit_vector[i] /= r ;

	/*cout << p0->index << " : " ;
	for (int d = 0 ; d < 3 ; d++ )
		cout << p0->x.at(d) << " " ;
	cout << "  " << p1->index << " : " ;
	for (int d = 0 ; d < 3 ; d++ )
		cout << p1->x.at(d) << " " ;
	cout << "  r2 = " << r2 ;
	cout << endl; */

	double energy = 0.0 ;

	for ( unsigned int i = 0  ; i < p0-> arms.size() ; i++ ) {

		for ( unsigned int j = 0 ; j < p1->arms.size() ; j++ ) {

			energy += get_arm_energy ( p0p1_unit_vector , &r2 , p0->arms.at(i) , p1->arms.at(j) );

		}

	}

	return energy ;

}

/**
 * Return the energy of the particle pinning of particle p
 */
double Potential::get_pinning_energy ( Particle * p ) {

	if ( pin_particles == 0 ) {

		return 0.0 ;

	} else if ( pin_particles == 1 ) {

		return 0.5 * p->pinning_kappa * p->r2_pinning();

	} else {

		cout << "Error: Unknown value of pin_particle=" <<  pin_particles << ". Exit 22585." <<endl;
		exit(22585);
		return 0.0 ;

	}
}

/**
 * Return the energy of the particle pinning all particles
 */
double Potential::get_pinning_energy ( ) {

	double out = 0.0 ;

	for(unsigned int i = 0 ; i < particles->size() ; i++ )
		out += get_pinning_energy ( particles->at(i) ) ;

	return out ;

}

/**
 * Enable a pinning in space potential of particles
 */
void Potential::setup_pinning ( string variables ) {

	pin_particles = get_variable( variables ,"pin_particles" , 0 ) ;

	if ( pin_particles == 0 ) {

		cout << "Note: Do not pin particles in space." << endl;

	} else if ( pin_particles == 1 ) {

		cout << "Note: Pin particles in space is enabled." << endl;

		double pin_kappa = get_variable( variables ,"pin_kappa" , 1.0 ) ;
		cout << "Note:     ... all particles are pinned to a harmonic-pin with spring constant of pin_kappa=" << pin_kappa << ". " << endl;

		for( unsigned int i = 0 ; i < particles->size() ; i++ )
			particles->at(i)->pinning_kappa = pin_kappa ;

	} else {

		cout << "Error: Unknown value of pin_particle=" <<  pin_particles << ". Exit 1845." <<endl;
		exit(1845);

	}
}

/**
 * Return (umbrella) energy of orientational order
 */
double Potential::get_orientational_order_energy ( ) {

	if(do_orientational_order_umbrella){
		double dQ = get_orientational_order() -  umbrella_orientational_order_center ;
		return 0.5*umbrella_orientational_order_kappa*dQ*dQ;
	} else {
		return 0.0 ;
	}
}

/**
 * Initializes parameters for the orientational order parameter
 */
void Potential::set_parameters_orientational_order(string variables){

	cout << "Note:   Selections of order parameters:                                                       " << endl ;
	cout << "Note:   type_of_orientational_order=0    ->    Do not calculate order parameter               " << endl ;
	cout << "Note:   type_of_orientational_order=1    ->    Use legendre4                                  " << endl ;
	cout << "Note:   type_of_orientational_order=2    ->    Use legendre4_tagged      orientation relative to particle 0 " << endl   ;
	cout << "Note:   type_of_orientational_order=3    ->    Use qubatic               as in [JPC-B (2011) 115,14205]       " << endl ;
	cout << "Note:   type_of_orientational_order=4    ->    Use legendre4_max_bonds   with Nb=4           " << endl ;
	cout << "Note:   type_of_orientational_order=6    ->    Use ...                   with Nb=6           " << endl ;
	cout << "Note:   type_of_orientational_order=8    ->    Use ...                   with Nb=8           " << endl ;
	cout << "Note:   type_of_orientational_order=14   ->    Use legendre4_fix_bonds   with Nb=4           " << endl ;
	cout << "Note:   type_of_orientational_order=16   ->    Use ...                   with Nb=6           " << endl ;
	cout << "Note:   type_of_orientational_order=18   ->    Use ...                   with Nb=8           " << endl ;
	cout << "Note:   type_of_orientational_order=24   ->    Use legendre4_neighbours  with Nb=4           " << endl ;
	cout << "Note:   type_of_orientational_order=26   ->    Use ...                   with Nb=6           " << endl ;
	cout << "Note:   type_of_orientational_order=28   ->    Use ...                   with Nb=8           " << endl ;
	cout << "Note:   type_of_orientational_order=106  ->    Use hexdia (optimized for hexagonal diamond)  " << endl ;

	type_of_orientational_order = get_variable( variables ,"type_of_orientational_order" , 1 ) ;

	//r2_cutoff_orientational_order_max = get_variable( variables ,"r2_cutoff_orientational_order_max" , r2_cutoff )       ;
	//r2_cutoff_orientational_order_min = get_variable( variables ,"r2_cutoff_orientational_order_min" , r2_cutoff / 4.0 ) ;

	// Make reference vectors along the Cartesian coordinates.
	orientational_order_vectors.clear();
	num_orientational_order_vectors = 3 ;
	for(unsigned int i = 0 ; i<num_orientational_order_vectors ; i++ ) {
		vector<double> tmp( 3 , 0.0 ) ;
		tmp.at(i) = 1.0 ;
		orientational_order_vectors.push_back(tmp) ;
		//cout << tmp.at(0) << " " << tmp.at(1) << " " << tmp.at(2) << endl;
	}

	// Set a treshold so that
	cout << "Note: If orientational_order_threshold is enables (1) then the order is unity if it is within the threshold value and zero otherwise." << endl;
	int enable_orientational_order_threshold = get_variable( variables ,"enable_orientational_order_threshold" , 0 ) ;
	if(enable_orientational_order_threshold > 0 ) {
		orientational_order_threshold_is_enabled = true ;
	} else {
		orientational_order_threshold_is_enabled = false ;
	}
	orientational_order_threshold_max = get_variable( variables ,"orientational_order_threshold_max" , 1.25 ) ;
	orientational_order_threshold_min = get_variable( variables ,"orientational_order_threshold_min" , 0.75 ) ;

	umbrella_orientational_order_kappa = get_variable( variables ,"umbrella_orientational_order_kappa" , 0.0 ) ;
	umbrella_orientational_order_center = get_variable( variables ,"umbrella_orientational_order_center" , 0.0 ) ;
	if( umbrella_orientational_order_kappa>0.0 ) {
		do_orientational_order_umbrella = true  ;
		cout << "Note: Orientational order umbrella is enabled (kappa>0). " << endl   ;
	} else {
		do_orientational_order_umbrella = false ;
		cout << "Note: Orientational order umbrella is disabled (kappa<=0). " << endl ;
	}
}

/**
 * Return local orientational order parameter
 */
double Potential::get_orientational_order ( ) {

	double out = 0.0 ;

	for ( unsigned int i = 0 ; i < particles->size() ; i++ ) {

		out += get_orientational_order( particles->at(i) );

	}

	out /= (double) particles->size() ;

	return out ;
}

double Potential::get_orientational_order ( Particle *p0 ) {

	double val = 0.0 ;

	if( type_of_orientational_order==0 ) {
		return 0.0 ;
	} else if ( type_of_orientational_order==1 ) {
		val = get_orientational_order_legendre4 ( p0 ) ;
	} else if ( type_of_orientational_order==2 ) {
		val = get_orientational_order_legendre4_tagged ( p0 ) ;
	} else if ( type_of_orientational_order==3 ) {
		val = get_orientational_order_qubaticJCP115 ( p0 ) ;
	} else if ( type_of_orientational_order==4 ) {
		val = get_orientational_order_legendre4_max_bonds ( 4 , p0 ) ;
	} else if ( type_of_orientational_order==6 ) {
		val = get_orientational_order_legendre4_max_bonds ( 6 , p0 ) ;
	} else if ( type_of_orientational_order==8 ) {
		val = get_orientational_order_legendre4_max_bonds ( 8 , p0 ) ;
	} else if ( type_of_orientational_order==14 ) {
		val = get_orientational_order_legendre4_fix_bonds ( 4 , p0 ) ;
	} else if ( type_of_orientational_order==16 ) {
		val = get_orientational_order_legendre4_fix_bonds ( 6 , p0 ) ;
	} else if ( type_of_orientational_order==18 ) {
		val = get_orientational_order_legendre4_fix_bonds ( 8 , p0 ) ;
	} else if ( type_of_orientational_order==24 ) {
		val = get_orientational_order_legendre4_neighbours ( 4 , p0 ) ;
	} else if ( type_of_orientational_order==26 ) {
		val = get_orientational_order_legendre4_neighbours ( 6 , p0 ) ;
	} else if ( type_of_orientational_order==28 ) {
		val = get_orientational_order_legendre4_neighbours ( 8 , p0 ) ;
	} else if ( type_of_orientational_order==106 ) {
		val = get_orientational_order_hexdia ( p0 ) ;
	} else {
		cout << "Error: unknown value of type_of_orientational_order=" << type_of_orientational_order << ". Exit." << endl;
		exit(12660);
		return 0.0 ;
	}


	if(orientational_order_threshold_is_enabled){
		if(val>orientational_order_threshold_min && val<orientational_order_threshold_max) {
			return 1.0 ;
		} else {
			return 0.0 ;
		}
	}else{
		return val;
	}

}

/**
 * Return local orientational order parameter of particle.
 *
 *   (cos^2(2*theta)-7/15)/(1-7/15)
 *
 *  Notes:       7/15  = 0.466666666666667
 *        1/(1 - 7/15) = 1.875
 */
double Potential::get_orientational_order_legendre4 ( Particle *p0 ) {

	//cout << sqrt(r2_cutoff_order_max) << " " << sqrt(r2_cutoff_order_min) << endl;

	int num_vector_pairs = 0    ;
	double sum = 0.0            ;

	vector<double> vector0( 3 , 0.0 ) ;
	vector<double> vector1( 3 , 0.0 ) ;

	for(unsigned int i = 0 ; i < p0->particles_bounded.size() ; i++ ) {

		Particle *p1 = p0 -> particles_bounded.at(i) ;

		//cout << p0->index << " - " <<  p1->index << ": ";

		//double r2 = p0 -> r2_mi ( p1 ) ;

		for ( unsigned int iv0 = 0 ; iv0<orientational_order_vectors.size() ; iv0++ ) {
		//{   unsigned int iv0 = 0 ;

			p0->get_vector_in_global_frame ( &vector0 , &orientational_order_vectors.at(iv0) ) ;

			for ( unsigned int iv1 = 0 ; iv1<orientational_order_vectors.size() ; iv1++ ){

				num_vector_pairs++;

				p1->get_vector_in_global_frame ( &vector1 , &orientational_order_vectors.at(iv1) ) ;

				// Note: it is assumed that that vector0 and vector1 are unit vectors;
				double cos_theta = 0.0 ;
				for (int i = 0 ; i < 3 ; i++) cos_theta+=vector0.at(i)*vector1.at(i) ;

				/*double theta;
				if(cos_theta>1.0 && cos_theta<1.00001 ){   // Take care of round-off error
					theta = 0.0 ;
				} else if ( cos_theta>-1.00001 && cos_theta<-1.0 ) {
					theta = 3.14159265358979 ;
				} else {
					theta = acos(cos_theta) ;
				}*/

				// Fourth legencre polynomial P_4(x=cos(theta)) = 1/8( 35 x^4 - 30 x^2 + 3 )
				double cos_theta2=cos_theta*cos_theta   ;
				double cos_theta4=cos_theta2*cos_theta2 ;
				sum += 4.375*cos_theta4 - 3.75*cos_theta2 + 0.375;
				//cout << cos_theta << " ";
				//cout << "cos_theta = " <<  cos_theta  << " theta = " << theta << endl ;
			}
		}
		//cout << endl ;
	}

	if( num_vector_pairs > 0 ) {
		// scale_factor = ((1+2*P_4(0))/3)^-1 = 1.71428571428571
		return sum/(double)num_vector_pairs*1.71428571428571 ;
	} else {
		return 0.0 ;
	}
}

/**
 * Return local orientational order parameter of particle.
 *
 *   max_bonds   -> mnormalized by max_bonds and zero if number of bonds is larger than max_bonds
 */
double Potential::get_orientational_order_legendre4_max_bonds ( unsigned int num_max_bonds , Particle *p0 ) {

	//cout << sqrt(r2_cutoff_order_max) << " " << sqrt(r2_cutoff_order_min) << endl;

	int num_vector_pairs = 0    ;
	double sum = 0.0            ;

	vector<double> vector0( 3 , 0.0 ) ;
	vector<double> vector1( 3 , 0.0 ) ;

	for(unsigned int i = 0 ; i < p0->particles_bounded.size() ; i++ ) {

		Particle *p1 = p0 -> particles_bounded.at(i) ;

		//cout << p0->index << " - " <<  p1->index << ": ";

		//double r2 = p0 -> r2_mi ( p1 ) ;

		for ( unsigned int iv0 = 0 ; iv0<orientational_order_vectors.size() ; iv0++ ) {
		//{   unsigned int iv0 = 0 ;

			p0->get_vector_in_global_frame ( &vector0 , &orientational_order_vectors.at(iv0) ) ;

			for ( unsigned int iv1 = 0 ; iv1<orientational_order_vectors.size() ; iv1++ ){

				num_vector_pairs++;

				p1->get_vector_in_global_frame ( &vector1 , &orientational_order_vectors.at(iv1) ) ;

				// Note: it is assumed that that vector0 and vector1 are unit vectors;
				double cos_theta = 0.0 ;
				for (int i = 0 ; i < 3 ; i++) cos_theta+=vector0.at(i)*vector1.at(i) ;

				/*double theta;
				if(cos_theta>1.0 && cos_theta<1.00001 ){   // Take care of round-off error
					theta = 0.0 ;
				} else if ( cos_theta>-1.00001 && cos_theta<-1.0 ) {
					theta = 3.14159265358979 ;
				} else {
					theta = acos(cos_theta) ;
				}*/

				// Fourth legencre polynomial P_4(x=cos(theta)) = 1/8( 35 x^4 - 30 x^2 + 3 )
				double cos_theta2=cos_theta*cos_theta   ;
				double cos_theta4=cos_theta2*cos_theta2 ;
				sum += 4.375*cos_theta4 - 3.75*cos_theta2 + 0.375;
				//cout << cos_theta << " ";
				//cout << "cos_theta = " <<  cos_theta  << " theta = " << theta << endl ;
			}
		}
		//cout << endl ;
	}

	if( num_vector_pairs > 0 && p0->particles_bounded.size()<=num_max_bonds ) {
		// scale_factor = ((1+2*P_4(0))/3)^-1 = 1.71428571428571
		return sum/(double)num_vector_pairs*1.71428571428571*(p0->particles_bounded.size()/(double)num_max_bonds);
	} else {
		return 0.0 ;
	}
}


/**
 * Return local orientational order parameter of particle.
 *
 *  fix_bonds   ->  Must have num_bonds connections to be non-zero.
 */
double Potential::get_orientational_order_legendre4_fix_bonds ( unsigned int num_bonds , Particle *p0 ) {

	//cout << sqrt(r2_cutoff_order_max) << " " << sqrt(r2_cutoff_order_min) << endl;

	int num_vector_pairs = 0    ;
	double sum = 0.0            ;

	vector<double> vector0( 3 , 0.0 ) ;
	vector<double> vector1( 3 , 0.0 ) ;

	if( p0->particles_bounded.size() == num_bonds ) {
		for(unsigned int i = 0 ; i < p0->particles_bounded.size() ; i++ ) {

			Particle *p1 = p0 -> particles_bounded.at(i) ;

			//cout << p0->index << " - " <<  p1->index << ": ";

			//double r2 = p0 -> r2_mi ( p1 ) ;

			//double val = 0.0 ;

			for ( unsigned int iv0 = 0 ; iv0<orientational_order_vectors.size() ; iv0++ ) {
			//{   unsigned int iv0 = 0 ;

				p0->get_vector_in_global_frame ( &vector0 , &orientational_order_vectors.at(iv0) ) ;

				for ( unsigned int iv1 = 0 ; iv1<orientational_order_vectors.size() ; iv1++ ) {

					num_vector_pairs++ ;

					p1->get_vector_in_global_frame ( &vector1 , &orientational_order_vectors.at(iv1) ) ;

					// Note: it is assumed that that vector0 and vector1 are unit vectors;
					double cos_theta = 0.0 ;
					for (int i = 0 ; i < 3 ; i++) cos_theta+=vector0.at(i)*vector1.at(i) ;

					/*double theta;
					if(cos_theta>1.0 && cos_theta<1.00001 ){   // Take care of round-off error
						theta = 0.0 ;
					} else if ( cos_theta>-1.00001 && cos_theta<-1.0 ) {
						theta = 3.14159265358979 ;
					} else {
						theta = acos(cos_theta) ;
					}*/

					// Fourth legencre polynomial P_4(x=cos(theta)) = 1/8( 35 x^4 - 30 x^2 + 3 )
					double cos_theta2=cos_theta*cos_theta   ;
					double cos_theta4=cos_theta2*cos_theta2 ;
					sum += 4.375*cos_theta4 - 3.75*cos_theta2 + 0.375;
					//val += 4.375*cos_theta4 - 3.75*cos_theta2 + 0.375;
					//cout << cos_theta << " ";
					//cout << 4.375*cos_theta4 - 3.75*cos_theta2 + 0.375 << " ";
					//cout << "cos_theta = " <<  cos_theta  << " theta = " << theta << endl ;
				}
			}
			//cout << endl << (-val/9.0)*3.5 << " " << val*1.71428571428571/9.0 << endl;
			//cout << endl ;
		}

		//cout << "           " << sum/(double)num_vector_pairs <<  endl ;

		return sum/(double)num_vector_pairs*1.71428571428571;

	} else {

		return 0.0 ;

	}

	/*if( num_vector_pairs > 0 && p0->particles_bounded.size()<=num_max_bonds ) {
		// scale_factor = ((1+2*P_4(0))/3)^-1 = 1.71428571428571
	} else {

	}*/
}


/**
 * Return local orientational order parameter of particle.
 *
 *  Sum over all orientations in neighbor cluster
 *
 *   max_bonds   -> mnormalized by max_bonds and zero if number of bonds is larger than max_bonds
 */
double Potential::get_orientational_order_legendre4_neighbours ( unsigned int num_max_bonds , Particle *p_center ) {

	//cout << sqrt(r2_cutoff_order_max) << " " << sqrt(r2_cutoff_order_min) << endl;

	int num_vector_pairs = 0    ;
	double sum = 0.0            ;

	vector<double> vector0( 3 , 0.0 ) ;
	vector<double> vector1( 3 , 0.0 ) ;

	if ( p_center->particles_bounded.size()>0 && p_center->particles_bounded.size()<=num_max_bonds ) {

		for(int j = -1 ; j < (int)p_center->particles_bounded.size()-1 ; j++ ) {

			Particle *p0      ;
			if ( j==-1 ) {
				p0 = p_center ;
			} else {
				p0 = p_center -> particles_bounded.at(j) ;
			}

			for(unsigned int i = j+1 ; i < p_center->particles_bounded.size() ; i++ ) {

				//cout << "  wee    " << p_center->particles_bounded.size() << endl ;


				Particle *p1 = p_center -> particles_bounded.at(i) ;

				//cout << p0->index << " - " <<  p1->index << ": ";

				//double r2 = p0 -> r2_mi ( p1 ) ;

				for ( unsigned int iv0 = 0 ; iv0<orientational_order_vectors.size() ; iv0++ ) {
				//{   unsigned int iv0 = 0 ;

					p0->get_vector_in_global_frame ( &vector0 , &orientational_order_vectors.at(iv0) ) ;

					for ( unsigned int iv1 = 0 ; iv1<orientational_order_vectors.size() ; iv1++ ){

						num_vector_pairs++;

						p1->get_vector_in_global_frame ( &vector1 , &orientational_order_vectors.at(iv1) ) ;

						// Note: it is assumed that that vector0 and vector1 are unit vectors;
						double cos_theta = 0.0 ;
						for (int i = 0 ; i < 3 ; i++) cos_theta+=vector0.at(i)*vector1.at(i) ;

						/*double theta;
						if(cos_theta>1.0 && cos_theta<1.00001 ){   // Take care of round-off error
							theta = 0.0 ;
						} else if ( cos_theta>-1.00001 && cos_theta<-1.0 ) {
							theta = 3.14159265358979 ;
						} else {
							theta = acos(cos_theta) ;
						}*/

						// Fourth legencre polynomial P_4(x=cos(theta)) = 1/8( 35 x^4 - 30 x^2 + 3 )
						double cos_theta2=cos_theta*cos_theta   ;
						double cos_theta4=cos_theta2*cos_theta2 ;
						sum += 4.375*cos_theta4 - 3.75*cos_theta2 + 0.375;
						//cout << cos_theta << " ";
						//cout << "cos_theta = " <<  cos_theta  << " theta = " << theta << endl ;
					}
				}
				//cout << endl ;
			}
		}
	}

	// Count maximum number of bonds in neighbour cluster for normalization
	unsigned int num_max_bonds_input = num_max_bonds ;
	num_max_bonds = 0 ;
	for(unsigned int i = num_max_bonds_input ; i > 0 ; i-- ){
		num_max_bonds += i ;
	}

	if( num_vector_pairs > 0 && p_center->particles_bounded.size()<=num_max_bonds_input ) {
		// scale_factor = ((1+2*P_4(0))/3)^-1 = 1.71428571428571
		return sum/(double)num_vector_pairs*1.71428571428571*(p_center->particles_bounded.size()/(double)num_max_bonds);
	} else {
		return 0.0 ;
	}

}



/**
 * Return local orientational order parameter of particle.
 *
 *  Optimized for hexagonal diamond
 *
 *  In the hexagonal diamond structure three bonding particles have `qubatic` cartesian vectors.
 *  The fourth Legendre polynomial (L4) is used to detect this.
 *
 *  The fourth particle, however, does not have this alignment.
 *  The sum of the L4's for this particle is negative.
 *
 *  Thus for a orderes particles we expect to find one with a negative value.
 *  The vectors is renormalised so a value of one is gained for perfect hexagonal diamond structure.
 *
 */
double Potential::get_orientational_order_hexdia ( Particle *p0 ) {

	int num_vector_pairs = 0    ;
	double sum = 0.0            ;

	unsigned int num_bonds = 4  ;

	vector<double> vector0( 3 , 0.0 ) ;
	vector<double> vector1( 3 , 0.0 ) ;

	int num_negative_vals = 0 ;

	if( p0->particles_bounded.size() <= num_bonds && p0->particles_bounded.size() > 0 ) {

		for(unsigned int i = 0 ; i < p0->particles_bounded.size() ; i++ ) {

			Particle *p1 = p0 -> particles_bounded.at(i) ;

			//cout << p0->index << " - " <<  p1->index << ": ";

			double val = 0.0 ;

			for ( unsigned int iv0 = 0 ; iv0<orientational_order_vectors.size() ; iv0++ ) {

				p0->get_vector_in_global_frame ( &vector0 , &orientational_order_vectors.at(iv0) ) ;

				for ( unsigned int iv1 = 0 ; iv1<orientational_order_vectors.size() ; iv1++ ) {

					num_vector_pairs++ ;

					p1->get_vector_in_global_frame ( &vector1 , &orientational_order_vectors.at(iv1) ) ;

					// Note: it is assumed that that vector0 and vector1 are unit vectors;
					double cos_theta = 0.0 ;
					for (int i = 0 ; i < 3 ; i++) cos_theta+=vector0.at(i)*vector1.at(i) ;

					// Fourth legencre polynomial P_4(x=cos(theta)) = 1/8( 35 x^4 - 30 x^2 + 3 )
					double cos_theta2=cos_theta*cos_theta   ;
					double cos_theta4=cos_theta2*cos_theta2 ;
					val += 4.375*cos_theta4 - 3.75*cos_theta2 + 0.375;

				}
			}

			if ( val<0.0 ) {
				num_negative_vals++ ;
				sum+=-val*3.5 ;
			} else {
				sum+=val*1.71428571428571;
			}

		}

		//cout << "           " << sum/(double)num_vector_pairs <<  endl ;
		if ( (p0->particles_bounded.size()  < 4 && num_negative_vals  < 2  ) ||
			 (p0->particles_bounded.size() == 4 && num_negative_vals == 1  )
		) {
			return sum/(double)num_vector_pairs*p0->particles_bounded.size()/num_bonds;
		} else {
			return 0.0 ;
		}

	} else {

		return 0.0 ;

	}
}


/**
 * Return local orientational order parameter of particle.
 *
 *   _tagged    ->    Calculate the legendre4 relative to the tagged particle
 */
double Potential::get_orientational_order_legendre4_tagged ( Particle *p0 ) {

	//cout << sqrt(r2_cutoff_order_max) << " " << sqrt(r2_cutoff_order_min) << endl;

	int num_vector_pairs = 0    ;
	double sum = 0.0            ;

	vector<double> vector0( 3 , 0.0 ) ;
	vector<double> vector1( 3 , 0.0 ) ;

	if ( p0->particles_bounded.size() > 0 ) {

		Particle *p1 = particles->at(0) ;

		//cout << p0->index << " - " <<  p1->index << ": ";

		//double r2 = p0 -> r2_mi ( p1 ) ;

		for ( unsigned int iv0 = 0 ; iv0<orientational_order_vectors.size() ; iv0++ ) {

			p0->get_vector_in_global_frame ( &vector0 , &orientational_order_vectors.at(iv0) ) ;

			for ( unsigned int iv1 = 0 ; iv1<orientational_order_vectors.size() ; iv1++ ){

				num_vector_pairs++;

				p1->get_vector_in_global_frame ( &vector1 , &orientational_order_vectors.at(iv1) ) ;

				// Note: it is assumed that that vector0 and vector1 are unit vectors;
				double cos_theta = 0.0 ;
				for (int i = 0 ; i < 3 ; i++) cos_theta+=vector0.at(i)*vector1.at(i) ;

				/*double theta;
				if(cos_theta>1.0 && cos_theta<1.00001 ){   // Take care of round-off error
					theta = 0.0 ;
				} else if ( cos_theta>-1.00001 && cos_theta<-1.0 ) {
					theta = 3.14159265358979 ;
				} else {
					theta = acos(cos_theta) ;
				}*/

				// Fourth legencre polynomial P_4(x=cos(theta)) = 1/8( 35 x^4 - 30 x^2 + 3 )
				double cos_theta2=cos_theta*cos_theta   ;
				double cos_theta4=cos_theta2*cos_theta2 ;
				sum += 4.375*cos_theta4 - 3.75*cos_theta2 + 0.375;
				//cout << cos_theta << " ";
				//cout << "cos_theta = " <<  cos_theta  << " theta = " << theta << endl ;
			}
		}
		//cout << endl ;
	}

	if( num_vector_pairs > 0 ) {
		// scale_factor = ((1+2*P_4(0))/3)^-1 = 1.71428571428571
		return sum/(double)num_vector_pairs*1.71428571428571 ;
	} else {
		return 0.0 ;
	}
}


/**
 * Return local orientational order parameter of particle.
 *
 *   (cos^2(2*theta)-7/15)/(1-7/15)
 *
 *  Notes:       7/15  = 0.466666666666667
 *        1/(1 - 7/15) = 1.875
 */
double Potential::get_orientational_order_qubaticJCP115 ( Particle *p0 ) {

	//cout << sqrt(r2_cutoff_order_max) << " " << sqrt(r2_cutoff_order_min) << endl;

	int num_vector_pairs = 0    ;
	double sum = 0.0            ;

	vector<double> vector0( 3 , 0.0 ) ;
	vector<double> vector1( 3 , 0.0 ) ;

	for(unsigned int i = 0 ; i < p0->particles_bounded.size() ; i++ ) {

		Particle *p1 = p0 -> particles_bounded.at(i) ;

		//cout << p0->index << " - " <<  p1->index << ": ";

		//double r2 = p0 -> r2_mi ( p1 ) ;

		for ( unsigned int iv0 = 0 ; iv0<orientational_order_vectors.size() ; iv0++ ) {

			p0->get_vector_in_global_frame ( &vector0 , &orientational_order_vectors.at(iv0) ) ;

			for ( unsigned int iv1 = 0 ; iv1<orientational_order_vectors.size() ; iv1++ ){

				num_vector_pairs++;

				p1->get_vector_in_global_frame ( &vector1 , &orientational_order_vectors.at(iv1) ) ;

				// Note: it is assumed that that vector0 and vector1 are unit vectors;
				double cos_theta = 0.0 ;
				for (int i = 0 ; i < 3 ; i++) cos_theta+=vector0.at(i)*vector1.at(i) ;

				double theta;
				if(cos_theta>1.0 && cos_theta<1.00001 ){   // Take care of round-off error
					theta = 0.0 ;
				} else if ( cos_theta>-1.00001 && cos_theta<-1.0 ) {
					theta = 3.14159265358979 ;
				} else {
					theta = acos(cos_theta) ;
				}

				double tmp = cos(2.0*theta) ;
				sum += (tmp*tmp)            ;
				//cout << cos_theta << " ";
				//cout << "cos_theta = " <<  cos_theta  << " theta = " << theta << endl ;
			}
		}
		//cout << endl ;
	}

	if( num_vector_pairs > 0 ) {
		return (sum/(double)num_vector_pairs-0.466666666666667)*1.875;
	} else {
		return 0.0 ;
	}
}


/**
 * Return unit vector and squard length of vector between partivle p0 and p1.
 */
void Potential::get_p0p1_unit_vector( double* p0p1_unit_vector , double *r2 , Particle *p0 , Particle *p1){
	p0p1_unit_vector[0] = p0->get_dx_img(p1,0);
	p0p1_unit_vector[1] = p0->get_dx_img(p1,1);
	p0p1_unit_vector[2] = p0->get_dx_img(p1,2);

	*r2 = 0.0 ;

	for ( int i = 0 ; i < 3 ; i++ ) *r2+= p0p1_unit_vector[i]*p0p1_unit_vector[i];

	double r = sqrt(*r2);

	for ( int i = 0 ; i < 3 ; i++ ) p0p1_unit_vector[i] /= r ;

}

//void Potential::reset_arm_allegiance_variables ( Particle * p , unsigned int ia  ) {
//}

/**
 * Reset the arm allegiance variables of Particle p
 */
void Potential::reset_arm_allegiance_variables ( Particle * p ) {

	p->particles_bounded.clear();

	// Loop arms
	for ( unsigned int ia = 0  ; ia < p -> arms.size() ; ia++ ) {

		p->arms.at(ia)->clear_allegiance() ;

		// Loop neighbor list particles
		for( unsigned int ipn = 0 ; ipn < p->pairs.size() ; ipn++ ) {

			unsigned int neighbour_index = p->pairs.at(ipn)->index ;

			double p0p1_unit_vector[] = {   p->get_dx_img(particles->at(neighbour_index),0) ,
											p->get_dx_img(particles->at(neighbour_index),1) ,
											p->get_dx_img(particles->at(neighbour_index),2) } ;
			double r2 = 0.0 ;
			for ( int i = 0 ; i < 3 ; i++ ) r2+= p0p1_unit_vector[i]*p0p1_unit_vector[i] ;
			double r = sqrt(r2) ;
			for ( int i = 0 ; i < 3 ; i++ ) p0p1_unit_vector[i] /= r ;

			// Loop arms of neighbor particle
			for ( unsigned int ian = 0  ; ian < particles->at(neighbour_index)->arms.size() ; ian++ ) {

				double arm_energy = get_arm_energy (    p0p1_unit_vector ,
														&r2 ,
														p->arms.at(ia) ,
														particles->at(neighbour_index)->arms.at(ian) ) ;
				if ( arm_energy < 0.0 ) {
					if ( only_bind_to_lowest_arm_energy == 1 || only_bind_to_lowest_arm_energy == 2 ) {
						p->arms.at(ia)->add_possible_allegiance( arm_energy , neighbour_index , ian , false ) ;
					} else {
						p->arms.at(ia)->add_possible_allegiance( arm_energy , neighbour_index , ian , true ) ;
					}
				}
			}
		}



		if ( only_bind_to_lowest_arm_energy == 0 ) {
			// Bind to all arms ... all is good.
		} else if ( only_bind_to_lowest_arm_energy == 2 ) {

			double lowest_energy = 0.0 ;
			int selected_arm = -1 ;

			for ( unsigned int i = 0 ; i < p->arms.at(ia)->possible_allegiances_energies.size() ; i++ ) {

				double next_energy = p->arms.at(ia)->possible_allegiances_energies.at(i) ;

				if ( next_energy < lowest_energy ) {

					lowest_energy = next_energy ;
					selected_arm = i;

				}

			}

			if( p->arms.at(ia)->possible_allegiances_energies.size() > 0 ){
				p->arms.at(ia)->possible_allegiances_is_bound.at(selected_arm) = true ;

				Particle * p_other = particles->at( p->arms.at(ia)->possible_allegiances_connecting_particle.at(selected_arm) );
				p->particles_bounded.push_back( p_other ) ;
			}

		} else if ( only_bind_to_lowest_arm_energy == 1 ) {
			// Loop possible allegiance and bind to the lowest energy one, and set ''is_bound'' variable

			double lowest_energy = 0.0 ;
			int selected_arm = -1 ;

			for ( unsigned int i = 0 ; i < p->arms.at(ia)->possible_allegiances_energies.size() ; i++ ) {

				double next_energy = p->arms.at(ia)->possible_allegiances_energies.at(i) ;

				if ( next_energy < lowest_energy ) {

					lowest_energy = next_energy ;
					selected_arm = i;

				}

			}

			// Test if this is also the lowest energy arm of the other particle.
			if ( selected_arm > -1 ) {

				Particle * p_other = particles->at( p->arms.at(ia)->possible_allegiances_connecting_particle.at(selected_arm) );
				Arm * a_other = p_other->arms.at( p->arms.at(ia)->possible_allegiances_connecting_arm.at( selected_arm ) ) ;

				double lowest_energy_neighbour = 0.0 ;
				//int selected_arm_neighbour = -1 ;

				for ( unsigned int p_other_pair = 0 ; p_other_pair < p_other->pairs.size() ; p_other_pair++ ) {
					for ( unsigned int yet_another_arm = 0 ; yet_another_arm < p_other->pairs.at(p_other_pair)->arms.size() ; yet_another_arm++ ){

						double r2 = 0.0;
						double p0p1_unit_vector[3] = { 0.0 , 0.0 , 0.0 };
						get_p0p1_unit_vector(p0p1_unit_vector , &r2 , p_other , p_other->pairs.at(p_other_pair) );

						double next_energy
						     = get_arm_energy    (
							   p0p1_unit_vector ,
							   &r2              ,
							   a_other          ,
							   p_other->pairs.at(p_other_pair)->arms.at(yet_another_arm)
						) ;

						if ( next_energy < lowest_energy_neighbour ) {
							lowest_energy_neighbour=next_energy ;
						}
					}
				}

				if( lowest_energy==lowest_energy_neighbour ) {
					p->arms.at(ia)->possible_allegiances_is_bound.at(selected_arm) = true ;
					p->particles_bounded.push_back( p_other ) ;
				}
			}
		} else {
			cout << "Unknown value of only_bind_to_lowest_arm_energy=" << only_bind_to_lowest_arm_energy << ". Exit 14541. " << endl ;
			exit(14541);
		}

	}





	/* OLD

	//reset_particle_indexes();

	//for( unsigned int ip = 0 ; ip < particles->size() ; ip++ ) {

		for ( unsigned int ia = 0  ; ia < particles->at(ip)-> arms.size() ; ia++ ) {

			particles->at(ip)->arms.at(ia)->clear_allegiance();

			for(  unsigned int ipn = 0 ; ipn < particles->at(ip)->pairs.size() ; ipn++ ) {

				unsigned int neighbour_index = particles->at(ip)->pairs.at(ipn)->index;

				double p0p1_unit_vector[] = {   particles->at(ip)->get_dx_img(particles->at(neighbour_index),0) ,
												particles->at(ip)->get_dx_img(particles->at(neighbour_index),1) ,
												particles->at(ip)->get_dx_img(particles->at(neighbour_index),2) } ;
				double r2 = 0.0 ;
				for ( int i = 0 ; i < 3 ; i++ ) r2+= p0p1_unit_vector[i]*p0p1_unit_vector[i];
				double r = sqrt(r2);
				for ( int i = 0 ; i < 3 ; i++ ) p0p1_unit_vector[i] /= r ;

				for ( unsigned int ian = 0  ; ian < particles->at(neighbour_index)->arms.size() ; ian++ ) {

					double arm_energy = get_arm_energy (    p0p1_unit_vector ,
															&r2 ,
															particles->at(ip)->arms.at(ia) ,
															particles->at(neighbour_index)->arms.at(ian) ) ;
					if ( arm_energy < 0.0 ) {

						particles->at(ip)->arms.at(ia)->add_possible_allegiance( arm_energy , neighbour_index , ian , true );

					}
				}
			}
		}
	//}
	 END OLD */
}

/**
 * Make sure that bindings are
 */
//void Potential::only_select_consistant_bonds ( Particle * p ) {

//}

