/*
 * Sim.cpp
 *
 *  Created on: Oct 6, 2011
 *      Author: urp
 */

#include "Sim.h"

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

#include "Potential.h"
#include "Particle.h"
#include "Arm.h"

Sim::Sim(string& variables_in) {

	variables = variables_in ;
	verbose   = get_variable ( variables , "verbose" , 9 ) ;

	seed      = get_variable ( variables , "seed" , 1 ) ;
	initialize_random_number ( ) ;

	double kT = get_variable ( variables , "kT" , 1.0 ) ;
	beta      = 1/kT ;

	chem_pot = get_variable ( variables , "chem_pot" , 0.0 ) ;

	time      = get_variable ( variables , "time" , 0.0 ) ;

	neighbour_list_updates = 0 ;

	// Set size of box
	bbox.clear ( ) ;
	bbox.push_back ( get_variable ( variables , "bboxX" , 6.0        ) ) ;
	bbox.push_back ( get_variable ( variables , "bboxY" , bbox.at(0) ) ) ;
	bbox.push_back ( get_variable ( variables , "bboxZ" , bbox.at(1) ) ) ;

	potential.set_parameters ( variables ) ;

	// Add particles
	//particle_counter        = 0 ;

	skin_distance           = get_variable ( variables , "skin_distance" , 0.3 ) ;
	int num_added_particles = get_variable ( variables , "num_added_particles" , 0 ) ;

	for ( int i = 0 ; i < num_added_particles ; i++ )  add_particle ( ) ;

	// Set-up potential;
	//potential.set_pair_eps(get_variable(variables,"pair_eps",1.0));
	//potential.set_pair_sigma(get_variable(variables,"pair_sigma",1.0));

	// Set-cutoff
	potential.particles = &particles;	// Assign particles to potential

	// For number of particles umbrella
	do_num_par_umbrella = false            ;
	num_par_center      = particles.size() ;
	num_par_kappa       = 0.1              ;

	// For changing volume changes
	do_volume_moves     = false ;
	type_of_volume_move = 0     ;
	delta_L_max         = 0.0   ;
	pressure            = 0.0   ;

}




Sim::~Sim() {
	// TODO Auto-generated destructor stub
}

/**
 * Use this function to (re)set kT. Ensure that the variable beta is also set.
void Sim::set_kT(double in_kT){

	kT   = in_kT ;
	beta = 1/kT  ;

}*/

/**
 * Update variables on the simulation. Should be run after particles and/or arms have been added or removed.
 */
void Sim::update(){

	update_particles_indexes();
	update_neighbour_list();
	update_arm_allegiance_variable();

}




/**
 * Set particles indexes
 */
void Sim::update_particles_indexes(){
	for( unsigned int p = 0 ; p < particles.size() ; p++ ) {
		particles.at(p)->index = p ;
	}
}




/**
 * Builds neighbor list
 */
void Sim::update_neighbour_list(){

	particles.at(0)->make_pair_list_all();

}




/**
 * Set arm-allegiance variables
 */
void Sim::update_arm_allegiance_variable(){

	for( unsigned int ip = 0 ; ip < particles.size() ; ip++ ) {
		potential.reset_arm_allegiance_variables( particles.at(ip) );
	}

}







/**
 * add particle with random coordinates
 */
void Sim::add_particle() {

	// Make random 4d unit vector for quaternion
	vector<double> q;
	double length_squared = 2.0;
	while (length_squared>1.0){
		q.clear();
		length_squared = 0.0;
		for( unsigned int i = 0 ; i < 4 ; i++ ){
			q.push_back(2.0*my_rand()-1.0);
			length_squared += q.at(i)*q.at(i);
		}
		//cout << length_squared << " ";
	}
	//cout << endl;
	double r = sqrt(length_squared);
	double q0 = q[0]/r ;
	double q1 = q[1]/r ;
	double q2 = q[2]/r ;
	double q3 = q[3]/r ;

	add_particle( 0         ,
				  my_rand() , my_rand() , my_rand() ,
				  0         , 0         , 0         ,
				  q0        , q1        , q2        ,        q3 );

	/*particles.push_back(
		new Particle(
			particle_counter++ , 0 ,
			&bbox,
			my_rand(),my_rand(),my_rand(),
			//tmp,0.0,0.0,
			0,0,0,
			//1.0,0.0,0.0,0.0,
			q0,q1,q2,q3,
			&particles,
			&potential.cut_distance,&skin_distance
		)
	);*/

}

/**
 * Add particle
 */
void Sim::add_particle(int i0,
			double x,double y,double z,
			int ii0,int ii1,int ii2,
			double dd0,double dd1,double dd2,double dd3) {

	unsigned int index = particles.size();

	particles.push_back(
		new Particle(
			index,i0,
			&bbox,
			x,y,z,
			ii0,ii1,ii2,
			dd0,dd1,dd2,dd3,
			&particles,
			&potential.cut_distance,&skin_distance,
			&neighbour_list_updates
		)
	);

}


void Sim::load_xyz(string& variables_in,string& filename){

	ifstream ifile (filename.c_str());

	if (ifile.is_open()) {
		string line;

		getline (ifile,line); // Number of atoms line
		int atoms_to_load = atoi (line.c_str());

		getline (ifile,line); // Header line

		cout << "Input file header: " << line << endl;
		cout << "atoms_to_load=" << atoms_to_load << endl;

		// Read time from header
		int read_time_from_input = get_variable ( variables , "read_time_from_input" , 0 ) ;
		if ( read_time_from_input!=0 ) {
			cout << "Note: Set time using input file header. Set read_time_from_input=0 to avoid this." << endl ;
			double in_time = get_header_variable ( line , "time" , time ) ;
			time = in_time ;
		} else {
			cout << "Note: Do not use time of input file header. Set read_time_from_input=1 to set time." << endl ;
		}

		// Read boundary box variables from header
		int read_bbox_from_input = get_variable ( variables , "read_bbox_from_input" , 0 ) ;
		if ( read_bbox_from_input!=0 ) {

			cout << "Note: Attempt to reset boundary box (bbox) to that of the input file header (Header example: bboxX=6.0 bboxY=6.0 bboxZ=6.0). Set read_bbox_from_input=0 to avoid this." << endl;

			double in_bboxX = get_header_variable ( line , "bboxX" , bbox.at(0) ) ;
			double in_bboxY = get_header_variable ( line , "bboxY" , bbox.at(1) ) ;
			double in_bboxZ = get_header_variable ( line , "bboxZ" , bbox.at(2) ) ;
			//cout << "Note: (Re)set bboxX=" << in_bboxX << endl ;
			bbox.at(0) = in_bboxX ;
			bbox.at(1) = in_bboxY ;
			bbox.at(2) = in_bboxZ ;

		} else {

			cout << "Note: Do not use boundary box of input file header. Set read_bbox_from_input=1 to avoid this." << endl ;

		}

		// Read temperature variable from header
		int read_kT_from_input = get_variable ( variables , "read_kT_from_input" , 0 ) ;
		if (read_kT_from_input!=0  ){

			cout << "Note: Set kT to that of input-file header. Set read_kT_from_input=0 to avoid this." << endl;
			double kT = get_header_variable ( line , "kT" , 1/beta ) ;
			beta = 1/kT;

		} else {

			cout << "Note: Do not use kT from input file header. Set read_kT_from_input=1 to set kT to that of input header." << endl ;

		}

		for (int i=0;i<atoms_to_load;i++){
			getline (ifile,line);

			char * pEnd;
			int    i0;
			double d0, d1, d2;
			int ii0,ii1,ii2;
			double dd0,dd1,dd2,dd3;

			i0  = strtol (line.c_str(),&pEnd,10);
			d0  = strtod (pEnd,&pEnd);
			d1  = strtod (pEnd,&pEnd);
			d2  = strtod (pEnd,&pEnd);
			ii0 = strtol (pEnd,&pEnd,10);
			ii1 = strtol (pEnd,&pEnd,10);
			ii2 = strtol (pEnd,&pEnd,10);
			dd0 = strtod (pEnd,&pEnd);
			dd1 = strtod (pEnd,&pEnd);
			dd2 = strtod (pEnd,&pEnd);
			dd3 = strtod (pEnd,&pEnd);

			double x_reduced = d0/bbox.at(0);
			double y_reduced = d1/bbox.at(1);
			double z_reduced = d2/bbox.at(2);

			/*vector<Arm> arms;
			arms.push_back ( Arm(  1.0 , 0.0 , 0.0 ) ) ;
			arms.push_back ( Arm( -1.0 , 0.0 , 0.0 ) ) ;
			arms.push_back ( Arm(  0.0 , 1.0 , 0.0 ) ) ;
			arms.push_back ( Arm(  0.0 ,-1.0 , 0.0 ) ) ;
			arms.push_back ( Arm(  0.0 , 0.0 , 1.0 ) ) ;
			arms.push_back ( Arm(  0.0 , 0.0 ,-1.0 ) ) ;*/

			add_particle( 0         ,
					      x_reduced , y_reduced , z_reduced ,
					      ii0       , ii1       , ii2       ,
					      dd0       , dd1       , dd2       ,        dd3 );

			/*particles.push_back(new Particle(
					particle_counter++,i0,
					&bbox,
					x_reduced,y_reduced,z_reduced,
					ii0,ii1,ii2,
					dd0,dd1,dd2,dd3,
					&particles,
					&potential.cut_distance,&skin_distance
					));
				*/
		}
	}else{
		cout << "Error: Unable to open file with input coordinates (" << filename << "). Exit." << endl;
		exit (27636);
	}
}



void Sim::print(){

	// Print box information
	cout << endl << "   ..:: Simulation variables ::.." << endl ;
	cout << "time = " << time << endl ;
	cout << "beta = " << beta << endl;
	cout << "boundary_box = [ ";
	for(int i = 0 ; i < 2 ; i++ ){
		cout << bbox.at(i) << " , ";
	}
	cout << bbox.at(2) << " ]" << endl;
	cout << "seed = " << seed << endl;
	cout << "skin_distance = " << skin_distance << endl;
	cout << "kT: " << 1/beta << endl;
	cout << "Potential energy: " << potential.get_internal_energy() << endl;
	cout << "Volume: " << bbox.at(0)*bbox.at(1)*bbox.at(2) << "\n";
	cout << "Number of particles (sim): " << particles.size() << "\n";
	cout << "Number of particles (potential): " << potential.particles->size() << "\n";
	cout << "Density: " << (double)particles.size()/(bbox.at(0)*bbox.at(1)*bbox.at(2)) << "\n";
	cout << "Average length of pair list: " << get_avg_pair_list_length() << endl;
	cout << endl;

	// Print particles
	/*for (unsigned int i = 0 ; i < particles.size() ; i++ ){
		particles.at(i)->print();
		//cout << "particle_energy = " << potential.particle_energy(particles.at(i)) << endl;
	}*/

}






/**
 * Print energy to standard output.
 */
void Sim::print_energy(){
	// Print box information
	cout << "energy: " << potential.get_internal_energy() << endl;
}




/**
 * Return string with energy (and other system variables)
 */
string Sim::get_ener_string(){
	stringstream out;

	out << time                            << " " ;      // 1
	out << potential.get_internal_energy() << " " ;    // 2
	out << particles.size()                << " " ;  // 3
	out << get_volume()				 	   << " " ;  // 4
	out << bbox.at(0)                      << " " ;  // 5
	out << bbox.at(1)                      << " " ;  // 6
	out << bbox.at(2)                      << " " ;  // 7
	out << potential.get_orientational_order()        << " " ;  // 8
	out << potential.get_orientational_order_energy() << " " ;  // 9
	out << potential.get_pinning_energy()  << " " ;  // 10
	out << potential.get_total_energy()    << " " ;  // 11
	return out.str ( ) ;
}





/**
 * Return string with coordinates in the xyz format
 */
string Sim::get_coordinates_xyz(){

	stringstream out;

	// Make header
	out << particles.size() << "\n";
	out << "MC simulation: time=" << time
		<< " bboxX=" << bbox.at(0)
		<< " bboxY=" << bbox.at(1)
		<< " bboxZ=" << bbox.at(2)
		<< " beta=" << beta
		<< " kT=" << 1/beta
		<< " V=" << bbox.at(0)*bbox.at(1)*bbox.at(2)
		<< " U=" << potential.get_internal_energy()
		<< " Upin=" << potential.get_pinning_energy()
		<< " Q=" << potential.get_orientational_order()
		<< " columns={type,x,y,z,ix,iy,iz,q0,q1,q2,q3,u,upin,Qi}"
	;

	// Write coordinates
	for(unsigned int i = 0 ; i < particles.size() ; i++ ){
		out<<"\n" << particles.at(i)->type
				  << " " << particles.at(i)->get_min_img_x()
				  << " " << particles.at(i)->get_min_img_y()
				  << " " << particles.at(i)->get_min_img_z()
				  << " " << particles.at(i)->img.at(0)
				  << " " << particles.at(i)->img.at(1)
				  << " " << particles.at(i)->img.at(2)
				  << " " << particles.at(i)->q.at(0)
				  << " " << particles.at(i)->q.at(1)
				  << " " << particles.at(i)->q.at(2)
				  << " " << particles.at(i)->q.at(3)
				  << " " << potential.get_particle_energy (    particles.at(i) )
				  << " " << potential.get_pinning_energy  (    particles.at(i) )
				  << " " << potential.get_orientational_order( particles.at(i) );
		;
	}
	return out.str();
}

/**
 * Return string with allegiance network
 */
string Sim::get_coordinates_network(){

	stringstream out;

	// Count number of contacts
	int num_bonds = 0;
	for( unsigned int p = 0 ; p < potential.particles->size() ; p++ ) {
		for ( unsigned int a = 0  ; a < potential.particles->at(p)-> arms.size() ; a++ ) {
			num_bonds += potential.particles->at(p)->arms.at(a)->possible_allegiances_connecting_particle.size();
		}
	}

	// Write header
	out << num_bonds << " " << potential.particles->size() << " " << time;
	out << endl << "# Above: [ number of connections ; number of nodes ; time ]. Below: [ particle ; arm ; energy ; connecting particle ; connecting arm ; is connected ]";

	// Loop arms and print bonds
	for( unsigned int p = 0 ; p < potential.particles->size() ; p++ ) {
		for ( unsigned int a = 0  ; a < potential.particles->at(p)-> arms.size() ; a++ ) {
			for ( unsigned int i = 0 ; i < potential.particles->at(p)->arms.at(a)->possible_allegiances_connecting_particle.size() ; i++) {
				out << endl << p
					<< " " << a
					<< " " << potential.particles->at(p)->arms.at(a)->possible_allegiances_energies.at(i)
					<< " " << potential.particles->at(p)->arms.at(a)->possible_allegiances_connecting_particle.at(i)
					<< " " << potential.particles->at(p)->arms.at(a)->possible_allegiances_connecting_arm.at(i)
				;
				if ( potential.particles->at(p)->arms.at(a)->possible_allegiances_is_bound.at(i) ) {
					out << " 1" ;
				} else {
					out << " 0" ;
				}
			}
		}
	}

	return out.str();

}


string Sim::get_coordinates_vmd ( int type_of_vmd_color ) {

	double max_energy =  0.0;
	double min_energy = -1.0*(double)particles.at(0)->arms.size() * potential.get_arm_eps();

	//cout << "min_energy= " <<  min_energy << endl;

	return get_coordinates_vmd(type_of_vmd_color, min_energy, max_energy);
}


string Sim::get_coordinates_vmd( int type_of_vmd_color , double min_energy , double max_energy ) {

	double max_energy_arm = 0.0;
	double min_energy_arm = - potential.get_arm_eps ( ) ;

	return get_coordinates_vmd (  type_of_vmd_color,
								  min_energy     ,
			                      max_energy     ,
			                      min_energy_arm ,
			                      max_energy_arm );
}


/**
 *
 * Usage:
 *
 * 	FILE *    ofile_vmd_end = fopen ( "end.vmd" , "w" ) ;
 *  fprintf ( ofile_vmd_end , "%s", sim.get_coordinates_vmd().c_str() ) ;
 *	fclose  ( ofile_vmd_end ) ;
 */
//string Sim::get_coordinates_vmd(double min_energy,double max_energy) {
string Sim::get_coordinates_vmd ( int type_of_vmd_color,
							      double min_energy    ,
							      double max_energy    ,
							      double min_energy_arm,
							      double max_energy_arm  ) {

	stringstream out;

	int resolution = 8 ;

	double arm_radius = 0.1 ;
	double arm_length = 0.5*potential.get_arm_length();

	for(unsigned int particle = 0 ; particle < particles.size() ; particle++ ){

		  // Set particle color

		//double max_energy =  0.0;
		//double min_energy = -6.0;

		double color = 0.0;
		if( type_of_vmd_color == 1 ) {
			color = potential.get_orientational_order( particles.at(particle) ) ;
			color *= 1023.;
			color += 33;
			if( color<33 ) color=33;
			if( color>1023 ) color=1023;
		} else { // The default is using energy as coloring (type_of_vmd_color = 0)
			color = potential.get_particle_energy ( particles.at(particle) );
			color -=	min_energy;
			color /= ( max_energy - min_energy );
			color *= 1023.;
			color += 33;
			if( color<33 ) color=33;
			if( color>1023 ) color=1023;
		}

	      // Print particle centers

		out << "graphics 0 color " << (int)floor(color) << endl;

		double particle_x = particles.at(particle)->get_min_img_x();
		double particle_y = particles.at(particle)->get_min_img_y();
		double particle_z = particles.at(particle)->get_min_img_z();

		out << "graphics 0 sphere {"  << particle_x << " "
									  << particle_y << " "
				                      << particle_z << "} "
			 << "radius " << 0.5*potential.get_particle_size() << " "
			 << "resolution " << resolution
			 << endl ;

		  // Print arms
		for ( unsigned int arm = 0 ; arm < particles.at(particle)->arms.size() ; arm++ ) {

			if( type_of_vmd_color == 1 ) {
				// Same as center sphere color = do nothing
			}else{
			   // Color arm according to arm energy
			   color = particles.at(particle)->arms.at(arm)->get_energy_from_possible_allegiance();
			   color -=	min_energy_arm;
			   color /= ( max_energy_arm - min_energy_arm );
			   color *= 1023.;
			   color += 33;
			   if( color<33 ) color=33;
			   if( color>1023 ) color=1023;
			}


		    out << "graphics 0 color " << (int)floor(color) << endl;

			double arm_x = particle_x + arm_length*(particles.at(particle)->arms.at(arm)->u_global.at(0)) ;
			double arm_y = particle_y + arm_length*(particles.at(particle)->arms.at(arm)->u_global.at(1)) ;
			double arm_z = particle_z + arm_length*(particles.at(particle)->arms.at(arm)->u_global.at(2)) ;

			out << "graphics 0 cylinder { " << particles.at(particle)->get_min_img_x() << " "
											<< particles.at(particle)->get_min_img_y() << " "
											<< particles.at(particle)->get_min_img_z() << " "
								            << "} { "   << arm_x
								            <<  " "     << arm_y
								            <<  " "     << arm_z
								            << " } radius " << arm_radius
								            << " resolution " << resolution
			                                << endl ;
			//if(particles.at(particle)->arms.at(arm)->arm_have_a_bond_in_possible_allegiance())
			//	out << " graphics 0 color gray " << endl;


			out << "graphics 0 sphere {"  << arm_x << " "
										  << arm_y << " "
					                      << arm_z << "} "
									      << "radius " << arm_radius << " "
				                          << "resolution " << resolution
				                          << endl ;

			//if(particles.at(particle)->arms.at(arm)->arm_have_a_bond_in_possible_allegiance())
			//	out << "graphics 0 color " << (int)floor(color) << endl;

		}

	}

	double bbox_radius = 0.2;

	// Print color bar.
	int len = 216;
	for ( int i = 0 ; i < len ; i++ ) {
		double frac = (double)i/(double)len ;
		out << "graphics 0 color " << (int)floor( frac*(double)1023 ) + 33 << endl;
		out << "graphics 0 sphere { " << frac*bbox.at(0) << " " << -bbox_radius << " " << -bbox_radius << " } "
				<< "radius " << bbox_radius << " "
				<< "resolution " << resolution
				<< endl ;
	}

	// Print periodic box
	out << "graphics 0 color blue" << endl;


	for ( double         x = 0.0 ; x <= 1.0 ; x += 1.0 ) {
		for ( double     y = 0.0 ; y <= 1.0 ; y += 1.0 ) {
			for ( double z = 0.0 ; z <= 1.0 ; z += 1.0 ) {
				out << "graphics 0 sphere { " << x*bbox.at(0) << " " << z*bbox.at(1) << " " << y*bbox.at(2) << " } "
					<< "radius " << bbox_radius << " "
					<< "resolution " << resolution
					<< endl ;
			}
		}
	}

	out << "graphics 0 cylinder "
			<<" { " << 0.0        << " " << 0.0        << " " << 0.0        << " } "
			<<" { " << bbox.at(0) << " " << 0.0        << " " << 0.0        << " } "
			<< "radius " << bbox_radius << " " << "resolution " << resolution << endl ;

	out << "graphics 0 cylinder "
			<<" { " << 0.0        << " " << 0.0        << " " << 0.0        << " } "
			<<" { " << 0.0        << " " << bbox.at(1) << " " << 0.0        << " } "
			<< "radius " << bbox_radius << " " << "resolution " << resolution << endl ;

	out << "graphics 0 cylinder "
			<<" { " << bbox.at(0) << " " << bbox.at(1) << " " << 0.0        << " } "
			<<" { " << bbox.at(0) << " " << 0.0        << " " << 0.0        << " } "
			<< "radius " << bbox_radius << " " << "resolution " << resolution << endl ;

	out << "graphics 0 cylinder "
			<<" { " << bbox.at(0) << " " << bbox.at(1) << " " << 0.0        << " } "
			<<" { " << 0.0        << " " << bbox.at(1) << " " << 0.0        << " } "
			<< "radius " << bbox_radius << " " << "resolution " << resolution << endl ;

	out << "graphics 0 cylinder "
			<<" { " << 0.0        << " " << 0.0        << " " << bbox.at(2) << " } "
			<<" { " << bbox.at(0) << " " << 0.0        << " " << bbox.at(2) << " } "
			<< "radius " << bbox_radius << " " << "resolution " << resolution << endl ;

	out << "graphics 0 cylinder "
			<<" { " << 0.0        << " " << 0.0        << " " << bbox.at(2) << " } "
			<<" { " << 0.0        << " " << bbox.at(1) << " " << bbox.at(2) << " } "
			<< "radius " << bbox_radius << " " << "resolution " << resolution << endl ;

	out << "graphics 0 cylinder "
			<<" { " << bbox.at(0) << " " << bbox.at(1) << " " << bbox.at(2) << " } "
			<<" { " << bbox.at(0) << " " << 0.0        << " " << bbox.at(2) << " } "
			<< "radius " << bbox_radius << " " << "resolution " << resolution << endl ;

	out << "graphics 0 cylinder "
			<<" { " << bbox.at(0) << " " << bbox.at(1) << " " << bbox.at(2) << " } "
			<<" { " << 0.0        << " " << bbox.at(1) << " " << bbox.at(2) << " } "
			<< "radius " << bbox_radius << " " << "resolution " << resolution << endl ;

	out << "graphics 0 cylinder "
			<<" { " << 0.0        << " " << 0.0        << " " << 0.0        << " } "
			<<" { " << 0.0        << " " << 0.0        << " " << bbox.at(2) << " } "
			<< "radius " << bbox_radius << " " << "resolution " << resolution << endl ;

	out << "graphics 0 cylinder "
			<<" { " << bbox.at(0) << " " << 0.0        << " " << 0.0        << " } "
			<<" { " << bbox.at(0) << " " << 0.0        << " " << bbox.at(2) << " } "
			<< "radius " << bbox_radius << " " << "resolution " << resolution << endl ;

	out << "graphics 0 cylinder "
			<<" { " << 0.0        << " " << bbox.at(1)  << " " << 0.0        << " } "
			<<" { " << 0.0        << " " << bbox.at(1)  << " " << bbox.at(2) << " } "
			<< "radius " << bbox_radius << " " << "resolution " << resolution << endl ;

	out << "graphics 0 cylinder "
			<<" { " << bbox.at(0) << " " << bbox.at(1)  << " " << 0.0        << " } "
			<<" { " << bbox.at(0) << " " << bbox.at(1)  << " " << bbox.at(2) << " } "
			<< "radius " << bbox_radius << " " << "resolution " << resolution << endl ;




	return out.str();

}

/**
 * Return volume of simulation box
 */
double Sim::get_volume(){
	return bbox.at(0)*bbox.at(1)*bbox.at(2) ;
}



double Sim::get_avg_pair_list_length(){

	double sum=0.0;

	for ( unsigned int i = 0 ; i < particles.size() ; i++ ) {

		sum += particles.at(i)->pairs.size();

	}

	return sum/(double)particles.size();

}


/**
 * Attempt a translation move.
 *
 * Return 1 if accepted, and 0 otherwise.
 *
 * Input:
 *   maximum step size
 */
int Sim::attempt_a_move_translation ( double *max_step_size ) {

	int out = 0;

	int p = (int)floor(my_rand()*particles.size());
	int dim = (int)floor(my_rand()*3);
	double step = 2.0*(my_rand()-0.5);

	// Test of selection of particle and dimension is valid (some random generators may give 1.0000 ... ).
	if(p<0)                  cout << "Warning in Sim::attempt_a_move_translation(double *max_step_size): invalid particle index, p=" << p << ". particles.size()=" << particles.size() << endl;
	if(p>particles.size()-1) cout << "Warning in Sim::attempt_a_move_translation(double *max_step_size): invalid particle index, p=" << p << ". particles.size()=" << particles.size() << endl;
	if(dim<0)                cout << "Warning in Sim::attempt_a_move_translation(double *max_step_size): dim=" << dim << endl;
	if(dim>3)                cout << "Warning in Sim::attempt_a_move_translation(double *max_step_size): dim=" << dim << endl;

	step *= *max_step_size ;

	// Move particle

	double current_energy = potential.get_particle_energy (  particles.at(p) )
			              + potential.get_orientational_order_energy ( )
						  + potential.get_pinning_energy  (  particles.at(p) ) ;

	particles.at(p)->translate_step(dim,step);

	double new_energy     = potential.get_particle_energy ( particles.at(p) )
			              + potential.get_orientational_order_energy ( )
			              + potential.get_pinning_energy  ( particles.at(p) ) ;

	double DeltaE = new_energy - current_energy ;

	//cout << step << endl;

	if ( my_rand() < exp( -1.0 * beta * DeltaE ) ) {

		// accept move
		out=1;

	} else {

		// reject move, move back
		out=0;
		particles.at(p)->translate_step(dim,-step);

	}

	return out;
}



/**
 *
 * Attempt a rotation move.
 *
 *  TODO when particles are moved, the arms_variable should be updated
 */
int Sim::attempt_a_move_rotation(double *max_step_size){

	int p = (int)floor(my_rand()*particles.size());

	// Test of selection of particle is valid (some random generators may give 1.0000... ).
	if(p<0) cout << "Warning in Sim::attempt_a_move_rotation(double *max_step_size): invalid particle index, p=" << p << ". particles.size()=" << particles.size() << endl;
	if(p>particles.size()-1) cout << "Warning in Sim::attempt_a_move_rotation(double *max_step_size):: invalid particle index, p=" << p << ". particles.size()=" << particles.size() << endl;

	//cout << endl << endl;
	//for(int i = 0 ; i < 4 ; i ++ ) cout << particles.at(p)->q.at(i) << " ";

	// Make random vector inside 4D sphere.
	vector<double> step;
	double length_squared  = 2.0;
	while ( length_squared > 1.0 ) {
		step.clear ( ) ;
		length_squared = 0.0 ;
		for( unsigned int i = 0 ; i < 4 ; i++ ){
			step.push_back(2.0*my_rand()-1.0);
			length_squared += step.at(i)*step.at(i);
		}
	}

	double length = sqrt(length_squared);
	double final_length = my_rand();
	// Put on unit circle and set length of vector
	for( int i = 0 ; i < 4 ; i++ ) {
		step.at(i) /= length;
		step.at(i) *= *max_step_size;
		step.at(i) *= final_length;
	}

	// Store old quaternian
	vector<double> q_old;
	q_old.clear();
	for(int i = 0 ; i < 4 ; i++ ) {
		q_old.push_back(particles.at(p)->q.at(i));
	}

	// Rotate
	double current_energy = potential.get_particle_energy(particles.at(p))
			              + potential.get_orientational_order_energy() ;
	//double current_energy = potential.get_total_energy();
	particles.at(p)->rotate_step(step);
	double new_energy = potential.get_particle_energy(particles.at(p))
			          + potential.get_orientational_order_energy() ;
	//double new_energy = potential.get_total_energy();

	double DeltaE=new_energy-current_energy;

	// Metropolis criteria
	int out;


	//cout << endl;
	//for(int i = 0 ; i < 4 ; i ++ ) cout << particles.at(p)->q.at(i) << " ";
	//cout << endl << exp(-beta*DeltaE) ;

	if(my_rand()<exp(-beta*DeltaE)){
		// accept move
		out=1;
	}else{
		// reject move, rotate back
		out=0;

		particles.at(p)
				-> q.at(0)
				= q_old.at(0);

		particles.at(p)
				-> q.at(1)
				= q_old.at(1);

		particles.at(p)
				-> q.at(2)
				= q_old.at(2);

		particles.at(p)
				-> q.at(3)
				= q_old.at(3);

		/*for( int i = 0 ; i < 4 ; i++ ) {
			particles.at(p)
					-> q.at(i)
					= q_old.at(i);
		}*/

	}

	//cout << endl;
	//for(int i = 0 ; i < 4 ; i ++ ) cout << particles.at(p)->q.at(i) << " ";

	return out;
}

/**
 * Attempt a volume change. To produce NpT ensemble.
 * Rebuild simulation
 */
int Sim::attempt_a_move_volume ( ) {

	// Make an isotropic volume move. Note is is assumed that bbox.at(0)=bbox.at(1)=bbox.at(2)
	if ( do_volume_moves && type_of_volume_move==1 ) {
		double V_current =  bbox.at(0)*bbox.at(1)*bbox.at(2) ;

		double current_energy = potential.get_internal_energy ( )
				              + potential.get_orientational_order_energy()
							  + potential.get_pinning_energy();

		double delta_L = ( 2.0 * my_rand() - 1.0) * delta_L_max ;
		double bbox_old = bbox.at(0) ;
		bbox.at(0)    += delta_L    ;
		bbox.at(1)     = bbox.at(0) ;
		bbox.at(2)     = bbox.at(0) ;

		if ( bbox.at(0)<potential.cut_distance*2.0 ) {
			cout << "Error in Sim::attempt_a_move_volume: Volume move (barostat) with yielded to short box length. bbox.at(0) = " << bbox.at(0) << " delta_L = " << delta_L << " potential.cut_distance*2.0=" << potential.cut_distance*2.0 << " . Exit." << endl;
			exit(11701) ;
		}

		update();
		double V_new      = bbox.at(0)*bbox.at(1)*bbox.at(2) ;
		double deltaV     = V_new - V_current                ;
		double new_energy = potential.get_internal_energy ( )
				          + potential.get_orientational_order_energy()
						  + potential.get_pinning_energy()   ;
		double deltaU     = new_energy - current_energy      ;
		double lnVV       = log( V_new/V_current )           ;

		if( my_rand() < exp ( - beta*deltaU - beta*pressure*deltaV + ( (double)particles.size( ) + 1.0 )*lnVV ) ) {

			//cout << "ACCEPTED" << endl ;

			return 1 ;

		} else {

			bbox.at(0) = bbox_old ;
			bbox.at(1) = bbox_old ;
			bbox.at(2) = bbox_old ;
			update ( ) ;

			//cout << "REJECTED" << endl ;

			return 0 ;
		}

	} else if ( do_volume_moves && type_of_volume_move==2 ) {	// Anisotropic move

		double V_current      = bbox.at(0)*bbox.at(1)*bbox.at(2) ;
		double current_energy = potential.get_internal_energy ( )
				              + potential.get_orientational_order_energy ( )
							  + potential.get_pinning_energy ( ) ;

		vector<double> bbox_current;
		for(unsigned int i = 0 ; i < 3 ; i++ ){

			bbox_current.push_back(bbox.at(i));

			bbox.at(i) += ( 2.0 * my_rand() - 1.0) * delta_L_max ;

			if ( bbox.at(i)<potential.cut_distance*2.0 ) {
				cout << "Error in Sim::attempt_a_move_volume: Volume move (barostat) with yielded to short box length. bbox.at("<< i <<") = " << bbox.at(i) << " potential.cut_distance*2.0=" << potential.cut_distance*2.0 << " . Exit." << endl;
				exit(11702) ;
			}
		}

		update();
		double V_new      = bbox.at(0)*bbox.at(1)*bbox.at(2) ;
		double deltaV     = V_new - V_current                ;
		double new_energy = potential.get_internal_energy ( )
				          + potential.get_orientational_order_energy ( )
						  + potential.get_pinning_energy ( ) ;
		double deltaU     = new_energy - current_energy      ;
		double lnVV       = log( V_new/V_current )           ;

		if( my_rand() < exp ( - beta*deltaU - beta*pressure*deltaV + ( (double)particles.size( ) + 1.0 )*lnVV ) ) {

			//cout << "ACCEPTED" << endl ;
			return 1 ;

		} else {

			bbox.at(0) = bbox_current.at(0) ;
			bbox.at(1) = bbox_current.at(1) ;
			bbox.at(2) = bbox_current.at(2) ;
			update ( ) ;

			//cout << "REJECTED" << endl ;
			return 0 ;
		}

	} else if ( do_volume_moves && type_of_volume_move==3 ) {	// Anisotropic move in tetragonal geometry (a*a*c)

		double V_current      = bbox.at(0)*bbox.at(1)*bbox.at(2) ;
		double current_energy = potential.get_internal_energy ( )
				              + potential.get_orientational_order_energy ( )
							  + potential.get_pinning_energy ( ) ;


		vector<double> bbox_current;
		for(unsigned int i = 0 ; i < 3 ; i++ ){

			bbox_current.push_back(bbox.at(i));

			bbox.at(i) += ( 2.0 * my_rand() - 1.0) * delta_L_max ;

			if ( bbox.at(i)<potential.cut_distance*2.0 ) {
				cout << "Error in Sim::attempt_a_move_volume: Volume move (barostat) with yielded to short box length. bbox.at("<< i <<") = " << bbox.at(i) << " potential.cut_distance*2.0=" << potential.cut_distance*2.0 << " . Exit." << endl;
				exit(11702) ;
			}
		}

		// Make the box in the tetragonal geometry
		bbox.at(1)=bbox.at(0);

		update();
		double V_new      = bbox.at(0)*bbox.at(1)*bbox.at(2) ;
		double deltaV     = V_new - V_current                ;
		double new_energy = potential.get_internal_energy ( )
				          + potential.get_orientational_order_energy()
						  + potential.get_pinning_energy ( ) ;
		double deltaU     = new_energy - current_energy      ;
		double lnVV       = log( V_new/V_current )           ;

		if( my_rand() < exp ( - beta*deltaU - beta*pressure*deltaV + ( (double)particles.size( ) + 1.0 )*lnVV ) ) {

			//cout << "ACCEPTED" << endl ;
			return 1 ;

		} else {

			bbox.at(0) = bbox_current.at(0) ;
			bbox.at(1) = bbox_current.at(1) ;
			bbox.at(2) = bbox_current.at(2) ;
			update ( ) ;

			//cout << "REJECTED" << endl ;
			return 0 ;
		}
	} else if ( do_volume_moves && type_of_volume_move==4 ) {	// Anisotropic move in tetragonal geometry (a*a*c)

			double V_current      = bbox.at(0)*bbox.at(1)*bbox.at(2) ;
			double current_energy = potential.get_internal_energy ( ) +
					                potential.get_orientational_order_energy()
								  + potential.get_pinning_energy ( ) ;


			/*vector<double> bbox_current;
			for(unsigned int i = 0 ; i < 3 ; i++ ){

				bbox_current.push_back(bbox.at(i));

				bbox.at(i) += ( 2.0 * my_rand() - 1.0) * delta_L_max ;

				if ( bbox.at(i)<potential.cut_distance*2.0 ) {
					cout << "Error in Sim::attempt_a_move_volume: Volume move (barostat) with yielded to short box length. bbox.at("<< i <<") = " << bbox.at(i) << " potential.cut_distance*2.0=" << potential.cut_distance*2.0 << " . Exit." << endl;
					exit(11702) ;
				}
			}

			// Make the box in the tetragonal geometry
			bbox.at(1)=bbox.at(0);
		    */
			double bboxZ_current = bbox.at(2) ;
			bbox.at(2) += ( 2.0 * my_rand() - 1.0) * delta_L_max ;

			update();
			double V_new      = bbox.at(0)*bbox.at(1)*bbox.at(2) ;
			double deltaV     = V_new - V_current                ;
			double new_energy = potential.get_internal_energy ( )
					          + potential.get_orientational_order_energy() ;
							  + potential.get_pinning_energy ( ) ;
			double deltaU     = new_energy - current_energy      ;
			double lnVV       = log( V_new/V_current )           ;

			if( my_rand() < exp ( - beta*deltaU - beta*pressure*deltaV + ( (double)particles.size( ) + 1.0 )*lnVV ) ) {

				//cout << "ACCEPTED" << endl ;
				return 1 ;

			} else {

				bbox.at(2) = bboxZ_current ;
				update ( ) ;

				//cout << "REJECTED" << endl ;
				return 0 ;
			}
	} else {

		return 0 ;

	}

}

void Sim::enable_volume_changes ( int in_type_of_volume_move , double in_delta_L_max , double in_pressure ) {

	cout << "Note:  Set type_of_volume_move, delta_L_max, pressure to enable NpT simulations:"      << endl ;
	cout << "Note:      type_of_volume_move = 0   ->   No volume moves"                             << endl ;
	cout << "Note:      type_of_volume_move = 1   ->   Isotropic volume moves                (cubic box       ; a*a*a ) "        << endl ;
	cout << "Note:      type_of_volume_move = 2   ->   Anisotropic-orthorhombic volume moves (orthorhombic box; a*b*c ) " << endl ;
	cout << "Note:      type_of_volume_move = 3   ->   Anisotropic-tetragonal volume moves   (tetragonal box  ; a*a*c ) " << endl ;
	cout << "Note:      type_of_volume_move = 4   ->   Z-only moves                          (orthorhombic box; a*b*c ) " << endl ;

	if ( in_type_of_volume_move == 0 ) {

		do_volume_moves     = false ;
		type_of_volume_move = 0     ;
		delta_L_max         = 0.0   ;
		pressure            = 0.0   ;

		cout << "Note: No barostat ( type_of_volume_move = " << in_type_of_volume_move << " )" << endl ;

	} else if(in_type_of_volume_move == 1 ) {

		do_volume_moves     = true                       ;
		type_of_volume_move = in_type_of_volume_move     ;
		delta_L_max         = in_delta_L_max             ;
		pressure            = in_pressure                ;

		cout << "Note: Isotropic barostat is enabled ( type_of_volume_move = " << in_type_of_volume_move << " ) assumes a cubic box (bboxX=bboxY=bboxZ). bboxX is used." << endl;

		bbox.at(1) = bbox.at(0) ;
		bbox.at(2) = bbox.at(0) ;

	} else if( in_type_of_volume_move == 2 ) {
		do_volume_moves     = true                       ;
		type_of_volume_move = in_type_of_volume_move     ;
		delta_L_max         = in_delta_L_max             ;
		pressure            = in_pressure                ;

		cout << "Note: Anisotropic-orthorhombic barostat is enabled ( type_of_volume_move = " << in_type_of_volume_move << " )" << endl ;

	} else if( in_type_of_volume_move == 3 ) {

		do_volume_moves     = true                       ;
		type_of_volume_move = in_type_of_volume_move     ;
		delta_L_max         = in_delta_L_max             ;
		pressure            = in_pressure                ;

		cout << "Note: Anisotropic-tetragonal barostat is enabled ( type_of_volume_move = " << in_type_of_volume_move << " )."
				" The box is set to tetragonal geometry (bbox.at(1) = bbox.at(0))"<< endl ;
		bbox.at(1) = bbox.at(0) ;

	} else if( in_type_of_volume_move == 4 ) {

		do_volume_moves     = true                       ;
		type_of_volume_move = in_type_of_volume_move     ;
		delta_L_max         = in_delta_L_max             ;
		pressure            = in_pressure                ;

		cout << "Note: Z-only barostat is enabled ( type_of_volume_move = " << in_type_of_volume_move << " )." << endl;

	} else {

		cout << "Error: Invalid type of volume move (type_of_volume_move = " << in_type_of_volume_move << "). Exit." << endl ;
		exit(5837);

	}
}


int Sim::attempt_a_move_num_particles ( ) {

	//cout << "end "<< particles.back()->index << endl;

	if(do_num_par_umbrella && potential.do_orientational_order_umbrella){ // TODO Impliment the possiblilty of doing particles removals/insertion togeather with orientational umbrella
		cout << "Error: Performing particles removals/insertion is not implimentet with orientational umbrella (this is a TODO). Exit." << endl;
		exit(21491);
	}

	if( do_num_par_umbrella ) {

		int delta_num_par = 2*(int)floor(my_rand()*2) - 1 ;

		//cout << " delta_num_par = " << delta_num_par << endl;

		if( delta_num_par == -1 ) {	// Remove a particle

			int p = (int)floor(my_rand()*particles.size());

			double dE = -1.0 * potential.get_particle_energy(particles.at(p));

			//cout << "N=" << particles.size() << endl;
			//cout << " u(removal) " << potential.particle_energy(particles.at(p)) << endl << endl;

			// Note: ( N + dN - N0 )^2 - ( N - N0 )^2
			//         = dN^2 + 2*dN*( N - N0 )
			//         = 1 - 2*( N - N0 )
			dE += 0.5*num_par_kappa * ( 1.0 - 2.0*( (double)particles.size() - (double)num_par_center ) );


			// prefactor = N/V
			double prefactor = (double)particles.size();
                   prefactor /= (double)bbox.at(0) ;
                   prefactor /= (double)bbox.at(1) ;
                   prefactor /= (double)bbox.at(2) ;

			if( my_rand() < prefactor*exp ( - beta * ( dE + chem_pot ) ) ) {

				if( particles.size() > 1 ) {
					particles.erase ( particles.begin() + p ) ;
					update();
				}

				return 1;

			} else {

				return 0;

			}

		} else {  // Add a particle

			// prefactor = V/(N+1)
			double prefactor  =  (double)bbox.at(0) ;
					prefactor *=  (double)bbox.at(1) ;
					prefactor *=  (double)bbox.at(2) ;
					prefactor /= ((double)particles.size()+1.0) ;


			// Note: ( N + dN - N0 )^2 - ( N - N0 )^2
			//         = dN^2 + 2*dN*( N - N0 )
			//         = 1 + 2*( N - N0 )
			double dE = 0.5*num_par_kappa * ( 1.0 + 2.0*( (double)particles.size() - (double)num_par_center ) );

			add_particle();

			particles.back()->make_pair_list_this();

			dE += 1.0 * potential.get_particle_energy( particles.back() );

			//cout << "  u(add)   " << potential.particle_energy( particles.back() ) << endl;
			//cout << " dE(total) " << dE << endl << endl;

			if( my_rand() < prefactor*exp( - beta * ( dE - chem_pot ) ) ) {

				update();

				return 1;

			} else {

				particles.erase ( particles.end()-1 ) ;
				update();

				//particles.back()->make_pair_list_all();

				return 0;

			}
		}
	}

	return 0 ;

}



/**
 *   This is the version of the function, where I forgot the N! in the acceptance.
 *   This function can properly be deleted.
 */
int Sim::attempt_a_move_num_particles_OLD ( ) {

	//cout << "end "<< particles.back()->index << endl;

	if(do_num_par_umbrella && potential.do_orientational_order_umbrella){ // TODO Impliment the possiblilty of doing particles removals/insertion togeather with orientational umbrella
		cout << "Error: Performing particles removals/insertion is not implimentet with orientational umbrella (this is a TODO). Exit." << endl;
		exit(21491);
	}

	if( do_num_par_umbrella ) {

		int delta_num_par = 2*(int)floor(my_rand()*2) - 1 ;

		//cout << " delta_num_par = " << delta_num_par << endl;

		if( delta_num_par == -1 ) {	// Remove a particle

			int p = (int)floor(my_rand()*particles.size());

			double dE = -1.0 * potential.get_particle_energy(particles.at(p));

			//cout << "N=" << particles.size() << endl;
			//cout << " u(removal) " << potential.particle_energy(particles.at(p)) << endl << endl;


			// Note: ( N + dN - N0 )^2 - ( N - N0 )^2
			//         = dN^2 + 2*dN*( N - N0 )
			//         = 1 - 2*( N - N0 )
			dE += 0.5*num_par_kappa * ( 1.0 - 2.0*( (double)particles.size() - (double)num_par_center ) );

			if( my_rand() < exp ( - beta * dE ) ) {

				if( particles.size() > 1 ) {
					particles.erase ( particles.begin() + p ) ;
					update();
				}

				return 1;

			} else {

				return 0;

			}

		} else {  // Add a particle

			// Note: ( N + dN - N0 )^2 - ( N - N0 )^2
			//         = dN^2 + 2*dN*( N - N0 )
			//         = 1 + 2*( N - N0 )
			double dE = 0.5*num_par_kappa * ( 1.0 + 2.0*( (double)particles.size() - (double)num_par_center ) );

			add_particle();

			particles.back()->make_pair_list_this();

			dE += 1.0 * potential.get_particle_energy( particles.back() );

			//cout << "  u(add)   " << potential.particle_energy( particles.back() ) << endl;
			//cout << " dE(total) " << dE << endl << endl;

			if( my_rand() < exp( - beta * dE ) ) {

				update();

				return 1;

			} else {

				particles.erase ( particles.end()-1 ) ;
				update();

				//particles.back()->make_pair_list_all();

				return 0;

			}
		}
	}

	return 0 ;

}







void Sim::enable_umbrella_num_particles( double in_num_par_center , double in_num_par_kappa ){
	num_par_center = in_num_par_center ;
	num_par_kappa = in_num_par_kappa;
	do_num_par_umbrella = true ;
}



void Sim::wrap_coordinates(){
	for ( unsigned int i = 0 ; i < particles.size() ; i++ ){
		particles.at(i)->wrap_coordinates();
	}
}

void Sim::initialize_random_number(){
	srand(seed);
}


/**
*    Return random number r,  0.0 <= r < 1.0  (1.0 not included).
*/
double Sim::my_rand(){

	int r = rand();

	while(r==RAND_MAX){
		r=rand();
	}

	return (double)r/((double)(RAND_MAX));

}

