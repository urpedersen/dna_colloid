/**
 * Particle.cpp
 *
 *  Created on: Oct 6, 2011
 *      Author: urp
 */

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
#include "Arm.h"

/**
 * Construct particles
 *
 * Input
 * 	u_local coordinate
 *  y coordinate
 *  u_local coordinate
 *
 * 	q0 quaternian coordinate
 *  q1 quaternian coordinate
 *  q2 quaternian coordinate
 *  q3 quaternian coordinate
 */
Particle::Particle(
		unsigned int in_index,
		int in_type,
		vector<double> *in_bbox,
		double in_x,
		double in_y,
		double in_z,
		int in_img_x,
		int in_img_y,
		int in_img_z,
		double in_q0,
		double in_q1,
		double in_q2,
		double in_q3,
		vector<Particle*> *in_particles,
		double *in_cut_distance,
		double *in_skin_distance,
		int *in_neighbour_list_updates
	) {

	index = in_index;
	type = in_type;

	bbox = in_bbox;

	// Coordinates
	x.clear();
	x.push_back(in_x); // Make this variable a static array
	x.push_back(in_y);
	x.push_back(in_z);

	// Coordinates
	img.clear();
	img.push_back(in_img_x); // Make this variable a static array
	img.push_back(in_img_y);
	img.push_back(in_img_z);

	// Quaternian Coordinates
	q.clear();
	q.push_back(in_q0);  // Make this variable a static array
	q.push_back(in_q1);
	q.push_back(in_q2);
	q.push_back(in_q3);

	rotation_matrix.assign(9,0.0);	// Make this variable a static array

	// Pair list set-up
	particles     = in_particles;
	cut_distance  = in_cut_distance;
	skin_distance = in_skin_distance;
	dx_since_last_pair_list_update.clear();
	dx_since_last_pair_list_update.assign(3,*skin_distance+1.0);
	neighbour_list_updates = in_neighbour_list_updates;

	// Set-up arms
	set_arms("arms.xyz");

	update_orientation();

	// Note. pair list should be build elsewhere
	//make_pair_list_all();

	// Setup parameters for pinning potential
	//enable_pinning = false ;

	pinning_x.clear ( ) ;
	pinning_x.push_back ( in_x ) ;
	pinning_x.push_back ( in_y ) ;
	pinning_x.push_back ( in_z ) ;

	pinning_img.clear ( ) ;
	pinning_img.push_back ( in_img_x ) ;
	pinning_img.push_back ( in_img_y ) ;
	pinning_img.push_back ( in_img_z ) ;

	pinning_kappa = 0.0 ;

}

Particle::~Particle() {
	// TODO Auto-generated destructor stub
}


/**
 *  Update pair list of this particle.
 *  Note: List of particles in neighborhood should be build elsewhere.
 */
void Particle::make_pair_list_this(){

	pairs.clear ( ) ;

	int nl_length = 0 ;

	for(unsigned int i = 0 ; i < particles->size()  ; i++ ) {
		double d = r_mi(particles->at(i)) ;

		if( d < (*cut_distance + *skin_distance) ) {

			//cout << particles->at(i)->index << " " << index << endl;

			if( particles->at(i)->index != index) {

				pairs.push_back(particles->at(i));

				nl_length++;

				//cout << pairs.size() << endl;

			}
		}
	}

	//cout << nl_length << endl;

	// Reset travel distance
	dx_since_last_pair_list_update.clear();
	dx_since_last_pair_list_update.assign( 3 , 0.0 );

}

void Particle::make_pair_list_all(){

	//cout << "Make pair list neighbour_list_updates=" << (*neighbour_list_updates)++ << endl;
	(*neighbour_list_updates)++;

	for(unsigned int i = 0 ; i < particles->size() ; i++ )
		particles->at(i)->make_pair_list_this();
}

/*
 * Pair list is old if particle have mode a half skin or more since last pair list update.
 */
bool Particle::pair_list_is_old ( ) {
	double d = dx_since_last_pair_list_update.at(0)*dx_since_last_pair_list_update.at(0) ;
	d += dx_since_last_pair_list_update.at(1)*dx_since_last_pair_list_update.at(1) ;
	d += dx_since_last_pair_list_update.at(2)*dx_since_last_pair_list_update.at(2) ;
	d = sqrt(d) ;
	return d >= (*skin_distance)*0.5 ;
}

/**
 * Assign arms to particles
 */
void Particle::set_arms(string filename){

	arms.clear();

	ifstream ifile (filename.c_str());

	if (ifile.is_open()) {
		string line;

		getline ( ifile , line ) ; // Number of atoms line
		int arms_to_load = atoi ( line.c_str() ) ;

		getline ( ifile , line ) ; // Header line

		for ( int i=0 ; i < arms_to_load ; i++ ) {

			getline (ifile,line);

			char * pEnd;
			int    i0;
			double d0, d1, d2;

			i0  = strtol (line.c_str(),&pEnd,10);
			d0  = strtod (pEnd,&pEnd);
			d1  = strtod (pEnd,&pEnd);
			d2  = strtod (pEnd,&pEnd);

			arms.push_back ( new Arm ( d0 , d1 , d2 ) );

		}
	}else{
		cout << "Error: Unable to open file with arm coordinates (" << filename << "). Exit." << endl;
		exit (1);
	}

	/*arms.push_back ( new Arm(  1.0 , 0.0 , 0.0 ) ) ;
	arms.push_back ( new Arm( -1.0 , 0.0 , 0.0 ) ) ;
	arms.push_back ( new Arm(  0.0 , 1.0 , 0.0 ) ) ;
	arms.push_back ( new Arm(  0.0 ,-1.0 , 0.0 ) ) ;
	arms.push_back ( new Arm(  0.0 , 0.0 , 1.0 ) ) ;
	arms.push_back ( new Arm(  0.0 , 0.0 ,-1.0 ) ) ; */

}

/**
 * Translate particle in dim dimension by step.
 * Note: image is not wrapped.
 */
void Particle::translate_step(int dim,double step){
	double step_reduced = step/bbox->at(dim);
	x.at(dim) += step_reduced;

	// Test Pair list
	dx_since_last_pair_list_update.at(dim)+=step ;
	if(pair_list_is_old()) make_pair_list_all()  ;
}

/**
 * Rotate
 */
void Particle::rotate_step(vector<double> step){

	for ( unsigned int i = 0 ; i < 4 ; i++ ) {
		q.at(i) += step.at(i);
	}

	update_orientation();

}

/**
 * Wrap the position vector, u_local, into unit cube. Information about image configuration is transfered to
 */
void Particle::wrap_coordinates ( ) {
	for(int i = 0 ; i < 3 ; i++ ) {
		//cout << "Before: u_local=" << u_local.at(i) << " img=" << img.at(i) << "\n";
		int di = (int)floor(x.at(i));
		img.at(i) += di;
		x.at(i) -= (double)di;
		//cout << "After: u_local=" << u_local.at(i) << " img=" << img.at(i) << "\n\n";
	}
}

/**
 * Return  position in the minimum image of dim dimentaion
 *
 * Input: dimention wher 0 = u_local, 1 = y and 2 = z .
 */
double Particle::get_min_img(int dim){
	return x.at(dim)*bbox->at(dim);
}

/**
 * Return X position in the minimum image
 */
double Particle::get_min_img_x(){
	return x.at(0)*bbox->at(0);
}

/**
 * Return Y position in the minimum image
 */
double Particle::get_min_img_y(){
	return x.at(1)*bbox->at(1);
}

/**
 * Return Z position in the minimum image
 */
double Particle::get_min_img_z(){
	return x.at(2)*bbox->at(2);
}

/**
 * Return position in the minimum image of arm in dimention din
 */
double Particle::get_min_img( int dim , Arm* arm ) {
	return x.at(dim)*bbox->at(dim) + arm->u_global.at(dim) ;
}

/**
 * Return X position in the minimum image of arm
 */
double Particle::get_min_img_x( Arm* arm ) {
	return x.at(0)*bbox->at(0) + arm->u_global.at(0) ;
}

/**
 * Return Y position in the minimum image of arm
 */
double Particle::get_min_img_y( Arm* arm ) {
	return x.at(1)*bbox->at(1) + arm->u_global.at(1) ;
}

/**
 * Return Y position in the minimum image of arm
 */
double Particle::get_min_img_z( Arm* arm ) {
	return x.at(2)*bbox->at(2) + arm->u_global.at(2) ;
}

double Particle::get_dx_img(Particle *p,int d){

	double dx = p->x.at(d) - x.at(d);

	dx -= round(dx);

	dx *= bbox->at(d);

	return dx;
}

/**
 * Squared distance using minimum image conversion
 */
double Particle::r2_mi(Particle *p_ref) {
	double out=0.0 ;
	for (int i = 0 ; i < 3 ; i ++) {
		double x_tmp = x.at(i) - p_ref->x.at(i);
		x_tmp = x_tmp-round(x_tmp);
		x_tmp = x_tmp*bbox->at(i);
		out += x_tmp*x_tmp;
	}
	return out;
}



/**
 * Distance using minimum image conversion
 */
double Particle::r_mi(Particle *p_ref){
	double out = sqrt( r2_mi(p_ref) );
	return out;
}


/**
 * Squared distance to pinning position (using minimum position)
 */
double Particle::r2_pinning ( ) {

	double out=0.0 ;

	for (int i = 0 ; i < 3 ; i++ ) {
		double x_tmp  = x.at(i)         + (double)img.at(i)         ;
		       x_tmp -= pinning_x.at(i) + (double)pinning_img.at(i) ;
		       x_tmp *= bbox->at(i) ;
		out += x_tmp*x_tmp ;
	}

	return out;

}


/**
 *
 */
void Particle::get_vector_in_global_frame ( vector<double>* o , vector<double>* i ) {

	// Consistency test
	/*if( o->size() !=3  or  i->size() !=3  or  q.size()!=4 ){
		cout << "Error: In void Particle::get_vector_in_global_frame ( vector<double>* o , vector<double>* i ) " << endl;
		cout <<	" Wrong size of vectors pased to void Particle::get_vector_in_global_frame ( vector<double>* o , vector<double>* i ) " << endl;
		cout << "        o->size() = " << o->size() << " (should be 3)" << endl;
		cout << "        i->size() = " << i->size() << " (should be 3)" << endl;
		cout << "        q.size()  = " << q.size() <<  " (should be 4)" << endl;
		cout << "Exit." << endl;
		exit(0);
	}*/

	//cout << "use rot " << endl;

	o->at(0)  = rotation_matrix[0] * i->at(0);
	o->at(0) += rotation_matrix[1] * i->at(1);
	o->at(0) += rotation_matrix[2] * i->at(2);

	o->at(1)  = rotation_matrix[3] * i->at(0);
	o->at(1) += rotation_matrix[4] * i->at(1);
	o->at(1) += rotation_matrix[5] * i->at(2);

	o->at(2)  = rotation_matrix[6] * i->at(0);
	o->at(2) += rotation_matrix[7] * i->at(1);
	o->at(2) += rotation_matrix[8] * i->at(2);

	/*
	o->at(0)  = rotation_matrix.at(0) * i->at(0);
	o->at(0) += rotation_matrix.at(1) * i->at(1);
	o->at(0) += rotation_matrix.at(2) * i->at(2);

	o->at(1)  = rotation_matrix.at(3) * i->at(0);
	o->at(1) += rotation_matrix.at(4) * i->at(1);
	o->at(1) += rotation_matrix.at(5) * i->at(2);

	o->at(2)  = rotation_matrix.at(6) * i->at(0);
	o->at(2) += rotation_matrix.at(7) * i->at(1);
	o->at(2) += rotation_matrix.at(8) * i->at(2);
	*/


	// Construct rotation matrix
	/*
	double q0q0=q.at(0)*q.at(0);
	double q1q1=q.at(1)*q.at(1);
	double q2q2=q.at(2)*q.at(2);
	double q3q3=q.at(3)*q.at(3);

	double q0q1=q.at(0)*q.at(1);
	double q0q2=q.at(0)*q.at(2);
	double q0q3=q.at(0)*q.at(3);

	double q1q2=q.at(1)*q.at(2);
	double q1q3=q.at(1)*q.at(3);

	double q2q3=q.at(2)*q.at(3);

	o->at(0)  =     ( q0q0 + q1q1 - q2q2 - q3q3 ) * i->at(0);
	o->at(0) += 2.0*( q1q2 + q0q3 )               * i->at(1);
	o->at(0) += 2.0*( q1q3 - q0q2 )               * i->at(2);

	o->at(1)  = 2.0*( q1q2 - q0q3 )               * i->at(0);
	o->at(1) +=     ( q0q0 - q1q1 + q2q2 - q3q3 ) * i->at(1);
	o->at(1) += 2.0*( q2q3 + q0q1 )               * i->at(2);

	o->at(2)  = 2.0*( q1q3 + q0q2 )               * i->at(0);
	o->at(2) += 2.0*( q2q3 - q0q1 )               * i->at(1);
	o->at(2) +=     ( q0q0 - q1q1 - q2q2 + q3q3 ) * i->at(2);
	*/

	// Testing implementation
	/*
	cout << "vector_old = [ " << i->at(0) << " " << i->at(1) << " " << i->at(2) << " ]" << endl;
	cout << "    [ " << rotation_matrix.at(0) << " " << rotation_matrix.at(1) << " " << rotation_matrix.at(2) << " ]" << endl;
	cout << "A = [ " << rotation_matrix.at(3) << " " << rotation_matrix.at(4) << " " << rotation_matrix.at(5) << " ]" << endl;
	cout << "    [ " << rotation_matrix.at(6) << " " << rotation_matrix.at(7) << " " << rotation_matrix.at(8) << " ]" << endl;
	cout << "vector_new = [ " << o->at(0) << " " << o->at(1) << " " << o->at(2) << " ]" << endl;
	*/
	/*
	double len2_i = i->at(0)*i->at(0) + i->at(1)*i->at(1) + i->at(2)*i->at(2);
	double len2_o = o->at(0)*o->at(0) + o->at(1)*o->at(1) + o->at(2)*o->at(2);
	if(len2_i-len2_o > 1e-14 ) {
		cout << len2_i-len2_o << " ";
		cout << "r2(i) = " << len2_i << " ";
		cout << "r2(o) = " << len2_o << " ";
		cout << "r2(q) = " << q.at(0)*q.at(0) + q.at(1)*q.at(1) + q.at(2)*q.at(2) + q.at(3)*q.at(3) << endl;
	}
	*/
}




void Particle::update_orientation ( ) {

	normalize_q ( ) ;

	construct_rotation_matrix ( ) ;

	set_arm_orientations ( ) ;

}

/*number_of_particles
 * Normalize quaternian vector
 */
void Particle::normalize_q() {

	double tmp = 0.0;

	for(int i=0;i<4;i++){
		tmp+=q.at(i)*q.at(i);
	}

	tmp = sqrt(tmp);

	for(int i=0;i<4;i++){
		q.at(i) = q.at(i)/tmp;
	}

}

/**
* Rotate in_vector with the quatenian rotation matrix.
*
* Input: vector in local particle frame (i vector, second input).
*
* Return: a vector in global frame (to o vector, first input).
*
* Both vectors are 3D.
*
* Note: This function use the quaternian rotation matrix, A, see e.q. [Allen & Tildesley, Computersimulations of Liquids, section 3.3.1, page 88]:
*
*   o = A*i, where
*
*       [  q0*q0 + q1*q1 - q2*q2 - q3*q3  ,  2*( q1*q2 + q0*q3 )            ,  2*( q1*q3 - q0*q2 )            ]
*       [                                                                                                     ]
*   A = [  2*( q1*q2 - q0*q3 )            ,  q0*q0 - q1*q1 + q2*q2 - q3*q3  ,  2*( q2*q3 + q0*q1 )            ]
* 		[                                                                                                     ]
*       [  2*( q1*q3 + q0*q2 )            ,  2*( q2*q3 - q0*q1 )            ,  q0*q0 - q1*q1 - q2*q2 - q3*q3  ]
*
*/
void Particle::construct_rotation_matrix(){

	//cout << "make rot " << endl;

	// Construct rotation matrix
	double q0q0=q.at(0)*q.at(0);
	double q1q1=q.at(1)*q.at(1);
	double q2q2=q.at(2)*q.at(2);
	double q3q3=q.at(3)*q.at(3);

	double q0q1=q.at(0)*q.at(1);
	double q0q2=q.at(0)*q.at(2);
	double q0q3=q.at(0)*q.at(3);

	double q1q2=q.at(1)*q.at(2);
	double q1q3=q.at(1)*q.at(3);

	double q2q3=q.at(2)*q.at(3);

	rotation_matrix.at(0) = (       q0q0 + q1q1 - q2q2 - q3q3  ) ;
	rotation_matrix.at(1) = ( 2.0*( q1q2 + q0q3 )              ) ;
	rotation_matrix.at(2) = ( 2.0*( q1q3 - q0q2 )              ) ;

	rotation_matrix.at(3) = ( 2.0*( q1q2 - q0q3 )              ) ;
	rotation_matrix.at(4) = (       q0q0 - q1q1 + q2q2 - q3q3  ) ;
	rotation_matrix.at(5) = ( 2.0*( q2q3 + q0q1 )              ) ;

	rotation_matrix.at(6) = ( 2.0*( q1q3 + q0q2 )              ) ;
	rotation_matrix.at(7) = ( 2.0*( q2q3 - q0q1 )              ) ;
	rotation_matrix.at(8) = (       q0q0 - q1q1 - q2q2 + q3q3  ) ;

	/*
	cout << "q = [ " << q.at(0) << " " << q.at(1) << " " << q.at(2) << " " << q.at(3) << " ]" << endl ;
	cout << "    [ " << rotation_matrix.at(0) << " " << rotation_matrix.at(1) << " " << rotation_matrix.at(2) << " ]" << endl;
	cout << "A = [ " << rotation_matrix.at(3) << " " << rotation_matrix.at(4) << " " << rotation_matrix.at(5) << " ]" << endl;
	cout << "    [ " << rotation_matrix.at(6) << " " << rotation_matrix.at(7) << " " << rotation_matrix.at(8) << " ]" << endl;
	*/
}

void Particle::set_arm_orientations ( ) {

	//cout << "   arms.size() = " << arms.size() << endl;

	for ( unsigned int i = 0 ; i < arms.size() ; i++ ){
		get_vector_in_global_frame ( &arms.at(i)->u_global , &arms.at(i)->u_local );
	}

}




/**
 * Print particle state to standard output.
 */
void Particle::print(){
	cout << "Particle " << index << " : ";

	cout << "type=" << type << " ";

	cout << "x=[";
	for(unsigned int i = 0 ; i<x.size()-1 ; i ++){
		cout << x.at(i) << ",";
	}
	cout << x.at(x.size()-1) << "] ";

	cout << "x_real=[";
	for(unsigned int i = 0 ; i<x.size()-1 ; i ++){
		cout << x.at(i)*bbox->at(i) << ",";
	}
	cout << x.at(x.size()-1)*bbox->at(x.size()-1) << "] ";

	cout << "img=[";
	for(unsigned int i = 0 ; i<img.size()-1 ; i ++){
		cout << img.at(i) << ",";
	}
	cout << img.at(img.size()-1) << "] ";

	cout << "q=[";
	for(unsigned int i = 0 ; i<q.size()-1 ; i ++){
		cout << q.at(i) << ",";
	}
	cout << q.at(q.size()-1) << "]" << endl;

}


