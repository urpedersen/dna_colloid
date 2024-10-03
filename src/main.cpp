//============================================================================
// Name        : mc_rigid.cpp
// Author      : Ulf R. Pedersen
// Version     :
// Copyright   :  
// Description :  
//============================================================================

#include <iostream>
using namespace std;

#include <cstdlib>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <sys/time.h>

#include "Sim.h"
#include "helpers_io.h"

int main(int argc, char **arg) {

	   // Initialise

	cout << "Running " << arg[0] << endl ;
	string variables = set_variable_string ( argc , arg ) ;

	   // Create simulation box

	Sim sim ( variables ) ;

	   // Print random numbers (for testing)
	int  print_random_numbers = get_variable ( variables , "print_random_numbers" , 0 ) ;
	if(print_random_numbers>0){
		FILE * ofile_rand     = fopen ( "random_numbers.dat" , "w" ) ;
		for(int i = 0 ; i < print_random_numbers ; i++ ){
			fprintf( ofile_rand , "%f\n" , sim.my_rand() );
		}
		fclose(ofile_rand);
	}

		// Print pair potential
	int print_pair_num_points  = get_variable(variables , "print_pair_num_points"  , 250 ) ;
	int print_angle_num_points = get_variable(variables , "print_angle_num_points" , 180 ) ;

	sim.potential.print_pair_energy     ( "potential_pair.dat"      , print_pair_num_points  ) ;
	sim.potential.print_arm_energy      ( "potential_arm.dat"       , print_pair_num_points  , print_angle_num_points ) ;
	sim.potential.print_arm_energy_angle( "potential_arm_angle.dat" , print_angle_num_points ) ;
	sim.potential.print_arm_energy_distance_angle ( "potential_arm_distance_angle.dat" , print_pair_num_points  , print_angle_num_points ) ;

	   // Load coordinates from file, if it exists.
	if ( get_variable( variables , "load_input_file" , 1 ) ) {
		string filename = "ini.xyz" ;
		sim.load_xyz ( variables , filename ) ;
	}

	   // Scale input coordinates
	{
		cout << "Note: Scale input coordinates (including volume) with coordinate_scaling, that is x_new = coordinate_scaling * x_old (coordinate_scaling=1.0  ->  no scaling)." << endl;
		double coordinate_scaling  = get_variable(variables , "coordinate_scaling"  , 1.0 ) ;
		for (int i = 0 ; i < 3 ; i++ ) sim.bbox.at(i)*=coordinate_scaling ;
	}

	   // Exit if number of particles is 0
	if ( sim.particles.size() == 0 ) {
		cout << "Error: No particles in simulation box, sim.particles.size() = " << sim.particles.size()
		     << ". Exit." << endl ;
		exit ( 7215 ) ;
	}

	   // Enable volume changes (NpT)
	int type_of_volume_move = get_variable ( variables , "type_of_volume_move" , 0   )  ;
	double delta_L_max      = get_variable ( variables , "delta_L_max"         , 0.01 ) ;
	double pressure         = get_variable ( variables , "pressure"            , 0.0 )  ;
	int num_volume_moves_per_step = get_variable ( variables , "num_volume_moves_per_step" , 1 ) ;
	sim.enable_volume_changes ( type_of_volume_move , delta_L_max , pressure );
	if ( sim.do_volume_moves ) {
		cout << "NpT MC simulation is enabled with delta_V_max = " << delta_L_max << " and pressure = " << pressure << endl ;
	}

	   // Umbrella on number of particles
	cout << "[Particle umbrella] = kappa*(N-N_{center})^2" << endl ;
	cout << "umb_num_par sets the number of attempts pr MC step. Set to zero for no particle swat attempts." << endl;
	int umb_num_par = get_variable(variables , "umb_num_par"  , 0 ) ;
	if( umb_num_par > 0 ) {
		sim.enable_umbrella_num_particles (
				get_variable(variables , "num_par_center" , 10.0 ) ,
				get_variable(variables , "num_par_kappa"  , 0.1 )
			) ;
	}

		// Enable pinning of particles
	sim.potential.setup_pinning(variables);

	   // Initialise variables for MC moves
	double max_step_size_translation = get_variable ( variables , "max_step_size_translation" , 0.1 ) ;
	double max_step_size_rotation    = get_variable ( variables , "max_step_size_rotation"    , 0.1 ) ;

	int    num_steps                 = get_variable ( variables , "num_steps"                 , 100000 ) ;
	int    moves_per_mc_step         = get_variable ( variables , "moves_per_mc_step"         , 1 ) ;

	int print_ener_freq              = get_variable ( variables , "print_ener_freq"           , num_steps/10000 ) ;
	if( print_ener_freq < 1 ) {
		print_ener_freq = 1 ;
		cout << "Warning: Reset to print_ener_freq=1" << endl ;
	}

	int print_conf_freq              = get_variable ( variables , "print_conf_freq"           , num_steps/1000  ) ;
	if( print_conf_freq < 1 ) {
		print_conf_freq = 1 ;
		cout << "Warning: Reset to print_conf_freq=1" << endl ;
	}

	int print_vmd_freq              = get_variable ( variables , "print_vmd_freq" , num_steps/10  ) ;
	if( print_vmd_freq < 1 ) {
		print_vmd_freq = 1 ;
		cout << "Warning: Reset to print_vmd_freq=1" << endl ;
	}
	//double vmd_color_min      = get_variable ( variables , "vmd_color_min"  , -6.0 ) ;
	//double vmd_color_max      = get_variable ( variables , "vmd_color_max"  ,  0.0 ) ;
	int print_vmd_counter = 0 ;
	cout << "Note:   type_of_vmd_color=0   ->   particle- & arm energy" << endl << "Note:   type_of_vmd_color=1   ->   orientational order" << endl ;
	int type_of_vmd_color = get_variable ( variables , "type_of_vmd_color" , 0 ) ;

	int accepted_moves_translation   = 0 ;
	int accepted_moves_rotation      = 0 ;
	int accepted_moves_volume        = 0 ;
	int accepted_moves_num_particles = 0 ;
	int next_status_print            = 1 ;

	   // Open files for writing configurations
	FILE * ofile_ener     = fopen ( "ener.dat"    , "w" ) ;
	FILE * ofile_conf_xyz = fopen ( "conf.xyz"    , "w" ) ;
	FILE * ofile_network  = fopen ( "network.dat" , "w" ) ;

	// Build (first) pair list
	sim.update();

	// Wrap coordinates into primary image
	sim.wrap_coordinates();

   // Write simulation state to standard output
	sim.print();


	   // Run MC moves, and write output
	int attempted_moves = 0;

	clock_t start;
	start = clock();

	// Print current Time
	{
		time_t secs=time(0);
		tm *t=localtime(&secs);
		printf("Begin run at %04d-%02d-%02d %02d:%02d:%02d\n",
				t->tm_year+1900,t->tm_mon+1,t->tm_mday,
				t->tm_hour,t->tm_min,t->tm_sec);
	}

	double initial_time = sim.time ;
	for ( int step = 0 ;  step < num_steps ; step++ ) {

		sim.time = initial_time + (double) (step * moves_per_mc_step);

		// Print summary
		if (  (step==next_status_print)  |  (step==0)  ) {

			cout << "Stat: step " << step
				 << " ("          << (double)step/(double)num_steps*100.0                  << "%), "
				 << "t="		  << sim.time << " "
				 << "U="          << sim.potential.get_internal_energy() << " "
				 << "V="          << sim.get_volume() << " "
				 << "N="          << sim.particles.size() << " "
				 << "Q="		  << sim.potential.get_orientational_order() << " "
				 << "UQ="		  << sim.potential.get_orientational_order_energy() << " "
				 << "Upin="		  << sim.potential.get_pinning_energy() << " "
				 << "moves="      << attempted_moves << " "
				 << "acc_trans="  << (double)accepted_moves_translation/(double)attempted_moves*100.0 << "%, "
				 << "acc_rot="    << (double)accepted_moves_rotation/(double)attempted_moves*100.0    << "%, "
				 << "acc_vol="    << (double)accepted_moves_volume/(double)step/(double)num_volume_moves_per_step*100.0    << "%, "
				 << "acc_par="    << (double)accepted_moves_num_particles/(double)step/(double)umb_num_par*100.0    << "%, "
				 << "nl_upd="     << sim.neighbour_list_updates << " "
				 << "t_wall="     << ( clock() - start ) / (double)CLOCKS_PER_SEC << " "
				 << "t_est="     << ( clock() - start ) / (double)CLOCKS_PER_SEC / (double)step * (double)num_steps << " "
				 << "vmd_prints=" << print_vmd_counter << " "
				 << endl;

			next_status_print *= 2 ;

		}

		// Print energy
		if ( step % print_ener_freq == 0 ) {
			fprintf ( ofile_ener , "%s %f\n", sim.get_ener_string().c_str() ,  ( clock() - start ) / (double)CLOCKS_PER_SEC ) ;
		}

		// Print configuration, xyz and restart and network
		if ( step % print_conf_freq == 0 ) {

			sim.wrap_coordinates() ;
			fprintf ( ofile_conf_xyz , "%s\n", sim.get_coordinates_xyz().c_str() ) ;
			fprintf ( ofile_network , "%s\n", sim.get_coordinates_network().c_str() ) ;

			FILE *    ofile_conf_restart = fopen ( "restart.xyz" , "w" ) ;
			fprintf ( ofile_conf_restart , "%s\n", sim.get_coordinates_xyz().c_str() ) ;
			fclose  ( ofile_conf_restart );

		}

		// Print configuration, vmd
		if ( step % print_vmd_freq == 0 ) {
			sim.wrap_coordinates () ;

			stringstream filename;
			filename << "conf" << print_vmd_counter << ".vmd" ;
			FILE * ofile_conf_vmd     = fopen ( filename.str().c_str() , "w" ) ;

			sim.wrap_coordinates () ;

			fprintf ( ofile_conf_vmd , "%s\n", sim.get_coordinates_vmd(type_of_vmd_color).c_str() ) ;

			print_vmd_counter++ ;

			fclose ( ofile_conf_vmd ) ;

		}

		// Make MC moves
		{
			int is_accepted;

			// Translation and rotation
			for( int i = 0 ; i < (int)sim.particles.size()*moves_per_mc_step ; i ++ ) {
				attempted_moves++;

				is_accepted = sim.attempt_a_move_translation(&max_step_size_translation);
				accepted_moves_translation += is_accepted;

				if ( sim.potential.get_arm_eps()!=0.0 )
					is_accepted = sim.attempt_a_move_rotation(&max_step_size_rotation);
				else
					is_accepted = 0 ;

				accepted_moves_rotation += is_accepted;
			}

			// Volume
			for ( int i = 0 ; i < num_volume_moves_per_step ; i++ ) {
				is_accepted = sim.attempt_a_move_volume();
				accepted_moves_volume += is_accepted;
			}

			// Number of particles
			for ( int i = 0 ; i<umb_num_par ; i++ ) {
				is_accepted = sim.attempt_a_move_num_particles();
				accepted_moves_num_particles += is_accepted;
			}
		}
	}
	if(num_steps>0) sim.time = initial_time + (double) ( num_steps * moves_per_mc_step ) ;

	cout << "Done with " << num_steps << " MC steps," <<  endl;
	cout << "   corresponding to dt=" << num_steps*moves_per_mc_step << "." <<  endl;

	// Print current Time
	{
		time_t secs=time(0);
		tm *t=localtime(&secs);
		printf("End run at %04d-%02d-%02d %02d:%02d:%02d\n",
				t->tm_year+1900,t->tm_mon+1,t->tm_mday,
				t->tm_hour,t->tm_min,t->tm_sec);
	}

	   // Report acceptance rate
	double run_time = ( clock() - start ) / (double)CLOCKS_PER_SEC;
	cout << "Run time of " << run_time << " seconds "
		 << " = " << run_time/60.0 << " minutes "
		 << " = " << run_time/3600.0 << " hours "
		 << " = " << run_time/3600.0/24 << " days "
		 << " = " << run_time/3600.0/24/7 << " weeks "
		 << " = " << run_time/3600.0/365.242199*12 << " months "
		 << " = " << run_time/3600.0/365.242199 << " years, "
		 << " = " << run_time/3600.0/365.242199/9814072356.0 << " earth life-times, " << endl
		 << "  corresponding to " << (double)num_steps/run_time << " MC steps per second," << endl
	     << "  or " << (double)attempted_moves/run_time << " particle MC steps per second." << endl;

	cout << "Accepted " << accepted_moves_translation
		 << " out of "  << attempted_moves
		 << " attempted translation moves (" << (double)accepted_moves_translation/(double)attempted_moves*100.00
		 <<  "%)" << endl;

	cout << "Accepted " << accepted_moves_rotation
		 << " out of "  << attempted_moves
		 << " attempted rotation moves ("    << (double)accepted_moves_rotation/(double)attempted_moves*100.00
		 <<  "%)" << endl;

	cout << "Accepted " << accepted_moves_volume
		 << " out of "  << num_steps*num_volume_moves_per_step
		 << " attempted volume moves ("    << (double)accepted_moves_volume/(double)num_steps/(double)num_volume_moves_per_step*100.00
		 <<  "%)" << endl;

	cout << "Accepted " << accepted_moves_num_particles
		 << " out of "  << num_steps*umb_num_par
		 << " attempted particle insertions/removals ("    << (double)accepted_moves_num_particles/(double)num_steps/(double)umb_num_par*100.00
		 <<  "%)" << endl;

	cout << "Made " << sim.neighbour_list_updates << " neighbour list updates " << endl
	     << "  corresponding to " << (double)sim.neighbour_list_updates/(double)num_steps << " per MC step." << endl
         << "  corresponding to " << (double)sim.neighbour_list_updates/(double)attempted_moves << " per MC particle step." << endl ;

	sim.print();

	// Write final coordinates to end.xyz

	cout << "Write final configuration to end.xyz and end.vmd." << endl ;

	sim.wrap_coordinates() ;

	FILE *    ofile_conf_end = fopen ( "end.xyz" , "w" ) ;
	fprintf ( ofile_conf_end , "%s", sim.get_coordinates_xyz().c_str() ) ;
	fclose  ( ofile_conf_end ) ;

	FILE *    ofile_vmd_end = fopen ( "end.vmd" , "w" ) ;
	fprintf ( ofile_vmd_end , "%s", sim.get_coordinates_vmd(type_of_vmd_color).c_str() ) ;
	fclose  ( ofile_vmd_end ) ;

	   // Write simulation state to standard output
	//sim.print();

	   // Finalize

	fclose ( ofile_ener     ) ;
	fclose ( ofile_conf_xyz ) ;
	fclose ( ofile_network ) ;

	cout << "Variable string: " << variables << endl ;
	cout << "Happy ending of  "  << arg[0]    << endl ;

	return 0 ;
}
