#!/bin/bash -e

if [ $# != 4 ]; then
  echo "Usage: $0 <Lattice input file> <NA> <NB> <NC>"
  exit 1
fi

ifile=$1
NA=$2
NB=$3
NC=$4

awk 'BEGIN{

	# Retieve number of particles in Box
	getline
	N=$1
	#print "N =",N
	print N*'$NA'*'$NB'*'$NC'
	
	getline
	#print "First comment line:",$0

	# Retieve particle coordinates
	for(n=0;n<N;n++){
		getline
		#print "Coordinate line:",$0
		x[n]=$2
		y[n]=$3
		z[n]=$4
		q0[n]=$5
		q1[n]=$6
		q2[n]=$7
		q3[n]=$8
	}
	
	# Retieve lattice vectores
	getline
	#print "Secound comment line:",$0
	for(n=0;n<3;n++){
		getline
		#print "Lattice vector line:",$0
		ux[n]=$1
		uy[n]=$2
		uz[n]=$3
	}
	
	# Write comment line
	printf "Made from %s with %s, bboxX=%f bboxY=%f bboxZ=%f\n","'$ifile'","'$0'",1,1,1

	# Print coordinates in extended xyz format
	for(A=0;A<'$NA';A++){for(B=0;B<'$NB';B++){for(C=0;C<'$NC';C++){
		#print "A B C =",A,B,C
		for(n=0;n<N;n++){
			printf "0 "
			printf "%f ",(x[n]+A)*ux[0]+(x[n]+B)*ux[1]+(x[n]+C)*ux[2]
			printf "%f ",(y[n]+A)*uy[0]+(y[n]+B)*uy[1]+(y[n]+C)*uy[2]
			printf "%f ",(z[n]+A)*uz[0]+(z[n]+B)*uz[1]+(z[n]+C)*uz[2]
			printf "0 0 0 %f %f %f %f \n",q0[n],q1[n],q2[n],q3[n]
		}
	}}}
	
}' $ifile

