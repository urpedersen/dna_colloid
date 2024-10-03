#!/bin/bash
#  By Ulf R. Pedersen (www.urp.dk)

if [ $# != 4 ]; then
   echo "Usage: $0 <filename of xyz-file> <nX> <nY> <nZ> <X> <Y> <Z> <scale factor>"
   exit 1
fi

# Input file
ifile=$1

# Times to copy
nx=$2
ny=$3
nz=$4

# Box bondaries
X=$5
Y=$6
Z=$7

# Scale factor
s=$8

head -n 1 $ifile | awk '{print $1*'$nx'*'$ny'*'$nz'}'

echo Copied configuration in $ifile nx=$nx ny=$ny nz=$nz times. Initial bboxX=$X bboxY=$Y bboxZ=$Z.

awk '{if(NR>2){
  for(ix=1;ix<='$nx';ix++){
    for(iy=1;iy<='$ny';iy++){
      for(iz=1;iz<='$nz';iz++){
        print $1,($2+(ix-1)*'$X')*'$s',($3+(iy-1)*'$Y')*'$s',($4+(iz-1)*'$Z')*'$s',0,0,0,$8,$9,$10,$11
      }
    }
  }
}}' $ifile
