#!/bin/bash

odir=movie
mkdir $odir

for i in `awk 'BEGIN{while(i<4){print ++i}}'`;do
	
	ifmt=`printf "%8.8d" $i`

	echo " "
	echo i=$i ifmt=$ifmt
	
	ln -fs conf$i.vmd vmd.tcl
	vmd -dispdev text -e macro.vmd
	
	convert vmdscene.pov.tga $odir/$ifmt.jpg
done

ffmpeg -i $odir/%08d.jpg movie.mp4

