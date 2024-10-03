#!/bin/bash

idir=vmd

odir=movie
mkdir $odir

for i in `awk 'BEGIN{while(i<540){print ++i}}'`;do
	
	ifmt=`printf "%8.8d" $i`

	echo " "
	echo i=$i ifmt=$ifmt
	
	ln -fs $idir/conf$i.vmd vmd.tcl
	vmd -dispdev text -e macro.vmd
	
	convert plot.tga $odir/$ifmt.jpg
done

ffmpeg -b 20M -i $odir/%08d.jpg movie.mp4

