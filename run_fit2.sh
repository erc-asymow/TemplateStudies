#!/bin/sh

GRID=$1

for CORRX in 10
do
    for CORRY in 4
    do
	for A0X in 2
	do		    
	    for A0Y in 2
	    do		    
		echo 'Now doing ', $CORRX, $CORRY
		./fit --nevents=100000000 --tag=$GRID --run=closure\
                --degs_corr_x=$CORRX --degs_corr_y=$CORRY\
		--degs_A0_x=$A0X --degs_A0_y=$A0Y\
                --j0
	    done
	done
    done
done
