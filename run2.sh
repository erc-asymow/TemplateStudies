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
		echo 'Now doing ', $CORRX, $CORRY, $A0X, $A0Y
		./jac2 --nevents=1000000000 --tag=$GRID --run=closure\
                --degs_corr_x=$CORRX --degs_corr_y=$CORRY\
		--degs_A0_x=$A0X --degs_A0_y=$A0Y\
                --toyTF2_corr
		#--fit_qt_y
		#--do_absy
	    done
	done
    done
done
