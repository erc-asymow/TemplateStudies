#!/bin/sh

GRID=$1

for X in 12
do
    for Y in 6
    do
	for CORRX in 4
	do
	    for CORRY in 4
	    do
		for A0X in 2
		do		    
		    for A0Y in 2
		    do		    
			echo 'Now doing ', $X, $Y, $CORRX, $CORRY
			./fit --nevents=-1 --tag=$GRID --run=closure --degs_pdf_x=$X --degs_pdf_y=$Y --degs_corr_x=$CORRX --degs_corr_y=$CORRY\
		        --degs_A0_x=$A0X --degs_A0_y=$A0Y\
                        --j0 --j1 --j2
			#--j1 --j2
                        #--j0
                        #--j0
		    done
		done
	    done
	done
    done
done
