#!/bin/sh

GRID=$1

for CORRX in 10
do
    for CORRY in 4 
    do
	for A0X in 3
	do		    
	    for A0Y in 4
	    do		    
		for A1X in 3
		do		    
		    for A1Y in 4
		    do		    
			for A2X in 3
 			do		    
			    for A2Y in 4
			    do		    
				for A3X in 2
				do		    
				    for A3Y in 4
				    do	
					for A4X in 2
					do		    
					    for A4Y in 4
					    do		    
						echo 'Now doing ', $CORRX, $CORRY, $A0X, $A0Y, $A1X, $A1Y, $A2X, $A2Y, $A3X, $A3Y, $A4X, $A4Y
						./jac2 --nevents=100000000 --tag=$GRID --run=closure\
                                                --degs_corr_x=$CORRX --degs_corr_y=$CORRY\
                                 		--degs_A0_x=$A0X --degs_A0_y=$A0Y\
                                 		--degs_A1_x=$A1X --degs_A1_y=$A1Y\
                                 		--degs_A2_x=$A2X --degs_A2_y=$A2Y\
                                 		--degs_A3_x=$A3X --degs_A3_y=$A3Y\
                                 		--degs_A4_x=$A4X --degs_A4_y=$A4Y\
                                                --toyTF2_corr
						#--fit_qt_y
					    done
					done
				    done
				done
			    done
			done
		    done
		done
	    done
	done
    done
done
