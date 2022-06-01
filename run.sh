#!/bin/sh

GRID=$1

for X in 12
do
    for Y in 6
    do
	for CORRX in 3
	do
	    for CORRY in 3
	    do
		for A0X in 2
		do		    
		    for A0Y in 2
		    do		    
			echo 'Now doing ', $X, $Y, $CORRX, $CORRY
			#./main  --nevents=100000000 --tag=$GRID --run=closure --degs_pdf_x=$X --degs_pdf_y=$Y --degs_corr_x=$CORRX --degs_corr_y=$CORRY
			./jac --nevents=10000000 --tag=$GRID --run=closure --degs_pdf_x=$X --degs_pdf_y=$Y --degs_corr_x=$CORRX --degs_corr_y=$CORRY\
		        --degs_A0_x=$A0X --degs_A0_y=$A0Y\
		        --toyTF2_corr\
                        --normalize_pdfx\
			--normalize_pdfy			
			#--do_absy 
		    done
		done
	    done
	done
    done
done
