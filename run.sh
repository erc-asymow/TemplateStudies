#!/bin/sh

GRID=$1

for X in 10
do
    for Y in 4
    do
	for CORRX in 4
	do
	    for CORRY in 4 
	    do
		echo 'Now doing ', $X, $Y, $CORRX, $CORRY
		./main      --nevents=-1 --tag=$GRID --run=closure --degs_pdf_x=$X --degs_pdf_y=$Y --degs_corr_x=$CORRX --degs_corr_y=$CORRY --do_absY=1
		./main --nevents=1000000 --tag=$GRID --run=closure --degs_pdf_x=$X --degs_pdf_y=$Y --degs_corr_x=$CORRX --degs_corr_y=$CORRY --do_absY=1
	    done
	done
    done
done
