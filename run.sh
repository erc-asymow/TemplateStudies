#!/bin/sh

GRID=$1

for X in 10
do
    for Y in 6
    do
	for CORRX in 6
	do
	    for CORRY in 6
	    do
		echo 'Now doing ', $X, $Y, $CORRX, $CORRY
		#./main      --nevents=-1 --tag=$GRID --run=closure --degs_pdf_x=$X --degs_pdf_y=$Y --degs_corr_x=$CORRX --degs_corr_y=$CORRY --do_absY=1
		./main --nevents=100000000 --tag=$GRID --run=closure --degs_pdf_x=$X --degs_pdf_y=$Y --degs_corr_x=$CORRX --degs_corr_y=$CORRY --toyTF2_corr --normalize_pdfx --do_absY
	    done
	done
    done
done
