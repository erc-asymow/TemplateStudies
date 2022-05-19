#!/bin/sh

GRID=$1


for X in 10
do
    for Y in 4 6 8
    do
	for CORRX in 6
	do
	    for CORRY in 4 6 8 
	    do
		echo 'Now doing ', $X, $Y, $CORRX, $CORRY
 		./main        -1 $GRID closure $X $Y $CORRX $CORRY
		./main  10000000 $GRID closure $X $Y $CORRX $CORRY
	    done
	done
    done
done
