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
 		./main        -1 $GRID closure $X $Y $CORRX $CORRY
		./main   100000 $GRID closure $X $Y $CORRX $CORRY
	    done
	done
    done
done
