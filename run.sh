#!/bin/sh

GRID=$1

for X in 9 10
do
    for Y in 6 7 8
    do
	echo 'Now doing ', $X, $Y
	./main -1 $GRID closure $X $Y
	./main 10000000 $GRID closure $X $Y
    done
done
