#!/bin/sh

ADD="_ADDMC_ULA0A1A2A3A4"
ADDflag="--j0 --j1 --j2 --j3 --j4 --add_MC_uncert"

#ADD="_ULA0A1A2A3A4"
#ADDflag="--j0 --j1 --j2 --j3 --j4"

NEVENTS=100000000

echo $ADD $ADDflag

#./fit --nevents=${NEVENTS} --tag=DEBUGNEWA3_100M --run=grid --post_tag=DEBUG${ADD} --degs_corr_x=8 --degs_corr_y=6 --degs_A0_x=1 --degs_A0_y=1 --degs_A1_x=1 --degs_A1_y=1 --degs_A2_x=1 --degs_A2_y=1 --degs_A3_x=1 --degs_A3_y=1 --degs_A4_x=1 --degs_A4_y=1 --jUL   $ADDflag
#./fit --nevents=${NEVENTS} --tag=NEWA3_1G --run=full --post_tag=DEBUG${ADD} --degs_corr_x=10 --degs_corr_y=4 --degs_A0_x=3 --degs_A0_y=2 --degs_A1_x=3 --degs_A1_y=3 --degs_A2_x=3 --degs_A2_y=2 --degs_A3_x=3 --degs_A3_y=3 --degs_A4_x=3 --degs_A4_y=3 --jUL   $ADDflag
#./fit --nevents=${NEVENTS} --tag=NEWA3_1G --run=full --post_tag=DEBUG${ADD} --degs_corr_x=10 --degs_corr_y=4 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=3 --degs_A1_y=3 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=3 --degs_A3_y=3 --degs_A4_x=3 --degs_A4_y=3 --jUL   $ADDflag
#./fit --nevents=${NEVENTS} --tag=NEWA3_1G --run=full --post_tag=DEBUG${ADD} --degs_corr_x=10 --degs_corr_y=4 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=3 --degs_A1_y=3 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=3 --degs_A3_y=3 --degs_A4_x=1 --degs_A4_y=3 --jUL   $ADDflag
#./fit --nevents=${NEVENTS} --tag=NEWA3_1G --run=full --post_tag=DEBUG${ADD} --degs_corr_x=10 --degs_corr_y=4 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=3 --degs_A1_y=3 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=3 --degs_A3_y=3 --degs_A4_x=2 --degs_A4_y=3 --jUL   $ADDflag
./fit --nevents=${NEVENTS} --tag=NEWA3_100M --run=full --post_tag=DEBUG${ADD} --degs_corr_x=10 --degs_corr_y=4 --degs_A0_x=3 --degs_A0_y=2 --degs_A1_x=3 --degs_A1_y=3 --degs_A2_x=3 --degs_A2_y=2 --degs_A3_x=3 --degs_A3_y=3 --degs_A4_x=3 --degs_A4_y=3 --jUL   $ADDflag


exit

echo "Doing CORR 3-2-2-2-2-2-2-2-2-2-2-2..."
./fit --nevents=${NEVENTS} --tag=SYMFIXCORR_10M --run=corr --post_tag=DEBUG${ADD} --degs_corr_x=3 --degs_corr_y=2 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=2 --jUL   $ADDflag
./fit --nevents=${NEVENTS} --tag=SYMFIXCORR_100M --run=corr --post_tag=DEBUG${ADD} --degs_corr_x=3 --degs_corr_y=2 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=2 --jUL   $ADDflag 
./fit --nevents=${NEVENTS} --tag=SYMFIXCORR_1G --run=corr --post_tag=DEBUG${ADD} --degs_corr_x=3 --degs_corr_y=2 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=2 --jUL   $ADDflag
./fit --nevents=${NEVENTS} --tag=SYMFIXCORR_4G --run=corr --post_tag=DEBUG${ADD} --degs_corr_x=3 --degs_corr_y=2 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=2 --jUL   $ADDflag
./fit --nevents=${NEVENTS} --tag=SYMFIXCORR_10G --run=corr --post_tag=DEBUG${ADD} --degs_corr_x=3 --degs_corr_y=2 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=2 --jUL   $ADDflag
./fit --nevents=${NEVENTS} --tag=SYMFIXCORR_40G --run=corr --post_tag=DEBUG${ADD} --degs_corr_x=3 --degs_corr_y=2 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=2 --jUL   $ADDflag

echo "Doing CORR 5-2-2-2-2-2-2-2-2-2-2-2..."
./fit --nevents=${NEVENTS} --tag=SYMFIXCORR_10M --run=corr --post_tag=DEBUG${ADD} --degs_corr_x=5 --degs_corr_y=2 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=2 --jUL   $ADDflag
./fit --nevents=${NEVENTS} --tag=SYMFIXCORR_100M --run=corr --post_tag=DEBUG${ADD} --degs_corr_x=5 --degs_corr_y=2 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=2 --jUL   $ADDflag 
./fit --nevents=${NEVENTS} --tag=SYMFIXCORR_1G --run=corr --post_tag=DEBUG${ADD} --degs_corr_x=5 --degs_corr_y=2 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=2 --jUL   $ADDflag
./fit --nevents=${NEVENTS} --tag=SYMFIXCORR_4G --run=corr --post_tag=DEBUG${ADD} --degs_corr_x=5 --degs_corr_y=2 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=2 --jUL   $ADDflag
./fit --nevents=${NEVENTS} --tag=SYMFIXCORR_10G --run=corr --post_tag=DEBUG${ADD} --degs_corr_x=5 --degs_corr_y=2 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=2 --jUL   $ADDflag
./fit --nevents=${NEVENTS} --tag=SYMFIXCORR_40G --run=corr --post_tag=DEBUG${ADD} --degs_corr_x=5 --degs_corr_y=2 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=2 --jUL   $ADDflag

echo "Doing CORR 4-2-2-2-2-2-2-2-2-2-2-2..."
./fit --nevents=${NEVENTS} --tag=SYMFIXCORR_10M --run=corr --post_tag=DEBUG${ADD} --degs_corr_x=4 --degs_corr_y=2 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=2 --jUL   $ADDflag
./fit --nevents=${NEVENTS} --tag=SYMFIXCORR_100M --run=corr --post_tag=DEBUG${ADD} --degs_corr_x=4 --degs_corr_y=2 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=2 --jUL   $ADDflag 
./fit --nevents=${NEVENTS} --tag=SYMFIXCORR_1G --run=corr --post_tag=DEBUG${ADD} --degs_corr_x=4 --degs_corr_y=2 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=2 --jUL   $ADDflag
./fit --nevents=${NEVENTS} --tag=SYMFIXCORR_4G --run=corr --post_tag=DEBUG${ADD} --degs_corr_x=4 --degs_corr_y=2 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=2 --jUL   $ADDflag
./fit --nevents=${NEVENTS} --tag=SYMFIXCORR_10G --run=corr --post_tag=DEBUG${ADD} --degs_corr_x=4 --degs_corr_y=2 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=2 --jUL   $ADDflag
./fit --nevents=${NEVENTS} --tag=SYMFIXCORR_40G --run=corr --post_tag=DEBUG${ADD} --degs_corr_x=4 --degs_corr_y=2 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=2 --jUL   $ADDflag

echo "Doing CORR 3-2-1-2-1-2-1-2-1-2-1-2..."
./fit --nevents=${NEVENTS} --tag=SYMFIXCORR_10M --run=corr --post_tag=DEBUG${ADD} --degs_corr_x=3 --degs_corr_y=2 --degs_A0_x=1 --degs_A0_y=2 --degs_A1_x=1 --degs_A1_y=2 --degs_A2_x=1 --degs_A2_y=2 --degs_A3_x=1 --degs_A3_y=2 --degs_A4_x=1 --degs_A4_y=2 --jUL  $ADDflag
./fit --nevents=${NEVENTS} --tag=SYMFIXCORR_100M --run=corr --post_tag=DEBUG${ADD} --degs_corr_x=3 --degs_corr_y=2 --degs_A0_x=1 --degs_A0_y=2 --degs_A1_x=1 --degs_A1_y=2 --degs_A2_x=1 --degs_A2_y=2 --degs_A3_x=1 --degs_A3_y=2 --degs_A4_x=1 --degs_A4_y=2 --jUL   $ADDflag 
./fit --nevents=${NEVENTS} --tag=SYMFIXCORR_1G --run=corr --post_tag=DEBUG${ADD} --degs_corr_x=3 --degs_corr_y=2 --degs_A0_x=1 --degs_A0_y=2 --degs_A1_x=1 --degs_A1_y=2 --degs_A2_x=1 --degs_A2_y=2 --degs_A3_x=1 --degs_A3_y=2 --degs_A4_x=1 --degs_A4_y=2 --jUL   $ADDflag
./fit --nevents=${NEVENTS} --tag=SYMFIXCORR_4G --run=corr --post_tag=DEBUG${ADD} --degs_corr_x=3 --degs_corr_y=2 --degs_A0_x=1 --degs_A0_y=2 --degs_A1_x=1 --degs_A1_y=2 --degs_A2_x=1 --degs_A2_y=2 --degs_A3_x=1 --degs_A3_y=2 --degs_A4_x=1 --degs_A4_y=2 --jUL   $ADDflag
./fit --nevents=${NEVENTS} --tag=SYMFIXCORR_10G --run=corr --post_tag=DEBUG${ADD} --degs_corr_x=3 --degs_corr_y=2 --degs_A0_x=1 --degs_A0_y=2 --degs_A1_x=1 --degs_A1_y=2 --degs_A2_x=1 --degs_A2_y=2 --degs_A3_x=1 --degs_A3_y=2 --degs_A4_x=1 --degs_A4_y=2 --jUL   $ADDflag
./fit --nevents=${NEVENTS} --tag=SYMFIXCORR_40G --run=corr --post_tag=DEBUG${ADD} --degs_corr_x=3 --degs_corr_y=2 --degs_A0_x=1 --degs_A0_y=2 --degs_A1_x=1 --degs_A1_y=2 --degs_A2_x=1 --degs_A2_y=2 --degs_A3_x=1 --degs_A3_y=2 --degs_A4_x=1 --degs_A4_y=2 --jUL   $ADDflag

echo "Doing CORR 3-2-1-1-1-1-1-1-1-1-1-1..."
./fit --nevents=${NEVENTS} --tag=SYMFIXCORR_10M --run=corr --post_tag=DEBUG${ADD} --degs_corr_x=3 --degs_corr_y=2 --degs_A0_x=1 --degs_A0_y=1 --degs_A1_x=1 --degs_A1_y=1 --degs_A2_x=1 --degs_A2_y=1 --degs_A3_x=1 --degs_A3_y=1 --degs_A4_x=1 --degs_A4_y=1 --jUL  $ADDflag
./fit --nevents=${NEVENTS} --tag=SYMFIXCORR_100M --run=corr --post_tag=DEBUG${ADD} --degs_corr_x=3 --degs_corr_y=2 --degs_A0_x=1 --degs_A0_y=1 --degs_A1_x=1 --degs_A1_y=1 --degs_A2_x=1 --degs_A2_y=1 --degs_A3_x=1 --degs_A3_y=1 --degs_A4_x=1 --degs_A4_y=1 --jUL   $ADDflag 
./fit --nevents=${NEVENTS} --tag=SYMFIXCORR_1G --run=corr --post_tag=DEBUG${ADD} --degs_corr_x=3 --degs_corr_y=2 --degs_A0_x=1 --degs_A0_y=1 --degs_A1_x=1 --degs_A1_y=1 --degs_A2_x=1 --degs_A2_y=1 --degs_A3_x=1 --degs_A3_y=1 --degs_A4_x=1 --degs_A4_y=1 --jUL   $ADDflag
./fit --nevents=${NEVENTS} --tag=SYMFIXCORR_4G --run=corr --post_tag=DEBUG${ADD} --degs_corr_x=3 --degs_corr_y=2 --degs_A0_x=1 --degs_A0_y=1 --degs_A1_x=1 --degs_A1_y=1 --degs_A2_x=1 --degs_A2_y=1 --degs_A3_x=1 --degs_A3_y=1 --degs_A4_x=1 --degs_A4_y=1 --jUL   $ADDflag
./fit --nevents=${NEVENTS} --tag=SYMFIXCORR_10G --run=corr --post_tag=DEBUG${ADD} --degs_corr_x=3 --degs_corr_y=2 --degs_A0_x=1 --degs_A0_y=1 --degs_A1_x=1 --degs_A1_y=1 --degs_A2_x=1 --degs_A2_y=1 --degs_A3_x=1 --degs_A3_y=1 --degs_A4_x=1 --degs_A4_y=1 --jUL   $ADDflag
./fit --nevents=${NEVENTS} --tag=SYMFIXCORR_40G --run=corr --post_tag=DEBUG${ADD} --degs_corr_x=3 --degs_corr_y=2 --degs_A0_x=1 --degs_A0_y=1 --degs_A1_x=1 --degs_A1_y=1 --degs_A2_x=1 --degs_A2_y=1 --degs_A3_x=1 --degs_A3_y=1 --degs_A4_x=1 --degs_A4_y=1 --jUL   $ADDflag

exit

echo "Doing GRID 8-6-1-1-1-1-1-1-1-1-1-1..."
./fit --nevents=${NEVENTS} --tag=SYMFIX_10M  --run=grid --post_tag=DEBUG${ADD} --degs_corr_x=8 --degs_corr_y=6 --degs_A0_x=1 --degs_A0_y=1 --degs_A1_x=1 --degs_A1_y=1 --degs_A2_x=1 --degs_A2_y=1 --degs_A3_x=1 --degs_A3_y=1 --degs_A4_x=1 --degs_A4_y=1 --jUL   $ADDflag 
./fit --nevents=${NEVENTS} --tag=SYMFIX_100M --run=grid --post_tag=DEBUG${ADD} --degs_corr_x=8 --degs_corr_y=6 --degs_A0_x=1 --degs_A0_y=1 --degs_A1_x=1 --degs_A1_y=1 --degs_A2_x=1 --degs_A2_y=1 --degs_A3_x=1 --degs_A3_y=1 --degs_A4_x=1 --degs_A4_y=1 --jUL   $ADDflag 
./fit --nevents=${NEVENTS} --tag=SYMFIX_1G   --run=grid --post_tag=DEBUG${ADD} --degs_corr_x=8 --degs_corr_y=6 --degs_A0_x=1 --degs_A0_y=1 --degs_A1_x=1 --degs_A1_y=1 --degs_A2_x=1 --degs_A2_y=1 --degs_A3_x=1 --degs_A3_y=1 --degs_A4_x=1 --degs_A4_y=1 --jUL   $ADDflag 
./fit --nevents=${NEVENTS} --tag=SYMFIX_4G   --run=grid --post_tag=DEBUG${ADD} --degs_corr_x=8 --degs_corr_y=6 --degs_A0_x=1 --degs_A0_y=1 --degs_A1_x=1 --degs_A1_y=1 --degs_A2_x=1 --degs_A2_y=1 --degs_A3_x=1 --degs_A3_y=1 --degs_A4_x=1 --degs_A4_y=1 --jUL   $ADDflag 
./fit --nevents=${NEVENTS} --tag=SYMFIX_10G  --run=grid --post_tag=DEBUG${ADD} --degs_corr_x=8 --degs_corr_y=6 --degs_A0_x=1 --degs_A0_y=1 --degs_A1_x=1 --degs_A1_y=1 --degs_A2_x=1 --degs_A2_y=1 --degs_A3_x=1 --degs_A3_y=1 --degs_A4_x=1 --degs_A4_y=1 --jUL   $ADDflag 

echo "Doing FULL 10-4-3-2-3-3-3-2-3-2-3-3..."
./fit --nevents=${NEVENTS} --tag=SYMFIX_10M --run=full --post_tag=DEBUG${ADD} --degs_corr_x=10 --degs_corr_y=4 --degs_A0_x=3 --degs_A0_y=2 --degs_A1_x=3 --degs_A1_y=3 --degs_A2_x=3 --degs_A2_y=2 --degs_A3_x=3 --degs_A3_y=2 --degs_A4_x=3 --degs_A4_y=3 --jUL   $ADDflag 
./fit --nevents=${NEVENTS} --tag=SYMFIX_100M --run=full --post_tag=DEBUG${ADD} --degs_corr_x=10 --degs_corr_y=4 --degs_A0_x=3 --degs_A0_y=2 --degs_A1_x=3 --degs_A1_y=3 --degs_A2_x=3 --degs_A2_y=2 --degs_A3_x=3 --degs_A3_y=2 --degs_A4_x=3 --degs_A4_y=3 --jUL   $ADDflag 
./fit --nevents=${NEVENTS} --tag=SYMFIX_1G --run=full --post_tag=DEBUG${ADD} --degs_corr_x=10 --degs_corr_y=4 --degs_A0_x=3 --degs_A0_y=2 --degs_A1_x=3 --degs_A1_y=3 --degs_A2_x=3 --degs_A2_y=2 --degs_A3_x=3 --degs_A3_y=2 --degs_A4_x=3 --degs_A4_y=3 --jUL   $ADDflag 
./fit --nevents=${NEVENTS} --tag=SYMFIX_4G --run=full --post_tag=DEBUG${ADD} --degs_corr_x=10 --degs_corr_y=4 --degs_A0_x=3 --degs_A0_y=2 --degs_A1_x=3 --degs_A1_y=3 --degs_A2_x=3 --degs_A2_y=2 --degs_A3_x=3 --degs_A3_y=2 --degs_A4_x=3 --degs_A4_y=3 --jUL   $ADDflag 
./fit --nevents=${NEVENTS} --tag=SYMFIX_10G --run=full --post_tag=DEBUG${ADD} --degs_corr_x=10 --degs_corr_y=4 --degs_A0_x=3 --degs_A0_y=2 --degs_A1_x=3 --degs_A1_y=3 --degs_A2_x=3 --degs_A2_y=2 --degs_A3_x=3 --degs_A3_y=2 --degs_A4_x=3 --degs_A4_y=3 --jUL   $ADDflag
./fit --nevents=${NEVENTS} --tag=SYMFIX_40G --run=full --post_tag=DEBUG${ADD} --degs_corr_x=10 --degs_corr_y=4 --degs_A0_x=3 --degs_A0_y=2 --degs_A1_x=3 --degs_A1_y=3 --degs_A2_x=3 --degs_A2_y=2 --degs_A3_x=3 --degs_A3_y=2 --degs_A4_x=3 --degs_A4_y=3 --jUL   $ADDflag 

echo "Doing FULL 10-6-3-2-3-3-3-2-3-2-3-3..."
./fit --nevents=${NEVENTS} --tag=SYMFIX_10M --run=full --post_tag=DEBUG${ADD} --degs_corr_x=10 --degs_corr_y=6 --degs_A0_x=3 --degs_A0_y=2 --degs_A1_x=3 --degs_A1_y=3 --degs_A2_x=3 --degs_A2_y=2 --degs_A3_x=3 --degs_A3_y=2 --degs_A4_x=3 --degs_A4_y=3 --jUL   $ADDflag 
./fit --nevents=${NEVENTS} --tag=SYMFIX_100M --run=full --post_tag=DEBUG${ADD} --degs_corr_x=10 --degs_corr_y=6 --degs_A0_x=3 --degs_A0_y=2 --degs_A1_x=3 --degs_A1_y=3 --degs_A2_x=3 --degs_A2_y=2 --degs_A3_x=3 --degs_A3_y=2 --degs_A4_x=3 --degs_A4_y=3 --jUL   $ADDflag 
./fit --nevents=${NEVENTS} --tag=SYMFIX_1G --run=full --post_tag=DEBUG${ADD} --degs_corr_x=10 --degs_corr_y=6 --degs_A0_x=3 --degs_A0_y=2 --degs_A1_x=3 --degs_A1_y=3 --degs_A2_x=3 --degs_A2_y=2 --degs_A3_x=3 --degs_A3_y=2 --degs_A4_x=3 --degs_A4_y=3 --jUL   $ADDflag 
./fit --nevents=${NEVENTS} --tag=SYMFIX_4G --run=full --post_tag=DEBUG${ADD} --degs_corr_x=10 --degs_corr_y=6 --degs_A0_x=3 --degs_A0_y=2 --degs_A1_x=3 --degs_A1_y=3 --degs_A2_x=3 --degs_A2_y=2 --degs_A3_x=3 --degs_A3_y=2 --degs_A4_x=3 --degs_A4_y=3 --jUL   $ADDflag 
./fit --nevents=${NEVENTS} --tag=SYMFIX_10G --run=full --post_tag=DEBUG${ADD} --degs_corr_x=10 --degs_corr_y=6 --degs_A0_x=3 --degs_A0_y=2 --degs_A1_x=3 --degs_A1_y=3 --degs_A2_x=3 --degs_A2_y=2 --degs_A3_x=3 --degs_A3_y=2 --degs_A4_x=3 --degs_A4_y=3 --jUL   $ADDflag 

echo "Doing FULL 12-4-3-2-3-3-3-2-3-2-3-3..."
./fit --nevents=${NEVENTS} --tag=SYMFIX_10M --run=full --post_tag=DEBUG${ADD} --degs_corr_x=12 --degs_corr_y=4 --degs_A0_x=3 --degs_A0_y=2 --degs_A1_x=3 --degs_A1_y=3 --degs_A2_x=3 --degs_A2_y=2 --degs_A3_x=3 --degs_A3_y=2 --degs_A4_x=3 --degs_A4_y=3 --jUL   $ADDflag 
./fit --nevents=${NEVENTS} --tag=SYMFIX_100M --run=full --post_tag=DEBUG${ADD} --degs_corr_x=12 --degs_corr_y=4 --degs_A0_x=3 --degs_A0_y=2 --degs_A1_x=3 --degs_A1_y=3 --degs_A2_x=3 --degs_A2_y=2 --degs_A3_x=3 --degs_A3_y=2 --degs_A4_x=3 --degs_A4_y=3 --jUL   $ADDflag 
./fit --nevents=${NEVENTS} --tag=SYMFIX_1G --run=full --post_tag=DEBUG${ADD} --degs_corr_x=12 --degs_corr_y=4 --degs_A0_x=3 --degs_A0_y=2 --degs_A1_x=3 --degs_A1_y=3 --degs_A2_x=3 --degs_A2_y=2 --degs_A3_x=3 --degs_A3_y=2 --degs_A4_x=3 --degs_A4_y=3 --jUL   $ADDflag 
./fit --nevents=${NEVENTS} --tag=SYMFIX_4G --run=full --post_tag=DEBUG${ADD} --degs_corr_x=12 --degs_corr_y=4 --degs_A0_x=3 --degs_A0_y=2 --degs_A1_x=3 --degs_A1_y=3 --degs_A2_x=3 --degs_A2_y=2 --degs_A3_x=3 --degs_A3_y=2 --degs_A4_x=3 --degs_A4_y=3 --jUL   $ADDflag 
./fit --nevents=${NEVENTS} --tag=SYMFIX_10G --run=full --post_tag=DEBUG${ADD} --degs_corr_x=12 --degs_corr_y=4 --degs_A0_x=3 --degs_A0_y=2 --degs_A1_x=3 --degs_A1_y=3 --degs_A2_x=3 --degs_A2_y=2 --degs_A3_x=3 --degs_A3_y=2 --degs_A4_x=3 --degs_A4_y=3 --jUL   $ADDflag 


exit

echo "Doing 3-2-1-1-..."
#./fit --nevents=${NEVENTS} --tag=SYMFIX_1M --run=corr --post_tag=DEBUG${ADD} --degs_corr_x=3 --degs_corr_y=2 --degs_A0_x=1 --degs_A0_y=1 --degs_A1_x=1 --degs_A1_y=1 --degs_A2_x=1 --degs_A2_y=1 --degs_A3_x=1 --degs_A3_y=1 --degs_A4_x=1 --degs_A4_y=1 --jUL   $ADDflag 
./fit --nevents=${NEVENTS} --tag=SYMFIX_10M --run=corr --post_tag=DEBUG${ADD} --degs_corr_x=3 --degs_corr_y=2 --degs_A0_x=1 --degs_A0_y=1 --degs_A1_x=1 --degs_A1_y=1 --degs_A2_x=1 --degs_A2_y=1 --degs_A3_x=1 --degs_A3_y=1 --degs_A4_x=1 --degs_A4_y=1 --jUL   $ADDflag 
./fit --nevents=${NEVENTS} --tag=SYMFIX_100M --run=corr --post_tag=DEBUG${ADD} --degs_corr_x=3 --degs_corr_y=2 --degs_A0_x=1 --degs_A0_y=1 --degs_A1_x=1 --degs_A1_y=1 --degs_A2_x=1 --degs_A2_y=1 --degs_A3_x=1 --degs_A3_y=1 --degs_A4_x=1 --degs_A4_y=1 --jUL   $ADDflag 
./fit --nevents=${NEVENTS} --tag=SYMFIX_1G --run=corr --post_tag=DEBUG${ADD} --degs_corr_x=3 --degs_corr_y=2 --degs_A0_x=1 --degs_A0_y=1 --degs_A1_x=1 --degs_A1_y=1 --degs_A2_x=1 --degs_A2_y=1 --degs_A3_x=1 --degs_A3_y=1 --degs_A4_x=1 --degs_A4_y=1 --jUL   $ADDflag 
./fit --nevents=${NEVENTS} --tag=SYMFIX_4G --run=corr --post_tag=DEBUG${ADD} --degs_corr_x=3 --degs_corr_y=2 --degs_A0_x=1 --degs_A0_y=1 --degs_A1_x=1 --degs_A1_y=1 --degs_A2_x=1 --degs_A2_y=1 --degs_A3_x=1 --degs_A3_y=1 --degs_A4_x=1 --degs_A4_y=1 --jUL   $ADDflag 
./fit --nevents=${NEVENTS} --tag=SYMFIX_10G --run=corr --post_tag=DEBUG${ADD} --degs_corr_x=3 --degs_corr_y=2 --degs_A0_x=1 --degs_A0_y=1 --degs_A1_x=1 --degs_A1_y=1 --degs_A2_x=1 --degs_A2_y=1 --degs_A3_x=1 --degs_A3_y=1 --degs_A4_x=1 --degs_A4_y=1 --jUL   $ADDflag 
./fit --nevents=${NEVENTS} --tag=SYMFIX_40G --run=corr --post_tag=DEBUG${ADD} --degs_corr_x=3 --degs_corr_y=2 --degs_A0_x=1 --degs_A0_y=1 --degs_A1_x=1 --degs_A1_y=1 --degs_A2_x=1 --degs_A2_y=1 --degs_A3_x=1 --degs_A3_y=1 --degs_A4_x=1 --degs_A4_y=1 --jUL   $ADDflag 
#./fit --nevents=${NEVENTS} --tag=SYMFIX_80G --run=corr --post_tag=DEBUG${ADD} --degs_corr_x=3 --degs_corr_y=2 --degs_A0_x=1 --degs_A0_y=1 --degs_A1_x=1 --degs_A1_y=1 --degs_A2_x=1 --degs_A2_y=1 --degs_A3_x=1 --degs_A3_y=1 --degs_A4_x=1 --degs_A4_y=1 --jUL   $ADDflag 

echo "Doing 3-4-2-2-..."
#./fit --nevents=${NEVENTS} --tag=SYMFIX_1M --run=corr --post_tag=DEBUG${ADD} --degs_corr_x=3 --degs_corr_y=4 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=3 --jUL   $ADDflag 
./fit --nevents=${NEVENTS} --tag=SYMFIX_10M --run=corr --post_tag=DEBUG${ADD} --degs_corr_x=3 --degs_corr_y=4 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=3 --jUL   $ADDflag 
./fit --nevents=${NEVENTS} --tag=SYMFIX_100M --run=corr --post_tag=DEBUG${ADD} --degs_corr_x=3 --degs_corr_y=4 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=3 --jUL   $ADDflag 
./fit --nevents=${NEVENTS} --tag=SYMFIX_1G --run=corr --post_tag=DEBUG${ADD} --degs_corr_x=3 --degs_corr_y=4 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=3 --jUL   $ADDflag 
./fit --nevents=${NEVENTS} --tag=SYMFIX_4G --run=corr --post_tag=DEBUG${ADD} --degs_corr_x=3 --degs_corr_y=4 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=3 --jUL   $ADDflag 
./fit --nevents=${NEVENTS} --tag=SYMFIX_10G --run=corr --post_tag=DEBUG${ADD} --degs_corr_x=3 --degs_corr_y=4 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=3 --jUL   $ADDflag 
./fit --nevents=${NEVENTS} --tag=SYMFIX_40G --run=corr --post_tag=DEBUG${ADD} --degs_corr_x=3 --degs_corr_y=4 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=3 --jUL   $ADDflag 
#./fit --nevents=${NEVENTS} --tag=SYMFIX_80G --run=corr --post_tag=DEBUG${ADD} --degs_corr_x=3 --degs_corr_y=4 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=3 --jUL   $ADDflag 
#./fit --nevents=${NEVENTS} --tag=SYMFIX_200G --run=corr --post_tag=DEBUG${ADD} --degs_corr_x=3 --degs_corr_y=4 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=3 --jUL   $ADDflag 

echo "Doing 8-6-1-1-..."
#./fit --nevents=${NEVENTS} --tag=SYMFIX_1M --run=grid --post_tag=DEBUG${ADD} --degs_corr_x=8 --degs_corr_y=6 --degs_A0_x=1 --degs_A0_y=1 --degs_A1_x=1 --degs_A1_y=1 --degs_A2_x=1 --degs_A2_y=1 --degs_A3_x=1 --degs_A3_y=1 --degs_A4_x=1 --degs_A4_y=1 --jUL   $ADDflag 
./fit --nevents=${NEVENTS} --tag=SYMFIX_10M --run=grid --post_tag=DEBUG${ADD} --degs_corr_x=8 --degs_corr_y=6 --degs_A0_x=1 --degs_A0_y=1 --degs_A1_x=1 --degs_A1_y=1 --degs_A2_x=1 --degs_A2_y=1 --degs_A3_x=1 --degs_A3_y=1 --degs_A4_x=1 --degs_A4_y=1 --jUL   $ADDflag 
./fit --nevents=${NEVENTS} --tag=SYMFIX_100M --run=grid --post_tag=DEBUG${ADD} --degs_corr_x=8 --degs_corr_y=6 --degs_A0_x=1 --degs_A0_y=1 --degs_A1_x=1 --degs_A1_y=1 --degs_A2_x=1 --degs_A2_y=1 --degs_A3_x=1 --degs_A3_y=1 --degs_A4_x=1 --degs_A4_y=1 --jUL   $ADDflag 
./fit --nevents=${NEVENTS} --tag=SYMFIX_1G --run=grid --post_tag=DEBUG${ADD} --degs_corr_x=8 --degs_corr_y=6 --degs_A0_x=1 --degs_A0_y=1 --degs_A1_x=1 --degs_A1_y=1 --degs_A2_x=1 --degs_A2_y=1 --degs_A3_x=1 --degs_A3_y=1 --degs_A4_x=1 --degs_A4_y=1 --jUL   $ADDflag 
./fit --nevents=${NEVENTS} --tag=SYMFIX_4G --run=grid --post_tag=DEBUG${ADD} --degs_corr_x=8 --degs_corr_y=6 --degs_A0_x=1 --degs_A0_y=1 --degs_A1_x=1 --degs_A1_y=1 --degs_A2_x=1 --degs_A2_y=1 --degs_A3_x=1 --degs_A3_y=1 --degs_A4_x=1 --degs_A4_y=1 --jUL   $ADDflag 
./fit --nevents=${NEVENTS} --tag=SYMFIX_10G --run=grid --post_tag=DEBUG${ADD} --degs_corr_x=8 --degs_corr_y=6 --degs_A0_x=1 --degs_A0_y=1 --degs_A1_x=1 --degs_A1_y=1 --degs_A2_x=1 --degs_A2_y=1 --degs_A3_x=1 --degs_A3_y=1 --degs_A4_x=1 --degs_A4_y=1 --jUL   $ADDflag 
#./fit --nevents=${NEVENTS} --tag=SYMFIX_40G --run=grid --post_tag=DEBUG${ADD} --degs_corr_x=8 --degs_corr_y=6 --degs_A0_x=1 --degs_A0_y=1 --degs_A1_x=1 --degs_A1_y=1 --degs_A2_x=1 --degs_A2_y=1 --degs_A3_x=1 --degs_A3_y=1 --degs_A4_x=1 --degs_A4_y=1 --jUL   $ADDflag 
#./fit --nevents=${NEVENTS} --tag=SYMFIX_80G --run=grid --post_tag=DEBUG${ADD} --degs_corr_x=8 --degs_corr_y=6 --degs_A0_x=1 --degs_A0_y=1 --degs_A1_x=1 --degs_A1_y=1 --degs_A2_x=1 --degs_A2_y=1 --degs_A3_x=1 --degs_A3_y=1 --degs_A4_x=1 --degs_A4_y=1 --jUL   $ADDflag 


exit

./fit --nevents=100000000 --tag=DEV_1M --run=closure --post_tag=DEBUG_ADDMC_ULA0A1A2A3 --degs_corr_x=3 --degs_corr_y=4 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=3 --jUL  --do_cheb_as_modifiers --add_MC_uncert --j0 --j1 --j2 --j3
./fit --nevents=100000000 --tag=DEV_10M --run=closure --post_tag=DEBUG_ADDMC_ULA0A1A2A3 --degs_corr_x=3 --degs_corr_y=4 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=3 --jUL  --do_cheb_as_modifiers --add_MC_uncert --j0 --j1 --j2 --j3
./fit --nevents=100000000 --tag=DEV_100M --run=closure --post_tag=DEBUG_ADDMC_ULA0A1A2A3 --degs_corr_x=3 --degs_corr_y=4 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=3 --jUL  --do_cheb_as_modifiers --add_MC_uncert --j0 --j1 --j2 --j3
./fit --nevents=100000000 --tag=DEV_1G --run=closure --post_tag=DEBUG_ADDMC_ULA0A1A2A3 --degs_corr_x=3 --degs_corr_y=4 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=3 --jUL  --do_cheb_as_modifiers --add_MC_uncert --j0 --j1 --j2 --j3
./fit --nevents=100000000 --tag=DEV_10G --run=closure --post_tag=DEBUG_ADDMC_ULA0A1A2A3 --degs_corr_x=3 --degs_corr_y=4 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=3 --jUL  --do_cheb_as_modifiers --add_MC_uncert --j0 --j1 --j2 --j3
./fit --nevents=100000000 --tag=DEV_40G --run=closure --post_tag=DEBUG_ADDMC_ULA0A1A2A3 --degs_corr_x=3 --degs_corr_y=4 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=3 --jUL  --do_cheb_as_modifiers --add_MC_uncert --j0 --j1 --j2 --j3
./fit --nevents=100000000 --tag=DEV_80G --run=closure --post_tag=DEBUG_ADDMC_ULA0A1A2A3 --degs_corr_x=3 --degs_corr_y=4 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=3 --jUL  --do_cheb_as_modifiers --add_MC_uncert --j0 --j1 --j2 --j3

./fit --nevents=100000000 --tag=DEV_1M --run=closure --post_tag=DEBUG_ADDMC_ULA0A1A2 --degs_corr_x=3 --degs_corr_y=4 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=3 --jUL  --do_cheb_as_modifiers --add_MC_uncert --j0 --j1 --j2
./fit --nevents=100000000 --tag=DEV_10M --run=closure --post_tag=DEBUG_ADDMC_ULA0A1A2 --degs_corr_x=3 --degs_corr_y=4 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=3 --jUL  --do_cheb_as_modifiers --add_MC_uncert --j0 --j1 --j2
./fit --nevents=100000000 --tag=DEV_100M --run=closure --post_tag=DEBUG_ADDMC_ULA0A1A2 --degs_corr_x=3 --degs_corr_y=4 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=3 --jUL  --do_cheb_as_modifiers --add_MC_uncert --j0 --j1 --j2
./fit --nevents=100000000 --tag=DEV_1G --run=closure --post_tag=DEBUG_ADDMC_ULA0A1A2 --degs_corr_x=3 --degs_corr_y=4 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=3 --jUL  --do_cheb_as_modifiers --add_MC_uncert --j0 --j1 --j2
./fit --nevents=100000000 --tag=DEV_10G --run=closure --post_tag=DEBUG_ADDMC_ULA0A1A2 --degs_corr_x=3 --degs_corr_y=4 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=3 --jUL  --do_cheb_as_modifiers --add_MC_uncert --j0 --j1 --j2
./fit --nevents=100000000 --tag=DEV_40G --run=closure --post_tag=DEBUG_ADDMC_ULA0A1A2 --degs_corr_x=3 --degs_corr_y=4 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=3 --jUL  --do_cheb_as_modifiers --add_MC_uncert --j0 --j1 --j2
./fit --nevents=100000000 --tag=DEV_80G --run=closure --post_tag=DEBUG_ADDMC_ULA0A1A2 --degs_corr_x=3 --degs_corr_y=4 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=3 --jUL  --do_cheb_as_modifiers --add_MC_uncert --j0 --j1 --j2

./fit --nevents=100000000 --tag=DEV_1M --run=closure --post_tag=DEBUG_ADDMC_ULA0A1 --degs_corr_x=3 --degs_corr_y=4 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=3 --jUL  --do_cheb_as_modifiers --add_MC_uncert --j0 --j1 
./fit --nevents=100000000 --tag=DEV_10M --run=closure --post_tag=DEBUG_ADDMC_ULA0A1 --degs_corr_x=3 --degs_corr_y=4 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=3 --jUL  --do_cheb_as_modifiers --add_MC_uncert --j0 --j1 
./fit --nevents=100000000 --tag=DEV_100M --run=closure --post_tag=DEBUG_ADDMC_ULA0A1 --degs_corr_x=3 --degs_corr_y=4 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=3 --jUL  --do_cheb_as_modifiers --add_MC_uncert --j0 --j1 
./fit --nevents=100000000 --tag=DEV_1G --run=closure --post_tag=DEBUG_ADDMC_ULA0A1 --degs_corr_x=3 --degs_corr_y=4 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=3 --jUL  --do_cheb_as_modifiers --add_MC_uncert --j0 --j1 
./fit --nevents=100000000 --tag=DEV_10G --run=closure --post_tag=DEBUG_ADDMC_ULA0A1 --degs_corr_x=3 --degs_corr_y=4 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=3 --jUL  --do_cheb_as_modifiers --add_MC_uncert --j0 --j1 
./fit --nevents=100000000 --tag=DEV_40G --run=closure --post_tag=DEBUG_ADDMC_ULA0A1 --degs_corr_x=3 --degs_corr_y=4 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=3 --jUL  --do_cheb_as_modifiers --add_MC_uncert --j0 --j1
./fit --nevents=100000000 --tag=DEV_80G --run=closure --post_tag=DEBUG_ADDMC_ULA0A1 --degs_corr_x=3 --degs_corr_y=4 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=3 --jUL  --do_cheb_as_modifiers --add_MC_uncert --j0 --j1 

./fit --nevents=100000000 --tag=DEV_1M --run=closure --post_tag=DEBUG_ADDMC_ULA0 --degs_corr_x=3 --degs_corr_y=4 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=3 --jUL  --do_cheb_as_modifiers --add_MC_uncert --j0 
./fit --nevents=100000000 --tag=DEV_10M --run=closure --post_tag=DEBUG_ADDMC_ULA0 --degs_corr_x=3 --degs_corr_y=4 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=3 --jUL  --do_cheb_as_modifiers --add_MC_uncert --j0  
./fit --nevents=100000000 --tag=DEV_100M --run=closure --post_tag=DEBUG_ADDMC_ULA0 --degs_corr_x=3 --degs_corr_y=4 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=3 --jUL  --do_cheb_as_modifiers --add_MC_uncert --j0  
./fit --nevents=100000000 --tag=DEV_1G --run=closure --post_tag=DEBUG_ADDMC_ULA0 --degs_corr_x=3 --degs_corr_y=4 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=3 --jUL  --do_cheb_as_modifiers --add_MC_uncert --j0  
./fit --nevents=100000000 --tag=DEV_10G --run=closure --post_tag=DEBUG_ADDMC_ULA0 --degs_corr_x=3 --degs_corr_y=4 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=3 --jUL  --do_cheb_as_modifiers --add_MC_uncert --j0  
./fit --nevents=100000000 --tag=DEV_40G --run=closure --post_tag=DEBUG_ADDMC_ULA0 --degs_corr_x=3 --degs_corr_y=4 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=3 --jUL  --do_cheb_as_modifiers --add_MC_uncert --j0
./fit --nevents=100000000 --tag=DEV_80G --run=closure --post_tag=DEBUG_ADDMC_ULA0 --degs_corr_x=3 --degs_corr_y=4 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=3 --jUL  --do_cheb_as_modifiers --add_MC_uncert --j0  

./fit --nevents=100000000 --tag=DEV_1M --run=closure --post_tag=DEBUG_ADDMC_UL --degs_corr_x=3 --degs_corr_y=4 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=3 --jUL  --do_cheb_as_modifiers --add_MC_uncert  
./fit --nevents=100000000 --tag=DEV_10M --run=closure --post_tag=DEBUG_ADDMC_UL --degs_corr_x=3 --degs_corr_y=4 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=3 --jUL  --do_cheb_as_modifiers --add_MC_uncert  
./fit --nevents=100000000 --tag=DEV_100M --run=closure --post_tag=DEBUG_ADDMC_UL --degs_corr_x=3 --degs_corr_y=4 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=3 --jUL  --do_cheb_as_modifiers --add_MC_uncert  
./fit --nevents=100000000 --tag=DEV_1G --run=closure --post_tag=DEBUG_ADDMC_UL --degs_corr_x=3 --degs_corr_y=4 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=3 --jUL  --do_cheb_as_modifiers --add_MC_uncert  
./fit --nevents=100000000 --tag=DEV_10G --run=closure --post_tag=DEBUG_ADDMC_UL --degs_corr_x=3 --degs_corr_y=4 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=3 --jUL  --do_cheb_as_modifiers --add_MC_uncert  
./fit --nevents=100000000 --tag=DEV_40G --run=closure --post_tag=DEBUG_ADDMC_UL --degs_corr_x=3 --degs_corr_y=4 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=3 --jUL  --do_cheb_as_modifiers --add_MC_uncert
./fit --nevents=100000000 --tag=DEV_80G --run=closure --post_tag=DEBUG_ADDMC_UL --degs_corr_x=3 --degs_corr_y=4 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=2 --degs_A1_y=2 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=2 --degs_A3_y=2 --degs_A4_x=2 --degs_A4_y=3 --jUL  --do_cheb_as_modifiers --add_MC_uncert  
