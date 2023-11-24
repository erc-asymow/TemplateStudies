#!/bin/sh

#old tag = NEWA3ZEROSMEARY3p5X0p5Ymax70GW1p0NEWTOYX



#python run.py --tag=NEWA3ZEROSMEARGW1p0Y3p5 --algo=jac2_systs_vsN --smear

#python run.py --tag=REALISTICWPX0p3Y3p0 --algo=jac2_systs_vsN #--smear
#python run.py --tag=CMP --algo=jac2_systs_vsN #--smear
#python run.py --tag=FULLPS2WPX0p4Y3p5 --algo=jac2_systs_vsN #--smear
#python run.py --tag=FULLPSWPX0p4Y3p5 --algo=jac2_systs_vsN #--smear
#python run.py --tag=REALISTICWPX0p3Y3p0 --algo=fit_systs_vsN --post_tag=DEBUG --nevents=100000000
#python run.py --tag=DEVRNDWPX0p3Y3p0 --algo=fit_systs_vsN --post_tag=DEBUG --nevents=100000000
#python run.py --tag=FULLPS2WPX0p4Y3p5 --algo=fit_systs_vsN --post_tag=DEBUG --nevents=100000000
python run.py --tag=CMP --algo=fit_systs_vsN --post_tag=DEBUG --nevents=100000000

#python run.py --tag=NEWA3ZEROSMEARGW1p0 --algo=fit_systs_vsN --post_tag=DEBUG --nevents=100000000
#python run.py --tag=NEWA3ZEROSMEARGW1p0 --algo=fit_systs_vsN --post_tag=DEBUGrebin23 --rebinX=2 --rebinY=3 --nevents=100000000

#python run.py --tag=NEWA3ZEROSMEARGW1p0NEWTOYY_4G --algo=fit_grid_fast --post_tag=DEBUG --nevents=100000000

#./jac2 --nevents=1000000000  --tag=NEWA3_1G --run=full --degs_corr_x=10 --degs_corr_y=4 --degs_A0_x=3 --degs_A0_y=2 --degs_A1_x=3 --degs_A1_y=3 --degs_A2_x=3 --degs_A2_y=2 --degs_A3_x=3 --degs_A3_y=3 --degs_A4_x=3 --degs_A4_y=3  --max_y=3.5
#./jac2 --nevents=1000000000  --tag=NEWA3_1G --run=full --degs_corr_x=10 --degs_corr_y=4 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=3 --degs_A1_y=3 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=3 --degs_A3_y=3 --degs_A4_x=3 --degs_A4_y=3  --max_y=3.5
#./jac2 --nevents=1000000000  --tag=NEWA3_1G --run=full --degs_corr_x=10 --degs_corr_y=4 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=3 --degs_A1_y=3 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=3 --degs_A3_y=3 --degs_A4_x=1 --degs_A4_y=3  --max_y=3.5
#./jac2 --nevents=1000000000  --tag=NEWA3_1G --run=full --degs_corr_x=10 --degs_corr_y=4 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=3 --degs_A1_y=3 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=3 --degs_A3_y=3 --degs_A4_x=2 --degs_A4_y=3  --max_y=3.5
#./jac2 --nevents=100000000  --tag=NEWA3_100M --run=full --degs_corr_x=10 --degs_corr_y=4 --degs_A0_x=2 --degs_A0_y=2 --degs_A1_x=3 --degs_A1_y=3 --degs_A2_x=2 --degs_A2_y=2 --degs_A3_x=3 --degs_A3_y=3 --degs_A4_x=1 --degs_A4_y=3  --max_y=3.5
#./jac2 --nevents=100000000  --tag=DEBUGNEWA3_100M --run=grid --degs_corr_x=8 --degs_corr_y=6 --degs_A0_x=1 --degs_A0_y=1 --degs_A1_x=1 --degs_A1_y=1 --degs_A2_x=1 --degs_A2_y=1 --degs_A3_x=1 --degs_A3_y=1 --degs_A4_x=1 --degs_A4_y=1  --max_y=2.5

