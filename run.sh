#!/bin/sh

#python run.py --algo=fit_grid --tag=addmass0 --post_tag=1M  --nevents=100000000
#python run.py --algo=fit_grid --tag=addmass1 --post_tag=2M  --nevents=100000000
#python run.py --algo=fit_grid --tag=addmass2 --post_tag=4M  --nevents=100000000
#python run.py --algo=fit_grid --tag=addmass3 --post_tag=8M  --nevents=100000000
#python run.py --algo=fit_grid --tag=addmass4 --post_tag=16M --nevents=100000000
#python run.py --algo=fit_grid --tag=addmass5 --post_tag=32M --nevents=100000000
#python run.py --algo=fit_grid --tag=addmass6 --post_tag=64M --nevents=100000000

#./fit_grid --nevents=100000000 --tag=addmass0 --run=grid --post_tag=SCALE_ALL_BUT_A4 --degs_corr_x=8 --degs_corr_y=4 --jUL --j0 --j1 --j2 --j3 --j4 --scale0 --scale1 --scale2 --scale3
#./fit_grid --nevents=100000000 --tag=addmass1 --run=grid --post_tag=SCALE_ALL_BUT_A4 --degs_corr_x=8 --degs_corr_y=4 --jUL --j0 --j1 --j2 --j3 --j4 --scale0 --scale1 --scale2 --scale3
#./fit_grid --nevents=100000000 --tag=addmass2 --run=grid --post_tag=SCALE_ALL_BUT_A4 --degs_corr_x=8 --degs_corr_y=4 --jUL --j0 --j1 --j2 --j3 --j4 --scale0 --scale1 --scale2 --scale3
#./fit_grid --nevents=100000000 --tag=addmass3 --run=grid --post_tag=SCALE_ALL_BUT_A4 --degs_corr_x=8 --degs_corr_y=4 --jUL --j0 --j1 --j2 --j3 --j4 --scale0 --scale1 --scale2 --scale3
#./fit_grid --nevents=100000000 --tag=addmass4 --run=grid --post_tag=SCALE_ALL_BUT_A4 --degs_corr_x=8 --degs_corr_y=4 --jUL --j0 --j1 --j2 --j3 --j4 --scale0 --scale1 --scale2 --scale3
#./fit_grid --nevents=100000000 --tag=addmass5 --run=grid --post_tag=SCALE_ALL_BUT_A4 --degs_corr_x=8 --degs_corr_y=4 --jUL --j0 --j1 --j2 --j3 --j4 --scale0 --scale1 --scale2 --scale3
#./fit_grid --nevents=100000000 --tag=addmass6 --run=grid --post_tag=SCALE_ALL_BUT_A4 --degs_corr_x=8 --degs_corr_y=4 --jUL --j0 --j1 --j2 --j3 --j4 --scale0 --scale1 --scale2 --scale3
#./fit_grid --nevents=100000000 --tag=addmass7 --run=grid --post_tag=SCALE_ALL_BUT_A4 --degs_corr_x=8 --degs_corr_y=4 --jUL --j0 --j1 --j2 --j3 --j4 --scale0 --scale1 --scale2 --scale3
#./fit_grid --nevents=100000000 --tag=addmass8 --run=grid --post_tag=SCALE_ALL_BUT_A4 --degs_corr_x=8 --degs_corr_y=4 --jUL --j0 --j1 --j2 --j3 --j4 --scale0 --scale1 --scale2 --scale3


#./fit --nevents=100000000 --tag=addmass2 --run=closure --post_tag=1G_SCALE_ALL_BUT_A4 --degs_corr_x=10 --degs_corr_y=4 --degs_A0_x=3 --degs_A0_y=3 --degs_A1_x=3 --degs_A1_y=3 --degs_A2_x=3 --degs_A2_y=3 --degs_A3_x=3 --degs_A3_y=3 --degs_A4_x=3 --degs_A4_y=3 --jUL --j0 --j1 --j2 --j3 --j4 --scale0 --scale1 --scale2 --scale3 --verbose
#./fit --nevents=100000000 --tag=addmass3 --run=closure --post_tag=4G_SCALE_ALL_BUT_A4 --degs_corr_x=10 --degs_corr_y=4 --degs_A0_x=3 --degs_A0_y=3 --degs_A1_x=3 --degs_A1_y=3 --degs_A2_x=3 --degs_A2_y=3 --degs_A3_x=3 --degs_A3_y=3 --degs_A4_x=3 --degs_A4_y=3 --jUL --j0 --j1 --j2 --j3 --j4 --scale0 --scale1 --scale2 --scale3 --verbose
#./fit --nevents=100000000 --tag=addmass4 --run=closure --post_tag=10G_SCALE_ALL_BUT_A4_rebin22 --degs_corr_x=10 --degs_corr_y=4 --degs_A0_x=3 --degs_A0_y=3 --degs_A1_x=3 --degs_A1_y=3 --degs_A2_x=3 --degs_A2_y=3 --degs_A3_x=3 --degs_A3_y=3 --degs_A4_x=3 --degs_A4_y=3 --jUL --j0 --j1 --j2 --j3 --j4 --scale0 --scale1 --scale2 --scale3 --verbose --rebinX=2 --rebinY=2

#python run.py --algo=fit --tag=addmass2 --post_tag=1G_SCALE_ALL_BUT_A4_FIXSCALE  --nevents=100000000 --scale0 --scale1 --scale2 --scale3
#python run.py --algo=fit --tag=addmass3 --post_tag=4G_SCALE_ALL_BUT_A4_FIXSCALE  --nevents=100000000 --scale0 --scale1 --scale2 --scale3
#python run.py --algo=fit --tag=addmass4 --post_tag=10G_SCALE_ALL_BUT_A4_FIXSCALE --nevents=100000000 --scale0 --scale1 --scale2 --scale3
python run.py --algo=fit --tag=addmass4 --post_tag=DEBUG --nevents=100000000 --scale0 --scale1 --scale2 --scale3

#python run.py --algo=fit --tag=addmass2 --post_tag=1G_SCALE_ALL_BUT_A0A4_FIXSCALE  --nevents=100000000  --scale1 --scale2 --scale3
#python run.py --algo=fit --tag=addmass3 --post_tag=4G_SCALE_ALL_BUT_A0A4_FIXSCALE  --nevents=100000000  --scale1 --scale2 --scale3
#python run.py --algo=fit --tag=addmass4 --post_tag=10G_SCALE_ALL_BUT_A0A4_FIXSCALE --nevents=100000000  --scale1 --scale2 --scale3

#python run.py --algo=fit --tag=addmass2 --post_tag=1G_SCALE_ALL_BUT_A0A1A4_FIXSCALE  --nevents=100000000  --scale2 --scale3
#python run.py --algo=fit --tag=addmass3 --post_tag=4G_SCALE_ALL_BUT_A0A1A4_FIXSCALE  --nevents=100000000  --scale2 --scale3
#python run.py --algo=fit --tag=addmass4 --post_tag=10G_SCALE_ALL_BUT_A0A1A4_FIXSCALE --nevents=100000000  --scale2 --scale3

#python run.py --algo=fit --tag=addmass2 --post_tag=1G_SCALE_ALL_BUT_A0A1A2A4_FIXSCALE  --nevents=100000000  --scale3
#python run.py --algo=fit --tag=addmass3 --post_tag=4G_SCALE_ALL_BUT_A0A1A2A4_FIXSCALE  --nevents=100000000  --scale3
#python run.py --algo=fit --tag=addmass4 --post_tag=10G_SCALE_ALL_BUT_A0A1A2A4_FIXSCALE --nevents=100000000  --scale3

#python run.py --algo=fit --tag=addmass2 --post_tag=1G_SCALE_ALL_BUT_A0A1A2A3A4_FIXSCALE  --nevents=100000000
#python run.py --algo=fit --tag=addmass3 --post_tag=4G_SCALE_ALL_BUT_A0A1A2A3A4_FIXSCALE  --nevents=100000000
#python run.py --algo=fit --tag=addmass4 --post_tag=10G_SCALE_ALL_BUT_A0A1A2A3A4_FIXSCALE --nevents=100000000
