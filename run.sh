#!/bin/sh

#./fit_grid --nevents=100000000 --tag=addmass7_seed0   --run=grid --post_tag=seed0 --degs_corr_x=10 --degs_corr_y=4 --jUL --j0 --j1 --j2 --j3 --j4
#./fit_grid --nevents=100000000 --tag=addmass7_merged  --run=grid --post_tag=merged --degs_corr_x=10 --degs_corr_y=4 --jUL --j0 --j1 --j2 --j3 --j4
#./fit_grid --nevents=100000000 --tag=addmass7_merged2 --run=grid --post_tag=merged2 --degs_corr_x=10 --degs_corr_y=4 --jUL --j0 --j1 --j2 --j3 --j4
#./fit_grid --nevents=100000000 --tag=addmass7_merged3 --run=grid --post_tag=merged3 --degs_corr_x=10 --degs_corr_y=4 --jUL --j0 --j1 --j2 --j3 --j4
#./fit_grid --nevents=100000000 --tag=addmass7_merged4 --run=grid --post_tag=merged4 --degs_corr_x=10 --degs_corr_y=4 --jUL --j0 --j1 --j2 --j3 --j4
#./fit_grid --nevents=100000000 --tag=addmass7_merged5 --run=grid --post_tag=merged5 --degs_corr_x=10 --degs_corr_y=4 --jUL --j0 --j1 --j2 --j3 --j4

#./fit_grid --nevents=100000000 --tag=addmass7_hadd --run=grid --post_tag=hadd --degs_corr_x=10 --degs_corr_y=4 --jUL --j0 --j1 --j2 --j3 --j4

#./fit_grid --nevents=100000000 --tag=TEST_10M --run=exact --post_tag=TEST_10M --degs_corr_x=10 --degs_corr_y=4 --jUL --j0 --j1 --j2 --j3 --j4
#./fit_grid --nevents=100000000 --tag=TEST_100M --run=exact --post_tag=TEST_100M --degs_corr_x=10 --degs_corr_y=4 --jUL --j0 --j1 --j2 --j3 --j4
#./fit_grid --nevents=100000000 --tag=TEST_1G --run=exact --post_tag=TEST_1G --degs_corr_x=10 --degs_corr_y=4 --jUL --j0 --j1 --j2 --j3 --j4

#./fit_grid --nevents=100000000 --tag=addmass7  --run=grid --post_tag=debug128M  --degs_corr_x=8 --degs_corr_y=4 --jUL --j0 --j1 --j2 --j3 --j4 --rebinX=2 --rebinY=2
#./fit_grid --nevents=100000000 --tag=addmass8  --run=grid --post_tag=debug256M  --degs_corr_x=8 --degs_corr_y=4 --jUL --j0 --j1 --j2 --j3 --j4 --rebinX=2 --rebinY=2
#./fit_grid --nevents=100000000 --tag=addmass11 --run=grid --post_tag=debug2048M --degs_corr_x=8 --degs_corr_y=4 --jUL --j0 --j1 --j2 --j3 --j4 --rebinX=2 --rebinY=2

#python run.py --algo=jac3 --tag=addmass0  --nevents=1000000
#python run.py --algo=jac3 --tag=addmass1  --nevents=2000000
#python run.py --algo=jac3 --tag=addmass2  --nevents=4000000
#python run.py --algo=jac3 --tag=addmass3  --nevents=8000000
#python run.py --algo=jac3 --tag=addmass4  --nevents=16000000
#python run.py --algo=jac3 --tag=addmass5  --nevents=32000000
#python run.py --algo=jac3 --tag=addmass6  --nevents=64000000
#python run.py --algo=jac3 --tag=addmass7  --nevents=128000000
#python run.py --algo=jac3 --tag=addmass8  --nevents=256000000
#python run.py --algo=jac3 --tag=addmass9  --nevents=512000000
#python run.py --algo=jac3 --tag=addmass10 --nevents=1024000000
#python run.py --algo=jac3 --tag=addmass11 --nevents=2048000000

#./jac3 --nevents=128000000 --tag=addmass7_seed0 --run=grid --degs_corr_x=10 --degs_corr_y=4 --toyTF2_corr --seed=4357
#./jac3 --nevents=128000000 --tag=addmass7_seed1 --run=grid --degs_corr_x=10 --degs_corr_y=4 --toyTF2_corr --seed=5835
#./jac3 --nevents=128000000 --tag=addmass7_seed2 --run=grid --degs_corr_x=10 --degs_corr_y=4 --toyTF2_corr --seed=1097
#./jac3 --nevents=128000000 --tag=addmass7_seed3 --run=grid --degs_corr_x=10 --degs_corr_y=4 --toyTF2_corr --seed=6724
#./jac3 --nevents=128000000 --tag=addmass7_seed4 --run=grid --degs_corr_x=10 --degs_corr_y=4 --toyTF2_corr --seed=9061
#./jac3 --nevents=128000000 --tag=addmass7_seed5 --run=grid --degs_corr_x=10 --degs_corr_y=4 --toyTF2_corr --seed=3487

#./fit --nevents=100000000 --tag=jacsvsmass100M_50shift --run=closure --post_tag=TEST --degs_corr_x=10 --degs_corr_y=4 --degs_A0_x=3 --degs_A0_y=3 --degs_A1_x=3 --degs_A1_y=3 --degs_A2_x=3 --degs_A2_y=3 --degs_A3_x=3 --degs_A3_y=3 --degs_A4_x=3 --degs_A4_y=3 --jUL --j0 --j1 --j2 --j3 --j4

#./fit --nevents=100000000 --tag=jacsvsmass1G_50shift --run=closure --post_tag=TEST --degs_corr_x=10 --degs_corr_y=4 --degs_A0_x=3 --degs_A0_y=3 --degs_A1_x=3 --degs_A1_y=3 --degs_A2_x=3 --degs_A2_y=3 --degs_A3_x=3 --degs_A3_y=3 --degs_A4_x=3 --degs_A4_y=3 --jUL --j0 --j1 --j2 --j3 --j4

#./fit --nevents=100000000 --tag=jacVsM_4G --run=closure --post_tag=TEST --degs_corr_x=10 --degs_corr_y=4 --degs_A0_x=3 --degs_A0_y=3 --degs_A1_x=3 --degs_A1_y=3 --degs_A2_x=3 --degs_A2_y=3 --degs_A3_x=3 --degs_A3_y=3 --degs_A4_x=3 --degs_A4_y=3 --jUL --j0 --j1 --j2 --j3 --j4 --verbose

#python run.py --algo=jac2_systs --tag=SYSTS_1G --nevents=1000000000

#./fit --nevents=100000000 --tag=jacsvsmass100M_50shift --run=closure --post_tag=idx0 --degs_corr_x=10 --degs_corr_y=4 --degs_A0_x=3 --degs_A0_y=3 --degs_A1_x=3 --degs_A1_y=3 --degs_A2_x=3 --degs_A2_y=3 --degs_A3_x=3 --degs_A3_y=3 --degs_A4_x=3 --degs_A4_y=3 --jUL --j0 --j1 --j2 --j3 --j4 --jacmass=0
#./fit --nevents=100000000 --tag=jacsvsmass100M_50shift --run=closure --post_tag=idx1 --degs_corr_x=10 --degs_corr_y=4 --degs_A0_x=3 --degs_A0_y=3 --degs_A1_x=3 --degs_A1_y=3 --degs_A2_x=3 --degs_A2_y=3 --degs_A3_x=3 --degs_A3_y=3 --degs_A4_x=3 --degs_A4_y=3 --jUL --j0 --j1 --j2 --j3 --j4 --jacmass=1
#./fit --nevents=100000000 --tag=jacsvsmass100M_50shift --run=closure --post_tag=idx2 --degs_corr_x=10 --degs_corr_y=4 --degs_A0_x=3 --degs_A0_y=3 --degs_A1_x=3 --degs_A1_y=3 --degs_A2_x=3 --degs_A2_y=3 --degs_A3_x=3 --degs_A3_y=3 --degs_A4_x=3 --degs_A4_y=3 --jUL --j0 --j1 --j2 --j3 --j4 --jacmass=2

#python run.py --algo=fit_grid --tag=addmass0 --post_tag=1M  --nevents=100000000
#python run.py --algo=fit_grid --tag=addmass1 --post_tag=2M  --nevents=100000000
#python run.py --algo=fit_grid --tag=addmass2 --post_tag=4M  --nevents=100000000
#python run.py --algo=fit_grid --tag=addmass3 --post_tag=8M  --nevents=100000000
#python run.py --algo=fit_grid --tag=addmass4 --post_tag=16M --nevents=100000000
#python run.py --algo=fit_grid --tag=addmass5 --post_tag=32M --nevents=100000000
#python run.py --algo=fit_grid --tag=addmass6 --post_tag=64M --nevents=100000000
#python run.py --algo=fit_grid --tag=addmass7 --post_tag=128M --nevents=100000000
#python run.py --algo=fit_grid --tag=addmass8 --post_tag=256M --nevents=100000000
#python run.py --algo=fit_grid --tag=addmass9 --post_tag=512M --nevents=100000000
#python run.py --algo=fit_grid --tag=addmass10 --post_tag=1024M --nevents=100000000
#python run.py --algo=fit_grid --tag=addmass11 --post_tag=2048M --nevents=100000000

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
#python run.py --algo=fit --tag=addmass4 --post_tag=DEBUG --nevents=100000000 --scale0 --scale1 --scale2 --scale3

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


## Oct17
#python run.py --algo=fit --tag=jacVsM_1G --post_tag=1G --nevents=100000000
#python run.py --algo=fit --tag=jacVsM_1G --post_tag=1G_SCALE_ALL_BUT_A4  --nevents=100000000 --scale0 --scale1 --scale2 --scale3
#python run.py --algo=fit --tag=jacVsM_4G --post_tag=4G --nevents=100000000
#python run.py --algo=fit --tag=jacVsM_4G --post_tag=4G_SCALE_ALL_BUT_A4  --nevents=100000000 --scale0 --scale1 --scale2 --scale3
#python run.py --algo=fit --tag=jacVsM_4G --post_tag=4G_SCALE_ALL_BUT_A4_rebin21  --nevents=100000000 --scale0 --scale1 --scale2 --scale3 --rebinX=2 --rebinY=1
#python run.py --algo=fit --tag=jacVsM_4G --post_tag=4G_SCALE_ALL_BUT_A4_rebin12  --nevents=100000000 --scale0 --scale1 --scale2 --scale3 --rebinX=1 --rebinY=2
#python run.py --algo=fit --tag=jacVsM_4G --post_tag=4G_SCALE_ALL_BUT_A4_rebin22  --nevents=100000000 --scale0 --scale1 --scale2 --scale3 --rebinX=2 --rebinY=2
#python run.py --algo=fit --tag=jacVsM_4G --post_tag=4G_SCALE_ALL_BUT_A4_rebin32  --nevents=100000000 --scale0 --scale1 --scale2 --scale3 --rebinX=3 --rebinY=2
#python run.py --algo=fit --tag=jacVsM_4G --post_tag=4G_SCALE_ALL_BUT_A4_rebin23  --nevents=100000000 --scale0 --scale1 --scale2 --scale3 --rebinX=2 --rebinY=3
#python run.py --algo=fit --tag=jacVsM_4G --post_tag=4G_SCALE_ALL_BUT_A4_rebin33  --nevents=100000000 --scale0 --scale1 --scale2 --scale3 --rebinX=3 --rebinY=3
#python run.py --algo=fit_systs --tag=SYSTS_1G --post_tag=1G --nevents=100000000
#python run.py --algo=fit_systs --tag=SYSTS_1G --post_tag=1G_SCALE_ALL_BUT_A4  --nevents=100000000 --scale0 --scale1 --scale2 --scale3

python run.py --algo=fit --tag=jacVsM_4G --post_tag=4G_rebin21  --nevents=100000000  --rebinX=2 --rebinY=1
python run.py --algo=fit --tag=jacVsM_4G --post_tag=4G_rebin12  --nevents=100000000  --rebinX=1 --rebinY=2
python run.py --algo=fit --tag=jacVsM_4G --post_tag=4G_rebin22  --nevents=100000000  --rebinX=2 --rebinY=2
python run.py --algo=fit --tag=jacVsM_4G --post_tag=4G_rebin32  --nevents=100000000  --rebinX=3 --rebinY=2
python run.py --algo=fit --tag=jacVsM_4G --post_tag=4G_rebin23  --nevents=100000000  --rebinX=2 --rebinY=3
python run.py --algo=fit --tag=jacVsM_4G --post_tag=4G_rebin33  --nevents=100000000  --rebinX=3 --rebinY=3


