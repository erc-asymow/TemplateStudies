#!/bin/sh

# this is done once per Ai
#python run_getparam.py --algo=run --mt --jac --xf_max=0.4 --yf_max=3.5
#python run_getparam.py --algo=run --mt --jac --xf_max=0.3 --yf_max=3.0
#python run_getparam.py --algo=run --mt --fit --xf_max=0.4 --yf_max=3.5 --systname=scale_Ai
#python run_getparam.py --algo=run --mt --fit --xf_max=0.3 --yf_max=3.0 --systname=scale_Ai
#python run_getparam.py --algo=plot_fitres --xf_max=0.4 --yf_max=3.5 --systname=scale_Ai --batch
#python run_getparam.py --algo=plot_fitres --xf_max=0.3 --yf_max=3.0 --systname=scale_Ai --batch

python run_getparam.py --algo=run --mt --fit --xf_max=0.4 --yf_max=3.5 --systname=altpdf
python run_getparam.py --algo=plot_fitres --xf_max=0.4 --yf_max=3.5 --systname=altpdf --batch
