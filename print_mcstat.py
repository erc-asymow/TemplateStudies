import argparse
import os
import sys
import copy
import math
import ROOT
import numpy as np

parser = argparse.ArgumentParser(description='run')

parser.add_argument('--none', action='store_true'  , help = 'none')
parser.add_argument('--dryrun', action='store_true'  , help = 'dry run')
parser.add_argument('-f', '--file', default='1_20_0p015_decorr'  , help = '')

args = parser.parse_args()

def error(entry, key):
    if 'FC' in key:
        n = getattr(entry, 'nFC')
    else:
        n = getattr(entry, 'n')
    cov = getattr(entry, key+'_cov')
    relerr = math.sqrt(cov*(1-cov)/n)
    return relerr

def print_one( fname = '1_20_0p015_decorr'):
    f = ROOT.TFile( './root/mcstat_' + fname + '.root', 'READ' )
    t = f.Get('treesum')
    lumiscale = float(fname.split('_')[0])
    nbins = int(fname.split('_')[1])
    asym = float((fname.split('_')[2]).split('p')[1])*1e-03
    print('Doing file '+fname)
    vals = {
        'data'                : "Gaus           & Analytic & Hessian     & 0          ",
        'data5s'              : "               &          &             & $5\\sigma$  ",
        'dataPois'            : "Gaus           & Numeric  & Hessian     & 0          ",
        'data5sPois'          : "               &          & Hessian     & $5\\sigma$  ",
        'dataBB'              : "Gaus + BB-lite & Analytic & Hessian     & 0          ",
        'data5sBB'            : "               &          & Hessian     & $5\\sigma$  ",
        'dataPoisBB'          : "Gaus + BB-lite & Numeric  & Hessian     & 0          ",
        'data5sPoisBB'        : "               &          & Hessian     & $5\\sigma$  ",
        'dataPoisBBfull'      : "Gaus + BB-full & Numeric  & Hessian     & 0          ",
        'data5sPoisBBfull'    : "               &          & Hessian     & $5\\sigma$  ",
        'dataPoisBBfullPLR'   : "Gaus + BB-full & Numeric  & PLR scan    & 0          ",
        'data5sPoisBBfullPLR' : "               &          & PLR scan    & $5\\sigma$  ",
        'dataPoisBBfullFC'    : "Gaus + BB-full & Numeric  & FC Profile  & 0          ",
        'data5sPoisBBfullFC'  : "               &          & FC Profile  & $5\\sigma$  ",
        }
    if "Poisson" in fname:   
        vals['dataPois']   =    "Poisson        & Numeric  & Hessian     & 0          "
        vals['data5sPois'] =    "               & Numeric  & Hessian     & $5\\sigma$  "
    if "Barlett" in fname:   
        vals['dataPoisBBfullFC']   =   "Gaus + BB-full & Numeric  & PLR+Barlett & 0          "
        vals['data5sPoisBBfullFC'] =   "               & Numeric  & PLR+Barlett & $5\\sigma$  "
    if "Cheat" in fname:
        vals['dataPoisBBfullFC']   =   "Gaus + BB-full & Numeric  & FC Cheat    & 0          "
        vals['data5sPoisBBfullFC'] =   "               & Numeric  & FC Cheat    & $5\\sigma$  "
    if "Jtilde" in fname:
        vals['dataPoisBBfullFC']   =   "Gaus + BB-full & Numeric  & FC+\\tilde{J}& 0          "
        vals['data5sPoisBBfullFC'] =   "               & Numeric  & FC+\\tilde{J}& $5\\sigma$  "    

    for entry in t:                
        print( '\\begin{tabular}{cccccc}' )
        print( 'Likelihood     & Minimim. & CI method   & $\\mu_{0}$ & coverage          & mean error (median)              \\\\')
        print( '\\hline')
        #print( "\\rho=%.5f, Cond(V_{true})=%.0f, \\sigma_{true}=%.4f, asymprotic coverage=%.4f" % (entry.asym_corr, entry.asym_cond, entry.asym_err, entry.asym_cov ) )
        for key in vals.keys():
            print( vals[key] + (' & $%.3f \\pm %.3f$' % (getattr(entry, key+'_cov'), error(entry, key) ) ) + ' & ' + '$%.5f \\pm %.5f \\; (%.3f)$ \\\\' % (getattr(entry, key+'_err'), getattr(entry, key+'_derr'), getattr(entry, key+'_med' ))  )
            if 'sigma' in vals[key]:
                print( '\\hline')
        print( '\\hline')
        print( '\\end{tabular}' )
        print('\\caption{ Results for $n_{bins}=%d$, mc/data=%.0f, $\\alpha=%.3f$, $\\rho=%.5f$, $Cond(V_{true})=%.0f$, $\\sigma_{true}=%.4f$, asymptotic coverage=%.4f}' % ( nbins, lumiscale, asym, entry.asym_corr, entry.asym_cond, entry.asym_err, entry.asym_cov))
    return

if __name__ == '__main__':
    print_one( args.file )
