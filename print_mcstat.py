import argparse
import os
import sys
import copy
import math
import ROOT
import numpy as np

def error(entry, key):
    if 'FC' in key:
        n = getattr(entry, 'nFC')
    else:
        n = getattr(entry, 'n')
    cov = getattr(entry, key+'_cov')
    relerr = math.sqrt(cov*(1-cov)/n)
    return relerr

def print_one( fname = 'mcstat_1_20_0p015_decorr.root'):
    f = ROOT.TFile( './root/' + fname, 'READ' )
    t = f.Get('treesum')
    lumiscale = float(fname.split('_')[1])
    nbins = int(fname.split('_')[2])
    asym = float((fname.split('_')[3]).split('p')[1])*1e-03
    print('Doing file '+fname)
    vals = {
        'data'                : "Gaus Hessian                        ",
        'data5s'              : "Gaus Hessian (5\\sigma)              ",
        'dataBB'              : "Gaus BB-lite Hessian                ",
        'data5sBB'            : "Gaus BB-lite Hessian (5\\sigma)      ",
        'dataPoisBB'          : "Pois BB-lite Hessian                ",
        'data5sPoisBB'        : "Pois BB-lite Hessian (5\\sigma)      ",
        'dataPoisBBfull'      : "Pois BB-full Hessian                ",
        'data5sPoisBBfull'    : "Pois BB-full Hessian (5\\sigma)      ",
        'dataPoisBBfullPLR'   : "Pois BB-full PLR                    ",
        'data5sPoisBBfullPLR' : "Pois BB-full PLR (5\\sigma)          ",
        'dataPoisBBfullFC'    : "Pois BB-full FC                     ",
        'data5sPoisBBfullFC'  : "Pois BB-full FC (5\\sigma)           ",
        }
    if "Barlett" in fname:   
        vals['dataPoisBBfullFC'] = "Pois BB-full FC + Barlet            "
        vals['data5sPoisBBfullFC'] = "Pois BB-full FC + Barlet (5\\sigma)  "
    if "Cheat" in fname:
        vals['dataPoisBBfullFC'] = "Pois BB-full FC Cheat               "
        vals['data5sPoisBBfullFC'] = "Pois BB-full FC Cheat (5\\sigma)     "
        
    for entry in t:        
        print( "nbins=%d, lumiscale=%.0f, asym=%.3f" % ( nbins, lumiscale, asym)  )
        print( "corr=%.3f, cond=%.0f, true error=%.3f" % (entry.asym_corr, entry.asym_cond, entry.asym_err ) )
        for key in vals.keys():
            print( vals[key] + (' | %.3f \\pm %.3f' % (getattr(entry, key+'_cov'), error(entry, key) ) ) + ' | ' + '%.4f \\pm %.4f (%.4f)' % (getattr(entry, key+'_err'), getattr(entry, key+'_derr'), getattr(entry, key+'_med' ))  )
        
    return

if __name__ == '__main__':
    print_one( 'mcstat_1_200_0p015_decorr_Barlett.root' )
