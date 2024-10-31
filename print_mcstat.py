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
        'data'   : "Data            ",
        'data5s' : "Data 5\\sigma    "
        }
    for entry in t:
        print( "nbins=%d, lumiscale=%.0f, asym=%.3f" % ( nbins, lumiscale, asym)  )
        print( "corr=%.3f, cond=%.0f, true error=%.3f" % (entry.asym_corr, entry.asym_cond, entry.asym_err ) )
        for key in vals.keys():
            print( vals[key] + ' | ' + '%.4f \\pm %.4f (%.4f)' % (getattr(entry, key+'_err'), getattr(entry, key+'_derr'), getattr(entry, key+'_med' )) + (' | %.3f \\pm %.3f' % (getattr(entry, key+'_cov'), error(entry, key) ) ) )
        
    return

if __name__ == '__main__':
    print_one( 'mcstat_1_20_0p015_decorr.root' )
