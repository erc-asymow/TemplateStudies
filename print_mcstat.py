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

def print_coverage( val, err):
    err *= 1e+03
    str = (' $%.3f(%.0f)$' % (val, err) )
    return str
def print_mean( val, err, med):
    err *= 1e+04
    str = (' $%.3f, \\; %.3f$' % (val, med) )
    return str

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
        vals['data5sPois'] =    "               &          & Hessian     & $5\\sigma$  "
    if "Barlett" in fname:   
        vals['dataPoisBBfullFC']   =   "Gaus + BB-full & Numeric  & PLR+Barlett & 0          "
        vals['data5sPoisBBfullFC'] =   "               &          & PLR+Barlett & $5\\sigma$  "
    if "Cheat" in fname:
        vals['dataPoisBBfullFC']   =   "Gaus + BB-full & Numeric  & FC Cheat    & 0          "
        vals['data5sPoisBBfullFC'] =   "               &          & FC Cheat    & $5\\sigma$  "
    if "Jtilde" in fname:
        vals['dataPoisBBfullFC']   =   "Gaus + BB-full & Numeric  & FC+\\tilde{J}& 0          "
        vals['data5sPoisBBfullFC'] =   "               &          & FC+\\tilde{J}& $5\\sigma$  "    

    for entry in t:                
        print( '\\begin{tabular}{cccccc}' )
        print( 'Likelihood     & Minimim. & CI method   & $\\mu_{0}$ & coverage          & $1\\sigma$ unc.: mean, median           \\\\')
        print( '\\hline')
        for key in vals.keys():
            print( vals[key] + (' & $%.3f \\pm %.3f$' % (getattr(entry, key+'_cov'), error(entry, key) ) ) + ' & ' + '$%.5f \\pm %.5f \\; (%.3f)$ \\\\' % (getattr(entry, key+'_err'), getattr(entry, key+'_derr'), getattr(entry, key+'_med' ))  )
            if 'sigma' in vals[key]:
                print( '\\hline')
        print( '\\hline')
        print( '\\end{tabular}' )
        print('\\caption{ Results for $n_{bins}=%d$, mc/data=%.0f, $\\alpha=%.3f$, $\\rho=%.5f$, $Cond(V_{true})=%.0f$, $\\sigma_{true}=%.4f$, asymptotic coverage=%.4f}' % ( nbins, lumiscale, asym, entry.asym_corr, entry.asym_cond, entry.asym_err, entry.asym_cov))
    return


def print_all( fname = '1_200_0p015_decorr'):
    f_nom = ROOT.TFile( './root/mcstat_' + fname + '.root', 'READ' )
    t_nom = f_nom.Get('treesum')
    lumiscale = float(fname.split('_')[0])
    nbins = int(fname.split('_')[1])
    asym = float((fname.split('_')[2]).split('p')[1])*1e-03
    print('Doing files '+fname)
    vals_nom = {
        'dataPois'            : "Gaus           & Numeric  & Hessian     & 0          ",
        'data5sPois'          : "               &          &             & $5\\sigma$  ",
        'data'                : "Gaus           & Analytic & Hessian     & 0          ",
        'data5s'              : "               &          &             & $5\\sigma$  ",
        'dataBB'              : "Gaus + BB-lite & Analytic & Hessian     & 0          ",
        'data5sBB'            : "               &          &             & $5\\sigma$  ",
        'dataPoisBB'          : "Gaus + BB-lite & Numeric  & Hessian     & 0          ",
        'data5sPoisBB'        : "               &          &             & $5\\sigma$  ",
        'dataPoisBBfull'      : "Gaus + BB-full & Numeric  & Hessian     & 0          ",
        'data5sPoisBBfull'    : "               &          &             & $5\\sigma$  ",
        'dataPoisBBfullPLR'   : "Gaus + BB-full & Numeric  & PLR scan    & 0          ",
        'data5sPoisBBfullPLR' : "               &          &             & $5\\sigma$  ",
        'dataPoisBBfullFC'    : "Gaus + BB-full & Numeric  & FC Profile  & 0          ",
        'data5sPoisBBfullFC'  : "               &          &             & $5\\sigma$  ",
        }

    f_Poisson = ROOT.TFile( './root/mcstat_' + fname + '_Poisson.root', 'READ' )
    t_Poisson = f_Poisson.Get('treesum')
    vals_Poisson = {
        'dataPois'            : "Poisson        & Numeric  & Hessian     & 0          ",
        'data5sPois'          : "               &          &             & $5\\sigma$  ",
        }

    f_Barlett = ROOT.TFile( './root/mcstat_' + fname + '_Barlett.root', 'READ' )
    t_Barlett = f_Barlett.Get('treesum')
    vals_Barlett = {
        'dataPoisBBfullFC'    : "Gaus + BB-full & Numeric  & PLR+Barlett & 0          ",
        'data5sPoisBBfullFC'  : "               &          &             & $5\\sigma$  "
        }

    f_FCCheat = ROOT.TFile( './root/mcstat_' + fname + '_FCCheat.root', 'READ' )
    t_FCCheat = f_FCCheat.Get('treesum')
    vals_FCCheat = {
        'dataPoisBBfullFC'    : "Gaus + BB-full & Numeric  & FC Cheat    & 0          ",
        'data5sPoisBBfullFC'  : "               &          &             & $5\\sigma$  "
        }

    f_Jtilde = ROOT.TFile( './root/mcstat_' + fname + '_Jtilde.root', 'READ' )
    t_Jtilde = f_Jtilde.Get('treesum')
    vals_Jtilde = {
        'dataPoisBBfullFC'    : "Gaus + BB-full & Numeric  & {\\it a posteriori} HC  & 0        ",
        'data5sPoisBBfullFC'  : "               &          &                        & $5\\sigma$"    
        }

    print( '\\begin{tabular}{cccccc}' )
    print( 'Likelihood     & Minimim. & CI method   & $\\mu_{0}^{\\rm true }$   & Coverage    &  $1\\sigma$ unc.: mean, median            \\\\')
    print( '\\hline')
    print( '\\hline')

    for entry in t_nom:                
        print( 'Infinite MC    & Analytic & Hessian     & 0           &  ' + ('%.3f' % entry.asym_cov) + '      & ' + (' $%.3f$' % entry.asym_err) + ' \\\\'  )
    print( '\\hline')
    
    for entry in t_Poisson:
        for key in vals_Poisson.keys():
            print( vals_Poisson[key] + ' & ' + print_coverage( getattr(entry, key+'_cov'), error(entry, key) ) + ' & ' + print_mean( getattr(entry, key+'_err'), getattr(entry, key+'_derr'), getattr(entry, key+'_med' )) + ' \\\\' )
            if 'sigma' in vals_Poisson[key]:
                print( '\\hline')

    for entry in t_nom:                
        for key in vals_nom.keys():
            print( vals_nom[key] + ' & ' + print_coverage( getattr(entry, key+'_cov'), error(entry, key) ) + ' & ' + print_mean( getattr(entry, key+'_err'), getattr(entry, key+'_derr'), getattr(entry, key+'_med' )) + ' \\\\' )
            if 'sigma' in vals_nom[key]:
                print( '\\hline')

    for entry in t_Barlett:
        for key in vals_Barlett.keys():
            print( vals_Barlett[key] + ' & ' + print_coverage( getattr(entry, key+'_cov'), error(entry, key) ) + ' & ' + print_mean( getattr(entry, key+'_err'), getattr(entry, key+'_derr'), getattr(entry, key+'_med' )) + ' \\\\' )
            if 'sigma' in vals_Barlett[key]:
                print( '\\hline')

    for entry in t_Jtilde:
        for key in vals_Jtilde.keys():
            print( vals_Jtilde[key] + ' & ' + print_coverage( getattr(entry, key+'_cov'), error(entry, key) ) + ' & ' + print_mean( getattr(entry, key+'_err'), getattr(entry, key+'_derr'), getattr(entry, key+'_med' )) + ' \\\\' )
            if 'sigma' in vals_Jtilde[key]:
                print( '\\hline')
                
    for entry in t_FCCheat:
        for key in vals_FCCheat.keys():
            print( vals_FCCheat[key] + ' & ' + print_coverage( getattr(entry, key+'_cov'), error(entry, key) ) + ' & ' + print_mean( getattr(entry, key+'_err'), getattr(entry, key+'_derr'), getattr(entry, key+'_med' )) + ' \\\\' )
            if 'sigma' in vals_FCCheat[key]:
                print( '\\hline')
                
    print( '\\hline')
    print( '\\end{tabular}' )
    if 'decorr' in fname:        
        print('\\caption{ Results for $n_{{\\rm bins}}=%d$, ${\\rm MC/data}=%.0f$, and $\\alpha=%.3f$. The asymptotic correlation between the transformed POI\'s is $\\rho_{ {\\rm asym.}}=%.2f$ and the condition number of the covariance matrix is $\\mathrm{Cond}(V_{{\\rm asym.}})=%.0f$}.' % ( nbins, lumiscale, asym, entry.asym_corr, entry.asym_cond))
    else:
        print('\\caption{ Results for $n_{{\\rm bins}}=%d$, ${\\rm MC/data}=%.0f$, and $\\alpha=%.3f$. The asymptotic correlation between the original POI\'s is $\\rho_{ {\\rm asym.}}=%.5f$ and the condition number of the covariance matrix is $\\mathrm{Cond}(V_{{\\rm asym.}})=%.0f$}.' % ( nbins, lumiscale, asym, entry.asym_corr, entry.asym_cond))
    return



if __name__ == '__main__':
    #print_one( args.file )
    print_all( args.file )
