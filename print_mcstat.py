import argparse
import os
import sys
import math
import copy
import math
import ROOT
import numpy as np
from pprint import pprint

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

def print_coverage( val, err, tolerance):
    err *= 1e+03
    str = (' $%.3f(%.0f)$' % (val, err) )    
    if tolerance>0. and math.fabs(val-0.683)<=tolerance*0.683:
        str = (' ${\\bf %.3f}(%.0f)$' % (val, err) )
    str = str.replace('(0)', '(1)')
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
        'dataPois'            : "Gauss            & Numeric  & Hessian        & 0          ",
        'data5sPois'          : "                 &          &                & $5\\sigma$  ",
        'data'                : "Gauss            & Analytic & Hessian        & 0          ",
        'data5s'              : "                 &          &                & $5\\sigma$  ",
        'dataBB'              : "Gauss + MC stat. & Analytic & Hessian        & 0          ",
        'data5sBB'            : "                 &          &                & $5\\sigma$  ",
        'dataPoisBB'          : "Gauss + BB-lite  & Numeric  & Hessian        & 0          ",
        'data5sPoisBB'        : "                 &          &                & $5\\sigma$  ",
        'dataPoisBBfull'      : "Gauss + BB       & Numeric  & Hessian        & 0          ",
        'data5sPoisBBfull'    : "                 &          &                & $5\\sigma$  ",
        'dataPoisBBfullPLR'   : "Gauss + BB       & Numeric  & PLR scan       & 0          ",
        'data5sPoisBBfullPLR' : "                 &          &                & $5\\sigma$  ",
        'dataPoisBBfullFC'    : "Gauss + BB       & Numeric  & All-profiled FC& 0          ",
        'data5sPoisBBfullFC'  : "                 &          &                & $5\\sigma$  ",
        }

    f_FCfixToTrue = ROOT.TFile( './root/mcstat_' + fname + '_FCfixToTrue.root', 'READ' )
    t_FCfixToTrue = f_FCfixToTrue.Get('treesum')
    vals_FCfixToTrue = {
        'dataPoisBBfullFC'    : "Gauss + BB       & Numeric  & NP-profiled FC & 0          ",
        'data5sPoisBBfullFC'  : "                 &          &                & $5\\sigma$  "
        }

    f_Poisson = ROOT.TFile( './root/mcstat_' + fname + '_Poisson.root', 'READ' )
    t_Poisson = f_Poisson.Get('treesum')
    vals_Poisson = {
        'dataPois'            : "Poisson          & Numeric  & Hessian        & 0          ",
        'data5sPois'          : "                 &          &                & $5\\sigma$  ",
        }

    f_Barlett = ROOT.TFile( './root/mcstat_' + fname + '_BarlettFCfixToTrue.root', 'READ' )
    t_Barlett = f_Barlett.Get('treesum')
    vals_Barlett = {
        'dataPoisBBfullFC'    : "Gauss + BB       & Numeric  & PLR+Barlett    & 0          ",
        'data5sPoisBBfullFC'  : "                 &          &                & $5\\sigma$  "
        }

    f_FCCheat = ROOT.TFile( './root/mcstat_' + fname + '_FCCheat.root', 'READ' )
    t_FCCheat = f_FCCheat.Get('treesum')
    vals_FCCheat = {
        'dataPoisBBfullFC'    : "Gauss + BB       & Numeric  & FC Cheat       & 0          ",
        'data5sPoisBBfullFC'  : "                 &          &                & $5\\sigma$  "
        }

    f_Jtilde = ROOT.TFile( './root/mcstat_' + fname + '_JtildeFCfixToTrue.root', 'READ' )
    t_Jtilde = f_Jtilde.Get('treesum')
    vals_Jtilde = {
        'dataPoisBBfullFC'    : "Gauss + BB       & Numeric  & {\\it A posteriori} HC     & 0        ",
        'data5sPoisBBfullFC'  : "                 &          &                           & $5\\sigma$"    
        }

    print('')
    print('')
    print( '\\begin{tabular}{cccccc}' )
    print( 'Likelihood     & Minimim. & CI method      & $\\Delta\\mu^{\\rm true }$ & Coverage    &  $1\\sigma$ error (mean, median)            \\\\')
    print( '\\hline')
    print( '\\hline')

    for entry in t_nom:                
        print(  'Gauss + infinite MC& Analytic & Hessian      & 0           &  ' + ('%.3f' % entry.asym_cov) + '      & ' + (' $%.3f, \\; $%.3f$' % (entry.asym_err, entry.asym_err)) + ' \\\\'  )
        #print( '               &          &             & $5\\sigma$   &  ' + ('%.3f' % entry.asym_cov) + '      & ' + (' $%.3f$' % entry.asym5s_err) + ' \\\\'  )
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
                
    for entry in t_FCfixToTrue:
        for key in vals_FCfixToTrue.keys():
            print( vals_FCfixToTrue[key] + ' & ' + print_coverage( getattr(entry, key+'_cov'), error(entry, key) ) + ' & ' + print_mean( getattr(entry, key+'_err'), getattr(entry, key+'_derr'), getattr(entry, key+'_med' )) + ' \\\\' )
            if 'sigma' in vals_FCfixToTrue[key]:
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
        print('\\caption{ Results for $n=%d$, $k=%.0f$, and $\\alpha=%.3f$. The correlation between the transformed POI\'s is $\\mathrm{Corr}(\\mu_{0}^{\\prime},\\mu_{1}^{\\prime})=%.2f$ and the condition number of the covariance matrix is $\\mathrm{Cond}(V)=%.0f$.}' % ( nbins, lumiscale, asym, entry.asym_corr, entry.asym_cond))
    else:
        print('\\caption{ Results for $n=%d$, $k=%.0f$, and $\\alpha=%.3f$. The correlation between the nominal POI\'s is $\\mathrm{Corr}{\\mu_{0}^{\\prime},\\mu_{1}^{\\prime})=%.5f$ and the condition number of the covariance matrix is $\\mathrm{Cond}(V)=%.0f$.}' % ( nbins, lumiscale, asym, entry.asym_corr, entry.asym_cond))

    print( '' )
    print( '' )
    for entry in t_nom:
        asym   = getattr(entry, 'asym_err')
        bblite = getattr(entry, 'dataBB_err')
        bbfull = getattr(entry, 'dataPoisBBfull_err')
        adhoc = asym * math.sqrt( 1 + 1./lumiscale ) * bbfull/bblite
        print('Ad hoc error: %.3f' % (adhoc))
    
    return

def print_all_pm( fname = '1_200_0p015_decorr', tolerance=0.03):

    f_nom = ROOT.TFile( './root/mcstat_' + fname + '.root', 'READ' )
    t_nom = f_nom.Get('treesum')

    fm5_nom = ROOT.TFile( './root/mcstat_' + fname + 'm5s.root', 'READ' )
    tm5_nom = fm5_nom.Get('treesum')

    lumiscale = float(fname.split('_')[0])
    nbins = int(fname.split('_')[1])
    asym = float((fname.split('_')[2]).split('p')[1])*1e-03
    print('Doing files '+fname)
    vals_nom = {
        'dataPois'            : "Gauss            & Numeric  & Hessian        & 0                    ",
        'data5sPois'          : "                 &          &                &$+5\\sigma_{\\infty }$  ",
        'data'                : "Gauss            & Analytic & Hessian        & 0                    ",
        'data5s'              : "                 &          &                &$+5\\sigma_{\\infty }$  ",
        'dataBB'              : "Gauss + MC stat. & Analytic & Hessian        & 0                    ",
        'data5sBB'            : "                 &          &                &$+5\\sigma_{\\infty }$  ",
        'dataPoisBB'          : "Gauss + BB-lite  & Numeric  & Hessian        & 0                    ",
        'data5sPoisBB'        : "                 &          &                &$+5\\sigma_{\\infty }$  ",
        'dataPoisBBfull'      : "Gauss + BB       & Numeric  & Hessian        & 0                    ",
        'data5sPoisBBfull'    : "                 &          &                &$+5\\sigma_{\\infty }$  ",
        'dataPoisBBfullPLR'   : "Gauss + BB       & Numeric  & PLR scan       & 0                    ",
        'data5sPoisBBfullPLR' : "                 &          &                &$+5\\sigma_{\\infty }$  ",
        'dataPoisBBfullFC'    : "Gauss + BB       & Numeric  & All-profiled FC& 0                    ",
        'data5sPoisBBfullFC'  : "                 &          &                &$+5\\sigma_{\\infty }$  ",
        }

    f_FCfixToTrue = ROOT.TFile( './root/mcstat_' + fname + '_FCfixToTrue.root', 'READ' )
    t_FCfixToTrue = f_FCfixToTrue.Get('treesum')
    fm5_FCfixToTrue = ROOT.TFile( './root/mcstat_' + fname + 'm5s_FCfixToTrue.root', 'READ' )
    tm5_FCfixToTrue = fm5_FCfixToTrue.Get('treesum')
    vals_FCfixToTrue = {
        'dataPoisBBfullFC'    : "Gauss + BB       & Numeric  & NP-profiled FC & 0                    ",
        'data5sPoisBBfullFC'  : "                 &          &                &$+5\\sigma_{\\infty }$  "
        }

    f_Poisson = ROOT.TFile( './root/mcstat_' + fname + '_Poisson.root', 'READ' )
    t_Poisson = f_Poisson.Get('treesum')
    fm5_Poisson = ROOT.TFile( './root/mcstat_' + fname + 'm5s_Poisson.root', 'READ' )
    tm5_Poisson = fm5_Poisson.Get('treesum')
    vals_Poisson = {
        'dataPois'            : "Poisson          & Numeric  & Hessian        & 0                    ",
        'data5sPois'          : "                 &          &                &$+5\\sigma_{\\infty }$  ",
        }

    f_Barlett = ROOT.TFile( './root/mcstat_' + fname + '_BarlettFCfixToTrue.root', 'READ' )
    t_Barlett = f_Barlett.Get('treesum')
    fm5_Barlett = ROOT.TFile( './root/mcstat_' + fname + 'm5s_BarlettFCfixToTrue.root', 'READ' )
    tm5_Barlett = fm5_Barlett.Get('treesum')
    vals_Barlett = {
        'dataPoisBBfullFC'    : "Gauss + BB       & Numeric  & PLR+Barlett    & 0                    ",
        'data5sPoisBBfullFC'  : "                 &          &                &$+5\\sigma_{\\infty }$  "
        }

    f_FCCheat = ROOT.TFile( './root/mcstat_' + fname + '_FCCheat.root', 'READ' )
    t_FCCheat = f_FCCheat.Get('treesum')
    fm5_FCCheat = ROOT.TFile( './root/mcstat_' + fname + 'm5s_FCCheat.root', 'READ' )
    tm5_FCCheat = fm5_FCCheat.Get('treesum')
    vals_FCCheat = {
        'dataPoisBBfullFC'    : "Gauss + BB       & Numeric  & FC Cheat       & 0                    ",
        'data5sPoisBBfullFC'  : "                 &          &                &$+5\\sigma_{\\infty }$  "
        }

    f_Jtilde = ROOT.TFile( './root/mcstat_' + fname + '_JtildeFCfixToTrue.root', 'READ' )
    t_Jtilde = f_Jtilde.Get('treesum')
    fm5_Jtilde = ROOT.TFile( './root/mcstat_' + fname + 'm5s_JtildeFCfixToTrue.root', 'READ' )
    tm5_Jtilde = fm5_Jtilde.Get('treesum')
    vals_Jtilde = {
        'dataPoisBBfullFC'    : "Gauss + BB       & Numeric  & {\\it A posteriori} HC     & 0                  ",
        'data5sPoisBBfullFC'  : "                 &          &                           &$+5\\sigma_{\\infty }$"    
        }

    f_adhoc = ROOT.TFile( './root/mcstat_' + fname + '.root', 'READ' )
    t_adhoc = f_adhoc.Get('treesum')
    fm5_adhoc = ROOT.TFile( './root/mcstat_' + fname + 'm5s.root', 'READ' )
    tm5_adhoc = fm5_adhoc.Get('treesum')
    vals_adhoc = {
        'adhocdata'    : "Gauss + BB       & Numeric  & {\\it Ad hoc}   & 0                  ",
        'adhocdata5s'  : "                 &          &                &$+5\\sigma_{\\infty }$"    
        }

    print('')
    print('')
    print( '\\begin{tabular}{cccccc}' )
    print( 'Likelihood     & Minimim. & CI method      & $\\Delta\\mu$ & Coverage    &  $\\hat{\\sigma}$ (mean, median)          \\\\')
    print( '\\hline')
    print( '\\hline')

    for entry in t_nom:                
        print(  'Gauss (asympt.)  & Analytic & Hessian        & 0                     &  ' + ('$%.3f$' % entry.asym_cov)  + '    & ' + (' $%.3f, \\; %.3f$' % (entry.asym_err, entry.asym_err)) + ' \\\\'  )
    for entry in t_nom:                
        print(  '                 &          &                &$+5\\sigma_{\\infty }$   &  ' + ('$%.3f$' % entry.asym_cov) + '    & ' + (' $%.3f, \\; %.3f$' % (entry.asym5s_err, entry.asym5s_err)) + ' \\\\'  )
    for entry in tm5_nom:                
        print(  '                 &          &                &$-5\\sigma_{\\infty }$   &  ' + ('$%.3f$' % entry.asym_cov) + '    & ' + (' $%.3f, \\; %.3f$' % (entry.asym5s_err, entry.asym5s_err)) + ' \\\\'  )
    print( '\\hline')
    
    for entry in t_Poisson:
        for key in vals_Poisson.keys():
            if 'sigma' in vals_Poisson[key]:
                print( vals_Poisson[key] + ' & ' + print_coverage( getattr(entry, key+'_cov'), error(entry, key), tolerance ) + ' & ' + print_mean( getattr(entry, key+'_err'), getattr(entry, key+'_derr'), getattr(entry, key+'_med' )) + ' \\\\' )
                for entrym5 in tm5_Poisson:
                    print( vals_Poisson[key].replace('$+5', '$-5') + ' & ' + print_coverage( getattr(entrym5, key+'_cov'), error(entrym5, key), tolerance ) + ' & ' + print_mean( getattr(entrym5, key+'_err'), getattr(entrym5, key+'_derr'), getattr(entrym5, key+'_med' )) + ' \\\\' )
                print( '\\hline')
            else:
                print( vals_Poisson[key] + ' & ' + print_coverage( getattr(entry, key+'_cov'), error(entry, key), tolerance ) + ' & ' + print_mean( getattr(entry, key+'_err'), getattr(entry, key+'_derr'), getattr(entry, key+'_med' )) + ' \\\\' )

    for entry in t_nom:
        for key in vals_nom.keys():
            if 'sigma' in vals_nom[key]:
                print( vals_nom[key] + ' & ' + print_coverage( getattr(entry, key+'_cov'), error(entry, key), tolerance ) + ' & ' + print_mean( getattr(entry, key+'_err'), getattr(entry, key+'_derr'), getattr(entry, key+'_med' )) + ' \\\\' )
                for entrym5 in tm5_nom:
                    print( vals_nom[key].replace('$+5', '$-5') + ' & ' + print_coverage( getattr(entrym5, key+'_cov'), error(entrym5, key), tolerance ) + ' & ' + print_mean( getattr(entrym5, key+'_err'), getattr(entrym5, key+'_derr'), getattr(entrym5, key+'_med' )) + ' \\\\' )  
                print( '\\hline')
            else:
                print( vals_nom[key] + ' & ' + print_coverage( getattr(entry, key+'_cov'), error(entry, key), tolerance ) + ' & ' + print_mean( getattr(entry, key+'_err'), getattr(entry, key+'_derr'), getattr(entry, key+'_med' )) + ' \\\\' )
            
                    
    for entry in t_FCfixToTrue:
        for key in vals_FCfixToTrue.keys():
            if 'sigma' in vals_FCfixToTrue[key]:
                print( vals_FCfixToTrue[key] + ' & ' + print_coverage( getattr(entry, key+'_cov'), error(entry, key), tolerance ) + ' & ' + print_mean( getattr(entry, key+'_err'), getattr(entry, key+'_derr'), getattr(entry, key+'_med' )) + ' \\\\' )
                for entrym5 in tm5_FCfixToTrue:
                    print( vals_FCfixToTrue[key].replace('$+5', '$-5') + ' & ' + print_coverage( getattr(entrym5, key+'_cov'), error(entrym5, key), tolerance ) + ' & ' + print_mean( getattr(entrym5, key+'_err'), getattr(entrym5, key+'_derr'), getattr(entrym5, key+'_med' )) + ' \\\\' )  
                print( '\\hline')
            else:
                print( vals_FCfixToTrue[key] + ' & ' + print_coverage( getattr(entry, key+'_cov'), error(entry, key), tolerance ) + ' & ' + print_mean( getattr(entry, key+'_err'), getattr(entry, key+'_derr'), getattr(entry, key+'_med' )) + ' \\\\' )
                

    for entry in t_Barlett:
        for key in vals_Barlett.keys():
            if 'sigma' in vals_Barlett[key]:
                print( vals_Barlett[key] + ' & ' + print_coverage( getattr(entry, key+'_cov'), error(entry, key), tolerance ) + ' & ' + print_mean( getattr(entry, key+'_err'), getattr(entry, key+'_derr'), getattr(entry, key+'_med' )) + ' \\\\' )
                for entrym5 in tm5_Barlett:
                    print( vals_Barlett[key].replace('$+5', '$-5') + ' & ' + print_coverage( getattr(entrym5, key+'_cov'), error(entrym5, key), tolerance ) + ' & ' + print_mean( getattr(entrym5, key+'_err'), getattr(entrym5, key+'_derr'), getattr(entrym5, key+'_med' )) + ' \\\\' )
                print( '\\hline')
            else:
                print( vals_Barlett[key] + ' & ' + print_coverage( getattr(entry, key+'_cov'), error(entry, key), tolerance ) + ' & ' + print_mean( getattr(entry, key+'_err'), getattr(entry, key+'_derr'), getattr(entry, key+'_med' )) + ' \\\\' )

    for entry in t_Jtilde:
        for key in vals_Jtilde.keys():
            if 'sigma' in vals_Jtilde[key]:
                print( vals_Jtilde[key] + ' & ' + print_coverage( getattr(entry, key+'_cov'), error(entry, key), tolerance ) + ' & ' + print_mean( getattr(entry, key+'_err'), getattr(entry, key+'_derr'), getattr(entry, key+'_med' )) + ' \\\\' )
                for entrym5 in tm5_Jtilde:
                    print( vals_Jtilde[key].replace('$+5', '$-5') + ' & ' + print_coverage( getattr(entrym5, key+'_cov'), error(entrym5, key), tolerance ) + ' & ' + print_mean( getattr(entrym5, key+'_err'), getattr(entrym5, key+'_derr'), getattr(entrym5, key+'_med' )) + ' \\\\' ) 
                print( '\\hline')
            else:
                print( vals_Jtilde[key] + ' & ' + print_coverage( getattr(entry, key+'_cov'), error(entry, key), tolerance ) + ' & ' + print_mean( getattr(entry, key+'_err'), getattr(entry, key+'_derr'), getattr(entry, key+'_med' )) + ' \\\\' )
                
    for entry in t_FCCheat:
        for key in vals_FCCheat.keys():
            if 'sigma' in vals_FCCheat[key]:
                print( vals_FCCheat[key] + ' & ' + print_coverage( getattr(entry, key+'_cov'), error(entry, key), tolerance ) + ' & ' + print_mean( getattr(entry, key+'_err'), getattr(entry, key+'_derr'), getattr(entry, key+'_med' )) + ' \\\\' )
                for entrym5 in tm5_FCCheat:
                    print( vals_FCCheat[key].replace('$+5', '$-5') + ' & ' + print_coverage( getattr(entrym5, key+'_cov'), error(entrym5, key), tolerance ) + ' & ' + print_mean( getattr(entrym5, key+'_err'), getattr(entrym5, key+'_derr'), getattr(entrym5, key+'_med' )) + ' \\\\' )
                print( '\\hline')
            else:
                print( vals_FCCheat[key] + ' & ' + print_coverage( getattr(entry, key+'_cov'), error(entry, key), tolerance ) + ' & ' + print_mean( getattr(entry, key+'_err'), getattr(entry, key+'_derr'), getattr(entry, key+'_med' )) + ' \\\\' )

    for entry in t_adhoc:
        for key in vals_adhoc.keys():
            if 'sigma' in vals_adhoc[key]:
                print( vals_adhoc[key] + ' & ' + print_coverage( getattr(entry, key+'_cov'), error(entry, key), tolerance ) + ' & ' + print_mean( getattr(entry, key+'_err'), getattr(entry, key+'_derr'), getattr(entry, key+'_med' )) + ' \\\\' )
                for entrym5 in tm5_adhoc:
                    print( vals_adhoc[key].replace('$+5', '$-5') + ' & ' + print_coverage( getattr(entrym5, key+'_cov'), error(entrym5, key), tolerance ) + ' & ' + print_mean( getattr(entrym5, key+'_err'), getattr(entrym5, key+'_derr'), getattr(entrym5, key+'_med' )) + ' \\\\' )
                print( '\\hline')
            else:
                print( vals_adhoc[key] + ' & ' + print_coverage( getattr(entry, key+'_cov'), error(entry, key), tolerance ) + ' & ' + print_mean( getattr(entry, key+'_err'), getattr(entry, key+'_derr'), getattr(entry, key+'_med' )) + ' \\\\' )

                
    print( '\\hline')
    print( '\\end{tabular}' )
    if 'decorr' in fname:        
        #print('\\caption{ Results for $n=%d$, $k=%.0f$, and $\\alpha=%.3f$. The correlation between the transformed POI\'s is $\\mathrm{Corr}(\\mu_{0}^{\\prime},\\mu_{1}^{\\prime})=%.2f$ and the condition number of the covariance matrix is $\\mathrm{Cond}(V)=%.0f$.}' % ( nbins, lumiscale, asym, entry.asym_corr, entry.asym_cond))
        #print('\\caption{ Results for $n=%d$, $k=%.0f$, and $\\alpha=%.3f$. The correlation between the transformed POI\'s is $\\mathrm{Corr}(\\mu_{0}^{\\prime},\\mu_{1}^{\\prime})=%.2f$ and the condition number of the covariance matrix is $\\mathrm{Cond}(V)=%.0f$.}' % ( nbins, lumiscale, asym, entry.asym_corr, entry.asym_cond))
        print(('\\caption{ Coverage, mean, and median of the $1\\sigma$ confidence interval for $\\mu=\\mu^{\\prime}_{0}$ for the toy model with $n=%d$, $k=%.0f$, and $\\alpha=%.3f$. Cases where the observed coverage agrees with the expectation of 0.683 within $%.0f' % ( nbins, lumiscale, asym, tolerance*1e+02) ) +'\\%$ are highlighted with a bold font.}' )
    else:
        print(('\\caption{ Coverage, mean, and median of the $1\\sigma$ confidence interval for $\\mu=\\mu_{0}$ for the toy model with $n=%d$, $k=%.0f$, and $\\alpha=%.3f$. Cases where the observed coverage agrees with the expectation of 0.683 within $%.0f' % ( nbins, lumiscale, asym, tolerance*1e+02) ) +'\\%$ are highlighted with a bold font.}' )

    print( '' )
    print( '' )
    for entry in t_nom:
        asym   = getattr(entry, 'asym_err')
        bblite = getattr(entry, 'dataBB_err')
        bbfull = getattr(entry, 'dataPoisBBfull_err')
        adhoc = asym * math.sqrt( 1 + 1./lumiscale ) * bbfull/bblite
        print('Ad hoc error: %.3f' % (adhoc))    
    return



if __name__ == '__main__':
    #print_one( args.file )
    print_all_pm( args.file, tolerance=0.05 )
