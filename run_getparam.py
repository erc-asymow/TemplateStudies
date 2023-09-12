import argparse
import os
import sys
import copy
import math
import ROOT

parser = argparse.ArgumentParser(description='run')

parser.add_argument('--none', action='store_true'  , help = 'none')
parser.add_argument('--dryrun', action='store_true'  , help = 'dry run')
parser.add_argument('--algo',   default='all'        , help = 'algo')

args = parser.parse_args()

procs = {
     'A0' : {
         'deg_x' : [2,3,4,5,6],
         'deg_y' : [2,4,6,8],
         'opts'   : {
             'opt1' : {
                 'cmd' : '--dULx=12 --dULy=6 --doA0 --cA0x=0',
                 'tag' : 'A0_ctrx0'
             },
             #'opt1' : {
             #    'cmd' : '--doA0 --cA0x=1 --extrabinsX=20 --extrabinsY=20',
             #    'tag' : 'A0_ctrx1_nBins20'
             #},
             #'opt2' : {
             #    'cmd' : '--doA0 --cA0x=0',
             #    'tag' : 'A0_ctrx0'
             #},
         },
     },
}

def plot_pvals(fnames=[]):
    histo = ROOT.TH2D('histo_pvals', '', 10, 0,10,10,0,10)
    for ifname,fname in enumerate(fnames):
        f = ROOT.TFile('fout_'+fname[1]+'.root', 'READ')
        print(f.GetName())
        h = f.Get(fname[0]+'/h_info_'+fname[0])
        pval = h.GetBinContent(3)
        print('pval %s' % pval)
        dx = h.GetBinContent(6)
        dy = h.GetBinContent(7)
        idx = histo.GetXaxis().FindBin(dx)
        idy = histo.GetYaxis().FindBin(dy)
        histo.SetBinContent(idx,idy,pval)
        f.Close()
    histo.Draw("colz")
    input()
        
fnames = []
counter = 0
for iproc in procs.keys():
    for dx in procs[iproc]['deg_x']:
        for dy in procs[iproc]['deg_y']:
            for opt in procs[iproc]['opts'].keys():
                fname = procs[iproc]['opts'][opt]['tag']+'_x'+str(dx)+'_y'+str(dy)
                counter += 1
                if args.algo=='run':
                    command = './getparam --tag='+fname+' '+procs[iproc]['opts'][opt]['cmd']+' --dA0x='+str(dx)+' --dA0y='+str(dy)
                    print(command)
                    if args.dryrun:
                        continue
                    else:                        
                        os.system(command)
                elif args.algo=='plot':
                    fnames.append([iproc,fname])

if args.algo=='plot':
    plot_pvals(fnames)
    
print('Submitted %d proc' % counter)
