import argparse
import os
import sys
import copy
import math
import ROOT

from multiprocessing import Process
   
parser = argparse.ArgumentParser(description='run')

parser.add_argument('--none', action='store_true'  , help = 'none')
parser.add_argument('--dryrun', action='store_true'  , help = 'dry run')
parser.add_argument('--algo',   default='all'        , help = 'algo')
parser.add_argument('--batch', action='store_true'  , help = 'bath')
parser.add_argument('--mt', action='store_true'  , help = 'mt')

args = parser.parse_args()

if args.batch:
    from sys import argv
    argv.append( '-b-' )
    ROOT.gROOT.SetBatch(True)
    argv.remove( '-b-' )

procs = {
     'A0' : {
         'deg_x' : [2,3,4,5,6,7,8,9,10,11],
         'deg_y' : [2,4,6,8,10,12,14],
         'opts'   : {
             'opt1' : {
                 'cmd' : '--run=wp --extrabinsX=5 --extrabinsY=5 --cULx=1 --dULx=20 --dULy=10 --doA0 --cA0x=0',
                 'tag' : 'wp_A0'
             },
             'opt2' : {
                 'cmd' : '--run=wm --extrabinsX=5 --extrabinsY=5 --cULx=1 --dULx=20 --dULy=10 --doA0 --cA0x=0',
                 'tag' : 'wm_A0'
             },
             'opt3' : {
                 'cmd' : '--run=z --extrabinsX=5 --extrabinsY=5 --cULx=1 --dULx=20 --dULy=10 --doA0 --cA0x=0',
                 'tag' : 'z_A0'
             },
         },
     },
    'A1' : {
         'deg_x' : [1,2,3,4,5,6,7,8,9,10,11],
         'deg_y' : [1,3,5,7,9,11],
         'opts'   : {
             'opt1' : {
                 'cmd' : '--run=wp --extrabinsX=5 --extrabinsY=5 --cULx=1 --dULx=20 --dULy=10 --doA1 --cA1x=1',
                 'tag' : 'wp_A1'
             },
             'opt2' : {
                 'cmd' : '--run=wm --extrabinsX=5 --extrabinsY=5 --cULx=1 --dULx=20 --dULy=10 --doA1 --cA1x=1',
                 'tag' : 'wm_A1'
             },
             'opt3' : {
                 'cmd' : '--run=z --extrabinsX=5 --extrabinsY=5 --cULx=1 --dULx=20 --dULy=10 --doA1 --cA1x=1',
                 'tag' : 'z_A1'
             },
         },
     },
     'A2' : {
         'deg_x' : [2,3,4,5,6,7,8,9,10,11],
         'deg_y' : [2,4,6,8,10,12,14],
         'opts'   : {
             'opt1' : {
                 'cmd' : '--run=wp --extrabinsX=5 --extrabinsY=5 --cULx=1 --dULx=20 --dULy=10 --doA2 --cA2x=1',
                 'tag' : 'wp_A2'
             },
             'opt2' : {
                 'cmd' : '--run=wm --extrabinsX=5 --extrabinsY=5 --cULx=1 --dULx=20 --dULy=10 --doA2 --cA2x=1',
                 'tag' : 'wm_A2'
             },
             'opt3' : {
                 'cmd' : '--run=z --extrabinsX=5 --extrabinsY=5 --cULx=1 --dULx=20 --dULy=10 --doA2 --cA2x=1',
                 'tag' : 'z_A2'
             },
         },
     },
    'A3' : {
         'deg_x' : [2,3,4,5,6,7,8,9,10,11],
         'deg_y' : [2,4,6,8,10,12,14],
         'opts'   : {
             'opt1' : {
                 'cmd' : '--run=wp --extrabinsX=5 --extrabinsY=5 --cULx=1 --dULx=20 --dULy=10 --doA3 --cA3x=1',
                 'tag' : 'wp_A3'
             },
             'opt2' : {
                 'cmd' : '--run=wm --extrabinsX=5 --extrabinsY=5 --cULx=1 --dULx=20 --dULy=10 --doA3 --cA3x=1',
                 'tag' : 'wm_A3'
             },
             'opt3' : {
                 'cmd' : '--run=z --extrabinsX=5 --extrabinsY=5 --cULx=1 --dULx=20 --dULy=10 --doA3 --cA3x=1',
                 'tag' : 'z_A3'
             },
         },
     },
    'A4' : {
         'deg_x' : [2,3,4,5,6,7,8,9,10,11],
         'deg_y' : [1,3,5,7,9,11],
         'opts'   : {
             'opt1' : {
                 'cmd' : '--run=wp --extrabinsX=5 --extrabinsY=5 --cULx=1 --dULx=20 --dULy=10 --doA4 --cA4x=0',
                 'tag' : 'wp_A4'
             },
             'opt2' : {
                 'cmd' : '--run=wm --extrabinsX=5 --extrabinsY=5 --cULx=1 --dULx=20 --dULy=10 --doA4 --cA4x=0',
                 'tag' : 'wm_A4'
             },
             'opt3' : {
                 'cmd' : '--run=z --extrabinsX=5 --extrabinsY=5 --cULx=1 --dULx=20 --dULy=10 --doA4 --cA4x=0',
                 'tag' : 'z_A4'
             },
         },
     },
    'UL' : {
        'deg_y' : [6,8,10,12,14],
        'deg_x' : [8,10,12,14,16,18,20,22,24,26],
        'opts'   : {
             'opt1' : {
                 'cmd' : '--cULx=1 --run=wp --extrabinsX=5 --extrabinsY=5',
                 'tag' : 'wp_UL'
             },
             'opt2' : {
                 'cmd' : '--cULx=1 --run=wm --extrabinsX=5 --extrabinsY=5',
                 'tag' : 'wm_UL'
             },
             'opt3' : {
                 'cmd' : '--cULx=1 --run=z --extrabinsX=5 --extrabinsY=5',
                 'tag' : 'z_UL'
             },
         },
     },
}

#allowed_procs = ["A1","A2","A3","A4"]
allowed_procs = ["A2","A3", "A4"]

def plot_pvals(fnames=[], metric="pvals", proc="wp"):
    ROOT.gStyle.SetPadRightMargin(0.15)
    histo = ROOT.TH2D('histo_pvals', '', 30, 0,30,20,0,20)
    histo.GetXaxis().SetTitle("deg_{x}")
    histo.GetYaxis().SetTitle("deg_{y}")
    histo.SetStats(0)
    for ifname,fname in enumerate(fnames):
        if proc not in fname[1]:
            continue
        if metric=="pvalue":
            histo.SetMaximum(1.0)
            histo.SetMinimum(0.01)
            histo.SetTitle(fname[0]+' '+fname[1].split('_')[0]+', FOM: p-value')
        elif metric=="ratio":
            histo.SetMaximum(0.02)
            histo.SetMinimum(0.00)
            histo.SetTitle(fname[0]+' '+fname[1].split('_')[0]+', FOM: max(fit/mc - 1)')
        elif metric=="delta":
            histo.SetMaximum(+0.05)
            histo.SetMinimum(0.00)
            histo.SetTitle(fname[0]+' '+fname[1].split('_')[0]+', FOM: max(fit - mc)')
        f = ROOT.TFile('fout_'+fname[1]+'.root', 'READ')        
        print(f.GetName())
        h = f.Get(fname[0]+'/h_info_'+fname[0])
        pval = h.GetBinContent(3)
        dx = h.GetBinContent(6)
        dy = h.GetBinContent(7)
        idx = histo.GetXaxis().FindBin(dx)
        idy = histo.GetYaxis().FindBin(dy)
        if metric=="ratio":
            h2 = f.Get(fname[0]+'/h_ratio_'+fname[0])
            pval = ROOT.TMath.Abs(h2.GetMaximum() - 1)
        elif metric=="delta":
            h2 = f.Get(fname[0]+'/h_delta_'+fname[0])
            pval = ROOT.TMath.Abs( ROOT.TMath.Max(h2.GetMaximum(), h2.GetMinimum() ) )
        histo.SetBinContent(idx,idy,pval)
        print('pval %s' % pval)
        f.Close()
    c = ROOT.TCanvas("c", "canvas", 600, 600);
    c.cd()
    histo.Draw("colz")    
    outname = 'closure_'+fnames[0][0]+'_'+proc+'_'+metric+'.png'
    if args.batch:
        print('Saving image as '+outname)
        c.SaveAs(outname)
    else:
        input()


def run_one_opt(procs,iproc,opt):
        for dx in procs[iproc]['deg_x']:
            for dy in procs[iproc]['deg_y']:
                fname = procs[iproc]['opts'][opt]['tag']+'_x'+str(dx)+'_y'+str(dy)
                command = './getparam --tag='+fname+' '+procs[iproc]['opts'][opt]['cmd']+' --d'+iproc+'x='+str(dx)+' --d'+iproc+'y='+str(dy)
                os.system(command)                

def run_all(procs):
    ps = []        
    fnames = []
    counter = 0
    for iproc in procs.keys():
        if iproc not in allowed_procs:
            continue
        for opt in procs[iproc]['opts'].keys():
            if args.mt:
                p = Process(target=run_one_opt, args=(procs,iproc,opt))
                p.start()
                ps.append(p)
                continue
            for dx in procs[iproc]['deg_x']:
                for dy in procs[iproc]['deg_y']:
                    fname = procs[iproc]['opts'][opt]['tag']+'_x'+str(dx)+'_y'+str(dy)
                    counter += 1
                    if args.algo=='run':
                        command = './getparam --tag='+fname+' '+procs[iproc]['opts'][opt]['cmd']+' --d'+iproc+'x='+str(dx)+' --d'+iproc+'y='+str(dy)
                        print(command)
                        if args.dryrun:
                            continue
                        else:                        
                            os.system(command)
                    elif args.algo=='plot':
                        fnames.append([iproc,fname])
    if args.algo=='run' and args.mt:    
        for p in ps: 
            p.join()
    print('Submitted %d proc' % counter)
    return fnames


if __name__ == '__main__':    
    fnames = run_all(procs)
    if args.algo=='plot':
        #plot_pvals(fnames, metric="pvals", proc="wp")
        for fom in ['pvalue', 'delta']:
            for iproc in ["wp","wm","z"]:
                plot_pvals(fnames, metric=fom, proc=iproc)
      
