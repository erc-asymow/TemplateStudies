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
parser.add_argument('--jac', action='store_true'  , help = 'jac')
parser.add_argument('--fit', action='store_true'  , help = 'fit')
parser.add_argument('--syst', action='store_true'  , help = 'syst')
parser.add_argument('--doA0', action='store_true'  , help = '')
parser.add_argument('--doA1', action='store_true'  , help = '')
parser.add_argument('--doA2', action='store_true'  , help = '')
parser.add_argument('--doA3', action='store_true'  , help = '')
parser.add_argument('--doA4', action='store_true'  , help = '')
parser.add_argument('--systname', default='scet'  , help = '')
parser.add_argument('--posttag', default=''  , help = '')
parser.add_argument('--xf_max', dest = 'xf_max'  , type = float,  default=0.4, help='')
parser.add_argument('--yf_max', dest = 'yf_max'  , type = float,  default=4.0, help='')

args = parser.parse_args()

if args.batch:
    from sys import argv
    argv.append( '-b-' )
    ROOT.gROOT.SetBatch(True)
    argv.remove( '-b-' )

procs = {
     'A0' : {
         'deg_x' : [2,3,4,5,6,7,8],
         'deg_y' : [2,4,6,8,10,12],
         #'fit_deg_y' : [2,4,6],
         #'fit_deg_x' : [1,2,3,4],
         'fit_deg_y' : [2],
         'fit_deg_x' : [1,2],         
         'opts'   : {
             'opt_wp' : {
                 'cmd' : '--run=wp --extrabinsX=10 --extrabinsY=10 --cULx=1 --dULx=20 --dULy=10 --doA0 --cA0x=0',
                 'tag' : 'wp_A0',
                 'nom_deg_x' : 8,
                 'nom_deg_y' : 10,
                 'syst_deg_x' : 3,
                 'syst_deg_y' : 2,
                 'cmd_syst': '--doA0 --cA0x=0 --shift_A0=0.5'
             },
             'opt_wm' : {
                 'cmd' : '--run=wm --extrabinsX=10 --extrabinsY=10 --cULx=1 --dULx=20 --dULy=10 --doA0 --cA0x=0',
                 'tag' : 'wm_A0',
                 'nom_deg_x' : 8,
                 'nom_deg_y' : 8,
                 'syst_deg_x' : 3,
                 'syst_deg_y' : 2,
                 'cmd_syst': '--doA0 --cA0x=0 --shift_A0=0.5'
             },
             'opt_z' : {
                 'cmd' : '--run=z --extrabinsX=10 --extrabinsY=10 --cULx=1 --dULx=20 --dULy=10 --doA0 --cA0x=0',
                 'tag' : 'z_A0',
                 'nom_deg_x' : 8,
                 'nom_deg_y' : 6,
                 'syst_deg_x' : 3,
                 'syst_deg_y' : 2,
                 'cmd_syst': '--doA0 --cA0x=0 --shift_A0=0.5'
             },
         },
     },
    'A1' : {
        'deg_x' : [1,2,3,4,5,6,7,8],
        'deg_y' : [1,3,5,7,9],
        #'fit_deg_y' : [2,4,6],
        #'fit_deg_x' : [1,2,3,4],
        'fit_deg_y' : [2],
        'fit_deg_x' : [1,2],
        'opts'   : {
             'opt_wp' : {
                 'cmd' : '--run=wp --extrabinsX=10 --extrabinsY=10 --cULx=1 --dULx=20 --dULy=10 --doA1 --cA1x=1',
                 'tag' : 'wp_A1',
                 'nom_deg_x' : 8,
                 'nom_deg_y' : 9,
                 'syst_deg_x' : 3,
                 'syst_deg_y' : 2,
                 'cmd_syst': '--doA1 --cA1x=1 --shift_A1=0.5'
             },
            'opt_wm' : {
                'cmd' : '--run=wm --extrabinsX=10 --extrabinsY=10 --cULx=1 --dULx=20 --dULy=10 --doA1 --cA1x=1',
                'tag' : 'wm_A1',
                'nom_deg_x' : 8,
                'nom_deg_y' : 7,
                'syst_deg_x' : 3,
                'syst_deg_y' : 2,
                'cmd_syst': '--doA1 --cA1x=1 --shift_A1=0.5'
             },
             'opt_z' : {
                 'cmd' : '--run=z --extrabinsX=10 --extrabinsY=10 --cULx=1 --dULx=20 --dULy=10 --doA1 --cA1x=1',
                 'tag' : 'z_A1',
                 'nom_deg_x' : 8,
                 'nom_deg_y' : 7,
                 'syst_deg_x' : 3,
                 'syst_deg_y' : 2,
                 'cmd_syst': '--doA1 --cA1x=1 --shift_A1=0.5'
             },
         },
     },
     'A2' : {
         'deg_x' : [2,3,4,5,6,7,8],
         'deg_y' : [2,4,6,8,10,12],
         #'fit_deg_y' : [2,4,6],
         #'fit_deg_x' : [1,2,3,4],
         'fit_deg_y' : [2],
         'fit_deg_x' : [1,2],
         'opts'   : {
             'opt_wp' : {
                 'cmd' : '--run=wp --extrabinsX=10 --extrabinsY=10 --cULx=1 --dULx=20 --dULy=10 --doA2 --cA2x=1',
                 'tag' : 'wp_A2',
                 'nom_deg_x' : 7,
                 'nom_deg_y' : 6,
                 'syst_deg_x' : 3,
                 'syst_deg_y' : 2,
                 'cmd_syst': '--doA2 --cA2x=1 --shift_A2=0.5'
             },
             'opt_wm' : {
                 'cmd' : '--run=wm --extrabinsX=10 --extrabinsY=10 --cULx=1 --dULx=20 --dULy=10 --doA2 --cA2x=1',
                 'tag' : 'wm_A2',
                 'nom_deg_x' : 7,
                 'nom_deg_y' : 6,
                 'syst_deg_x' : 3,
                 'syst_deg_y' : 2,
                 'cmd_syst': '--doA2 --cA2x=1 --shift_A2=0.5'
             },
             'opt_z' : {
                 'cmd' : '--run=z --extrabinsX=10 --extrabinsY=10 --cULx=1 --dULx=20 --dULy=10 --doA2 --cA2x=1',
                 'tag' : 'z_A2',
                 'nom_deg_x' : 7,
                 'nom_deg_y' : 6,
                 'syst_deg_x' : 3,
                 'syst_deg_y' : 2,
                 'cmd_syst': '--doA2 --cA2x=1 --shift_A2=0.5'
             },
         },
     },
    'A3' : {
        'deg_x' : [2,3,4,5,6,7,8],
        'deg_y' : [2,4,6,8,10,12],
        #'fit_deg_y' : [2,4,6,8,10],
        #'fit_deg_x' : [1,2,3,4],
        'fit_deg_y' : [4],
        'fit_deg_x' : [1,2],
        'opts'   : {
             'opt_wp' : {
                 'cmd' : '--run=wp --extrabinsX=10 --extrabinsY=10 --cULx=1 --dULx=20 --dULy=10 --doA3 --cA3x=1 --cA3y=0',
                 'tag' : 'wp_A3',
                 'nom_deg_x' : 6,
                 'nom_deg_y' : 6,
                 'syst_deg_x' : 3,
                 'syst_deg_y' : 4,
                 #'syst_deg_y' : 6,
                 #'cmd_syst': '--doA3 --cA3x=1 --cA3y=0 --shift_A3=0.5'
                 'cmd_syst': '--doA3 --cA3x=1 --cA3y=0 --shift_A3=0.1 --syst_as_additive_A3'
             },
             'opt_wm' : {
                 'cmd' : '--run=wm --extrabinsX=10 --extrabinsY=10 --cULx=1 --dULx=20 --dULy=10 --doA3 --cA3x=1 --cA3y=0',
                 'tag' : 'wm_A3',
                 'nom_deg_x' : 6,
                 'nom_deg_y' : 6,
                 'syst_deg_x' : 3,
                 'syst_deg_y' : 4,
                 #'syst_deg_y' : 6,
                 #'cmd_syst': '--doA3 --cA3x=1 --cA3y=0 --shift_A3=0.5'
                 'cmd_syst': '--doA3 --cA3x=1 --cA3y=0 --shift_A3=0.1 --syst_as_additive_A3'
             },
             'opt_z' : {
                 'cmd' : '--run=z --extrabinsX=10 --extrabinsY=10 --cULx=1 --dULx=20 --dULy=10 --doA3 --cA3x=1 --cA3y=0',
                 'tag' : 'z_A3',
                 'nom_deg_x' : 6,
                 'nom_deg_y' : 6,
                 'syst_deg_x' : 3,
                 #'syst_deg_y' : 6,
                 #'cmd_syst': '--doA3 --cA3x=1 --cA3y=0 --shift_A3=0.5'
                 'syst_deg_y' : 4,
                 'cmd_syst': '--doA3 --cA3x=1 --cA3y=0 --shift_A3=0.1 --syst_as_additive_A3'
             },
         },
     },
    'A4' : {
        'deg_x' : [2,3,4,5,6,7,8],
        'deg_y' : [1,3,5,7,9,11,13],
        #'fit_deg_y' : [2,4,6,8],
        #'fit_deg_x' : [1,2,3,4],
        'fit_deg_y' : [4,6],
        'fit_deg_x' : [1,2,3],
        'opts'   : {
             'opt_wp' : {
                 'cmd' : '--run=wp --extrabinsX=10 --extrabinsY=10 --cULx=1 --dULx=20 --dULy=10 --doA4 --cA4x=0',
                 'tag' : 'wp_A4',
                 'nom_deg_x' : 7,
                 'nom_deg_y' : 13,
                 'syst_deg_x' : 4,
                 'syst_deg_y' : 8,
                 'cmd_syst': '--doA4 --cA4x=0 --shift_A4=0.2',
             },
             'opt_wm' : {
                 'cmd' : '--run=wm --extrabinsX=10 --extrabinsY=10 --cULx=1 --dULx=20 --dULy=10 --doA4 --cA4x=0',
                 'tag' : 'wm_A4',
                 'nom_deg_x' : 7,
                 'nom_deg_y' : 11,
                 'syst_deg_x' : 4,
                 'syst_deg_y' : 8,
                 'cmd_syst': '--doA4 --cA4x=0 --shift_A4=0.2',
             },
             'opt_z' : {
                 'cmd' : '--run=z --extrabinsX=10 --extrabinsY=10 --cULx=1 --dULx=20 --dULy=10 --doA4 --cA4x=0',
                 'tag' : 'z_A4',
                 'nom_deg_x' : 7,
                 'nom_deg_y' : 11,
                 'syst_deg_x' : 4,
                 'syst_deg_y' : 8,
                 'cmd_syst': '--doA4 --cA4x=0 --shift_A4=0.2',
             },
         },
     },
    'UL' : {
        'deg_y' : [6,8,10,12,14],
        'deg_x' : [8,10,12,14,16,18,20,22,24,26],
        #'fit_deg_y' : [2,4,6],
        #'fit_deg_x' : [6,7,8,9,10,11,12],
        'fit_deg_y' : [4],
        'fit_deg_x' : [5,6],
        'opts'   : {
             'opt_wp' : {
                 'cmd' : '--cULx=1 --run=wp --extrabinsX=10 --extrabinsY=10',
                 'tag' : 'wp_UL',
                 'nom_deg_x' : 20,
                 'nom_deg_y' : 10,
                 'syst_deg_x' : 8,
                 'syst_deg_y' : 6,
                 'cmd_syst' : '--run=wp --extrabinsX=10 --extrabinsY=10 --cULx=1 --shift_UL=0.1',
             },
             'opt_wm' : {
                 'cmd' : '--cULx=1 --run=wm --extrabinsX=10 --extrabinsY=10',
                 'tag' : 'wm_UL',
                 'nom_deg_x' : 20,
                 'nom_deg_y' : 10,
                 'syst_deg_x' : 8,
                 'syst_deg_y' : 6,
                 'cmd_syst' : '--run=wm --extrabinsX=10 --extrabinsY=10 --cULx=1 --shift_UL=0.1',
             },
             'opt_z' : {
                 'cmd' : '--cULx=1 --run=z --extrabinsX=10 --extrabinsY=10',
                 'tag' : 'z_UL',
                 'nom_deg_x' : 20,
                 'nom_deg_y' : 10,
                 'syst_deg_x' : 8,
                 'syst_deg_y' : 6,
                 'cmd_syst' : '--run=z --extrabinsX=10 --extrabinsY=10 --cULx=1 --shift_UL=0.1',
             },
        },
     },
}

allowed_procs = ["UL"]
if args.doA0:
    allowed_procs = ["A0"]
if args.doA1:
    allowed_procs = ["A1"]
if args.doA2:
    allowed_procs = ["A2"]
if args.doA3:
    allowed_procs = ["A3"]
if args.doA4:
    allowed_procs = ["A4"]

if args.syst:
    allowed_procs = ["UL","A0","A1","A2","A3","A4"]
    
xf_max_str = ("%.2f" % args.xf_max).replace('.', 'p')
yf_max_str = ("%.2f" % args.yf_max).replace('.', 'p')
#print(xf_max_str)
#print(yf_max_str)

def plot_pvals(fnames=[], metric="pvalue", proc="wp"):
    ROOT.gStyle.SetPadRightMargin(0.15)
    max_deg = 15
    if 'UL' in fnames[0]:
        max_deg = 30
    histo = ROOT.TH2D('histo_pvals', '', max_deg, 0,max_deg,max_deg,0, max_deg)
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

def plot_pulls(ifnames):
    #ROOT.gStyle.SetPadRightMargin(0.15)
    dx = ifnames[1].split('_')[2].replace('x', 'd_{x}=')
    dy = ifnames[1].split('_')[3].replace('y', 'd_{y}=')
    #print(dx+' '+dy)
    f = ROOT.TFile('fout_'+ifnames[1]+'.root', 'READ')        
    print(f.GetName())
    h1 = f.Get(ifnames[0]+'/h_pull_'+ifnames[0])
    h1.SetTitle(ifnames[0]+' '+ifnames[1].split('_')[0]+', 2D pulls, p-value='+('%.2f' % f.Get(ifnames[0]+'/h_info_'+ifnames[0]).GetBinContent(3) ))
    h1.SetStats(0)
    h1.SetMinimum(-4.0)
    h1.SetMaximum(+4.0)
    h1unrol = ROOT.TH1D("h1unrol", ";pull", 40,-4,4)
    h1unrol.SetTitle(ifnames[0]+' '+ifnames[1].split('_')[0]+', unrolled pulls')
    h1unrol.SetStats(0)
    h1unrol.SetLineWidth(3)
    for ix in range(1,h1.GetXaxis().GetNbins()+1):
        for iy in range(1,h1.GetYaxis().GetNbins()+1):
            h1unrol.Fill( h1.GetBinContent(ix,iy))            
    h1unrol.Fit("gaus", "0Q")
    leg1 = ROOT.TLegend(0.10, 0.68, 0.70, 0.90, "","brNDC")
    leg1.SetFillStyle(0)
    leg1.SetBorderSize(0)
    leg1.SetTextSize(0.035)
    leg1.SetFillColor(10)
    fitf = h1unrol.GetFunction("gaus")
    mu = fitf.GetParameter(1)
    rms = fitf.GetParameter(2)
    leg1.AddEntry(fitf, "(%.2f #pm %.2f)" % (mu,rms) ,"L")    
    h2 = f.Get(ifnames[0]+'/h_data_'+ifnames[0]).Clone('h_data_'+ifnames[0]+'_scaled')
    h2.SetTitle(ifnames[0]+' '+ifnames[1].split('_')[0]+', 2D input')
    h2.SetStats(0)
    if ifnames[0]=='UL':
        for ix in range(1,h2.GetXaxis().GetNbins()+1):
            for iy in range(1,h2.GetYaxis().GetNbins()+1):
                h2.SetBinContent( ix,iy, h2.GetBinContent(ix,iy) / h2.GetXaxis().GetBinWidth(ix) / h2.GetYaxis().GetBinWidth(iy) ) 
    h3 = f.Get(ifnames[0]+'/h_pdf_'+ifnames[0])
    h3.SetTitle(ifnames[0]+' '+ifnames[1].split('_')[0]+', pdf: '+dx+', '+dy)
    h3.GetXaxis().SetTitle( h2.GetXaxis().GetTitle() )
    h3.GetYaxis().SetTitle( h2.GetYaxis().GetTitle() )
    h3.SetStats(0)
    if ifnames[0]=='UL':
        ibin = h3.FindBin(0.06, 1.0)
        pdfval = h3.GetBinContent(ibin)
        h3.SetMinimum(0)
        h3.SetMaximum(pdfval*2)
    he = f.Get(ifnames[0]+'/h_exp_'+ifnames[0])
    he.SetStats(0)
    hd = f.Get(ifnames[0]+'/h_data_'+ifnames[0])
    hd1x = hd.ProjectionX("hd1x")
    he1x = he.ProjectionX("he1x")
    hd1x.Divide(he1x)
    hd1x.SetLineWidth(2)
    hd1x.SetStats(0)
    hd1x.SetTitle(ifnames[0]+' '+ifnames[1].split('_')[0]+', projection x')
    hd1x.GetYaxis().SetTitle("fit/mc")
    hd1x.SetMaximum(1.1)
    hd1x.SetMinimum(0.9)
    if ifnames[0]=='UL':
        hd1x.SetMaximum(1.005)
        hd1x.SetMinimum(0.995)
    if ifnames[0]=='A4':
        hd1x.SetMaximum(1.02)
        hd1x.SetMinimum(0.98)
    if ifnames[0]!='UL' and 'z' in ifnames[1]:
        hd1x.SetMaximum(1.2)
        hd1x.SetMinimum(0.8)        
    hd1y = hd.ProjectionY("hd1y")
    he1y = he.ProjectionY("he1y")
    hd1y.Divide(he1y)
    hd1y.SetLineWidth(2)
    hd1y.SetStats(0)
    hd1y.GetYaxis().SetTitle("fit/mc")
    hd1y.SetTitle(ifnames[0]+' '+ifnames[1].split('_')[0]+', projection y')
    hd1y.SetMaximum(1.1)
    hd1y.SetMinimum(0.9)
    if ifnames[0]=='UL':
        hd1y.SetMaximum(1.005)
        hd1y.SetMinimum(0.995)
    if ifnames[0]=='A4':
        hd1y.SetMaximum(1.025)
        hd1y.SetMinimum(0.975)
    if ifnames[0]!='UL' and 'z' in ifnames[1]:
        hd1y.SetMaximum(1.3)
        hd1y.SetMinimum(0.7)        
    c = ROOT.TCanvas("c", "canvas", 1200, 800);
    c.Divide(3,2)
    c.cd(1)
    h2.Draw("colz")    
    c.cd(2)
    h1.Draw("colz")    
    c.cd(3)
    h3.Draw("colz")    
    c.cd(4)
    h1unrol.Draw("histe")
    fitf.Draw("same")
    leg1.Draw()
    c.cd(5)
    hd1x.Draw("histe")
    c.cd(6)
    hd1y.Draw("histe")    
    outname = 'summary_'+ifnames[1].split('_')[0]+'_'+ifnames[0]+'.png'
    if args.batch:
        print('Saving image as '+outname)
        c.SaveAs(outname)
    else:
        input()
    return


def plot_fitres(ifnames, isyst="scale"):
    #print(ifnames[0],ifnames[1])
    dx = ifnames[1].split('_')[2].replace('x', 'd_{x}=')
    dy = ifnames[1].split('_')[3].replace('y', 'd_{y}=')
    #print(dx+' '+dy)
    fin = ROOT.TFile('fout_fit_'+isyst+'_'+ifnames[0]+'_x'+xf_max_str+'_y'+yf_max_str+'_'+ifnames[1]+'.root', 'READ')        
    #print(fin.GetName())
    hinfo = fin.Get("h_info")
    hstart = fin.Get("h_start")
    ns = hinfo.GetBinContent(3)
    c = ROOT.TCanvas("c", "canvas", 800, 800)
    c.Divide(2,2)
    #do_pull = False
    #if ifnames[0] in ["UL","A0", "A1", "A2", "A3", "A4"]:
    #    do_pull = True
    #if ifnames[0]=="UL" and isyst=="scet":
    #    do_pull = False
    for i in range(0, int(ns)):        
        print("Doing %s" % i)
        #if do_pull:
        pe_i = fin.Get(("h_pull_%s" % i))
        pd_i = fin.Get(("h_data_%s" % i))            
        pd_i_clone = pd_i.Clone(("h_data_%s_copy"%i))
        for ix in range(1, pd_i_clone.GetXaxis().GetNbins()+1):
            for iy in range(1, pd_i_clone.GetYaxis().GetNbins()+1):
                pull = (pd_i.GetBinContent(ix,iy) - hstart.GetBinContent(ix,iy))/hstart.GetBinError(ix,iy)
                pd_i_clone.SetBinContent(ix,iy, pull )
        pd_i_clone.SetMinimum(-3)
        pd_i_clone.SetMaximum(+3)
        pd_i_clone.SetStats(0)
        pd_i_clone.SetTitle(ifnames[0]+", "+isyst+" "+pd_i_clone.GetTitle()+": pre-fit")
        c.cd(1)	
        pd_i_clone.Draw("colz")				  
        pe_i.SetStats(0)
        pe_i.SetMinimum(-3)
        pe_i.SetMaximum(+3)
        c.cd(2)
        ROOT.gPad.SetRightMargin(0.15)
        pe_i.SetTitle(ifnames[0]+", "+isyst+" "+pe_i.GetTitle()+": post-fit "+dx+", "+dy)
        pe_i.Draw("COLZ")
        #else:
        re_i = fin.Get(("h_exp_%s" % i))
        rd_i = fin.Get(("h_data_%s" % i))
        rd_i_clone = rd_i.Clone(("h_data_%s_copy" % i))
        rd_i_clone.Divide(hstart)
        rd_i_clone.SetMinimum(0.8)
        rd_i_clone.SetMaximum(1.2)
        rd_i_clone.SetStats(0)
        rd_i_clone.SetTitle(ifnames[0]+", "+isyst+" "+rd_i_clone.GetTitle()+": pre-fit")
        c.cd(3)
        rd_i_clone.Draw("colz")
        re_i.Divide(rd_i)
        re_i.SetStats(0)
        re_i.SetMinimum(0.97)
        re_i.SetMaximum(1.03)
        c.cd(4)
        ROOT.gPad.SetRightMargin(0.15)
        re_i.SetTitle(ifnames[0]+", "+isyst+" "+re_i.GetTitle()+": post-fit "+dx+", "+dy)
        re_i.Draw("COLZ")
        outname = "plots/"+isyst+"_var"+str(i)+"_"+ifnames[0]+'_x'+xf_max_str+'_y'+yf_max_str+'_'+ifnames[1]+'.png'
        print(outname)
        if args.batch:
            print('Saving image as '+outname)
            c.SaveAs(outname)
        else:
            input()
    fin.Close()
    return
        
def run_one_opt(procs,iproc,opt):
        for dx in procs[iproc]['deg_x']:
            for dy in procs[iproc]['deg_y']:
                fname = procs[iproc]['opts'][opt]['tag']+'_x'+str(dx)+'_y'+str(dy)
                fname += ('_'+args.posttag if args.posttag!='' else '')
                command = './getparam --outtag='+fname+' '+procs[iproc]['opts'][opt]['cmd']+' --d'+iproc+'x='+str(dx)+' --d'+iproc+'y='+str(dy)
                if args.dryrun:
                    print(command)
                else:
                    os.system(command)

                    
def run_one_opt_jac(procs,iproc,opt):
        for dx in procs[iproc]['fit_deg_x']:
            for dy in procs[iproc]['fit_deg_y']:
                fname = procs[iproc]['opts'][opt]['tag']+'_x'+str(dx)+'_y'+str(dy)
                fname += ('_'+args.posttag if args.posttag!='' else '')
                command = './getparam --outtag=jac_x'+xf_max_str+'_y'+yf_max_str+'_'+fname+' '+procs[iproc]['opts'][opt]['cmd']+\
                    ' --f'+iproc+'x='+str(dx)+' --f'+iproc+'y='+str(dy)+\
                    ' --d'+iproc+'x='+str(procs[iproc]['opts'][opt]['nom_deg_x'])+' --d'+iproc+'y='+str(procs[iproc]['opts'][opt]['nom_deg_y'])+\
                    ' --xf_max='+str(args.xf_max)+' --yf_max='+str(args.yf_max)+\
                    ' --savePdf2data --saveJac'
                if args.dryrun:
                    print(command)
                else:
                    os.system(command)

def run_one_opt_syst(procs, iopt):
    fname = ('_'+args.posttag if args.posttag!='' else '')
    command = './getparam --outtag=syst_'+iopt+'_x'+xf_max_str+'_y'+yf_max_str+fname+\
        ' --xf_max='+str(args.xf_max)+' --yf_max='+str(args.yf_max)+\
        ' --saveSyst --savePdf '
    for iproc in allowed_procs:
        command += procs[iproc]['opts']['opt_'+iopt]['cmd_syst'] +\
            ' --d'+iproc+'x='+str( procs[iproc]['opts']['opt_'+iopt]['nom_deg_x'])+' --d'+iproc+'y='+str(procs[iproc]['opts']['opt_'+iopt]['nom_deg_y'])+\
            ' --f'+iproc+'x='+str( procs[iproc]['opts']['opt_'+iopt]['syst_deg_x'] )+' --f'+iproc+'y='+str( procs[iproc]['opts']['opt_'+iopt]['syst_deg_y'] )+' '
    if args.dryrun:
        print(command)
    else:
        os.system(command)     

def run_one_opt_fit(procs,iproc,opt):
        for dx in procs[iproc]['fit_deg_x']:
            for dy in procs[iproc]['fit_deg_y']:
                fname = procs[iproc]['opts'][opt]['tag']+'_x'+str(dx)+'_y'+str(dy)
                fname += ('_'+args.posttag if args.posttag!='' else '')
                command = './getparam --intag=jac_x'+xf_max_str+'_y'+yf_max_str+'_'+fname+' --outtag=x'+xf_max_str+'_y'+yf_max_str+'_'+fname+' '+procs[iproc]['opts'][opt]['cmd']+\
                    ' --xf_max='+str(args.xf_max-0.01)+' --yf_max='+str(args.yf_max-0.01)+\
                    ' --runfit --syst_'+args.systname
                if args.dryrun:
                    print(command)
                else:
                    os.system(command)     

def run_all(procs):
    ps = []        
    fnames = []
    counter = 0

    if args.algo=='run' and args.mt and args.syst:
        for iopt in ['wp','wm','z']:
            p = Process(target=run_one_opt_syst, args=(procs,iopt))
            p.start()
            ps.append(p)
            fnames.append(iopt+('_'+args.posttag if args.posttag!='' else ''))
            counter += 1
                                    
    for iproc in procs.keys():
        if args.syst: continue
        if iproc not in allowed_procs:
            continue
        for opt in procs[iproc]['opts'].keys():
            if args.algo=='run' and args.mt:
                if args.jac:
                    p = Process(target=run_one_opt_jac, args=(procs,iproc,opt))
                    p.start()
                    ps.append(p)
                elif args.fit:
                    p = Process(target=run_one_opt_fit, args=(procs,iproc,opt))
                    p.start()
                    ps.append(p)
                else:
                    p = Process(target=run_one_opt, args=(procs,iproc,opt))
                    p.start()
                    ps.append(p)
                #continue
            deg_x_name = 'deg_x'
            if (args.jac or args.fit or args.algo=='plot_fitres'):
                deg_x_name = 'fit_deg_x'
            elif args.syst:
                deg_x_name = 'syst_deg_x'
            deg_y_name = 'deg_y'
            if (args.jac or args.fit or args.algo=='plot_fitres'):
                deg_y_name = 'fit_deg_y'
            elif args.syst:
                deg_y_name = 'syst_deg_y'
            for dx in procs[iproc][deg_x_name]:
                for dy in procs[iproc][deg_y_name]:
                    fname = procs[iproc]['opts'][opt]['tag']+'_x'+str(dx)+'_y'+str(dy)
                    fname += ('_'+args.posttag if args.posttag!='' else '')
                    counter += 1
                    fnames.append([iproc,fname])

    if args.algo=='run' and args.mt:    
        for p in ps: 
            p.join()

    print('Submitted %d proc' % counter)
    return fnames


if __name__ == '__main__':    
    fnames = run_all(procs)
    if 'plot' in args.algo:
        if args.algo=='plot_pvals':
            plot_pvals(fnames, 'ratio', 'wp')
            plot_pvals(fnames, 'ratio', 'wm')
            plot_pvals(fnames, 'ratio', 'z')
        elif args.algo=='plot_pulls':
            for iproc in procs.keys():
                if iproc not in allowed_procs:
                    continue
                for opt in procs[iproc]['opts'].keys():
                    fname = procs[iproc]['opts'][opt]['tag']+'_x'+str(procs[iproc]['opts'][opt]['nom_deg_x'])+'_y'+str(procs[iproc]['opts'][opt]['nom_deg_y'])
                    fname += ('_'+args.posttag if args.posttag!='' else '')
                    proc = iproc
                    plot_pulls([proc,fname])
        elif args.algo=='plot_fitres':
            for [proc,fname] in fnames:
                #print(fname, proc)
                plot_fitres([proc,fname], args.systname)
