import ROOT

import os.path
from sys import argv
argv.append( '-b-' )
ROOT.gROOT.SetBatch(True)
argv.remove( '-b-' )

def plot_fitopt(var='corrxy', tag='testall_UL_10_4_A0_3_4_A1_3_4_A2_3_4_A3_3_4_A4_3_4_closure', post_tag='test', legend='', offset=0.1, xmin=0, xmax=1.2):
    fit_opts        = ['jUL', 'jUL_j0', 'jUL_j0_j1', 'jUL_j0_j1_j2', 'jUL_j0_j1_j2_j3',  'jUL_j0_j1_j2_j4', 'jUL_j0_j1_j2_j3_j4']
    #fit_opts        = ['jUL_j0_j1_j2_j3_j4']
    fit_opts_labels = ['UL', 'UL, A_{0}', 'UL, A_{0}, A_{1}', 'UL, A_{0}, A_{1}, A_{2}', 'UL, A_{0}, A_{1}, A_{2}, A_{3}',  'UL, A_{0}, A_{1}, A_{2}, A_{4}', 'UL, A_{0}, A_{1}, A_{2}, A_{3}, A_{4}' ]
    histos = {}
    for fit_opt in fit_opts:
        if 'A' in var and var[-1] not in fit_opt: continue
        fin = ROOT.TFile('root/fit_'+tag+'_'+fit_opt+'_'+post_tag+'.root', 'READ')
        histos[fit_opt] = {
            'fit' : fin.Get('fit_'+var).Clone('fit_'+var+'_'+fit_opt),
            'MC'  : fin.Get('fitMC_'+var).Clone('fit_'+var+'_'+fit_opt),            
        }
    #print(histos)
    c = ROOT.TCanvas("c", "canvas", 1200, 600)
    c.cd()
    leg1 = ROOT.TLegend(0.65 if var[0]!='A' else 0.15 ,0.60, 0.85 if var[0]!='A' else 0.45, 0.88, "","brNDC")
    leg1.SetFillStyle(0)
    leg1.SetBorderSize(0)
    leg1.SetTextSize(0.03)
    leg1.SetFillColor(10)
    leg1.SetHeader(legend)
    count = 0
    mg = ROOT.TMultiGraph()
    for k,key in enumerate(fit_opts):
        if 'A' in var and var[-1] not in key: continue
        h = histos[key]['fit']
        hMC = histos[key]['MC']
        hMC.SetMarkerColor(1)
        hMC.SetMarkerStyle(3)        
        for i in range(0, h.GetN() ):
            h.SetPointX(i, h.GetPointX(i)+offset*count)
        h.SetMarkerColor(count+1)
        h.SetMarkerStyle(2)
        h.SetLineColor(count+1)
        h.SetLineWidth(2)
        if count==0:
            leg1.AddEntry(hMC, 'Truth', 'LP')
            mg.Add(hMC)
        mg.Add(h)
        leg1.AddEntry(h, fit_opts_labels[k], 'LP')
        count += 1
    mg.GetYaxis().SetRangeUser(xmin,xmax)
    mg.GetYaxis().SetTitle(var)
    mg.GetXaxis().SetTitle('POI index')
    mg.Draw('ap')
    leg1.Draw()    
    #raw_input()
    c.SaveAs('root/plot_'+var+'_'+tag+'_'+post_tag+'.png')

plot_x_inty = True
plot_Ai = True

if plot_x_inty:
    plot_fitopt(var='x_inty', tag='dev0_UL_10_4_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure', post_tag='rebin11', legend='100M events, 36/120', offset=0.003)
    plot_fitopt(var='x_inty', tag='dev0_UL_10_4_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure', post_tag='rebin12', legend='100M events, 36/60', offset=0.003)
    plot_fitopt(var='x_inty', tag='dev0_UL_10_4_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure', post_tag='rebin13', legend='100M events, 36/40', offset=0.003)
    plot_fitopt(var='x_inty', tag='dev0_UL_10_4_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure', post_tag='rebin14', legend='100M events, 36/30', offset=0.003)
    plot_fitopt(var='x_inty', tag='dev0_UL_10_4_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure', post_tag='rebin16', legend='100M events, 36/20', offset=0.003)
    plot_fitopt(var='x_inty', tag='dev0_UL_10_4_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure', post_tag='rebin21', legend='100M events, 18/120', offset=0.003)
    plot_fitopt(var='x_inty', tag='dev0_UL_10_4_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure', post_tag='rebin31', legend='100M events, 12/120', offset=0.003)
if plot_Ai:
    for Ai in ['A0', 'A1', 'A2', 'A3', 'A4']:
        plot_fitopt(var=Ai, tag='dev0_UL_10_4_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure', post_tag='rebin11', legend='100M events, 36/120', offset=0.1, xmin=-1.5, xmax=1.5)
        plot_fitopt(var=Ai, tag='dev0_UL_10_4_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure', post_tag='rebin12', legend='100M events, 36/60', offset=0.1, xmin=-1.5, xmax=1.5)
        plot_fitopt(var=Ai, tag='dev0_UL_10_4_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure', post_tag='rebin13', legend='100M events, 36/40', offset=0.1, xmin=-1.5, xmax=1.5)
        plot_fitopt(var=Ai, tag='dev0_UL_10_4_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure', post_tag='rebin14', legend='100M events, 36/30', offset=0.1, xmin=-1.5, xmax=1.5)
        plot_fitopt(var=Ai, tag='dev0_UL_10_4_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure', post_tag='rebin16', legend='100M events, 36/20', offset=0.1, xmin=-1.5, xmax=1.5)
        plot_fitopt(var=Ai, tag='dev0_UL_10_4_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure', post_tag='rebin21', legend='100M events, 18/120', offset=0.1, xmin=-1.5, xmax=1.5)
        plot_fitopt(var=Ai, tag='dev0_UL_10_4_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure', post_tag='rebin31', legend='100M events, 12/120', offset=0.1, xmin=-1.5, xmax=1.5)

#plot_fitopt(var='A4', tag='polAifix_UL_10_4_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure')
