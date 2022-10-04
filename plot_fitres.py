import ROOT

import os.path
from sys import argv
argv.append( '-b-' )
ROOT.gROOT.SetBatch(True)
argv.remove( '-b-' )
import math

def plot_fitopt(var='corrxy', tag='testall_UL_10_4_A0_3_4_A1_3_4_A2_3_4_A3_3_4_A4_3_4_closure', post_tag='test', legend='', offset=0.1, xmin=0, xmax=1.2, view_err=False):
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
            if view_err: 
                h.SetPointY(i, 0.0)
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
    c.SaveAs('root/plot_'+var+'_'+tag+'_'+post_tag+('_err' if view_err else '')+'.png')

def plot_chi2(tag='testall_UL_10_4_A0_3_4_A1_3_4_A2_3_4_A3_3_4_A4_3_4_closure', post_tag='test', legend='', offset=80.0, xmin=0.0, xmax=10):
    fit_opts        = ['jUL', 
                       'jUL_j0', 
                       'jUL_j0_j1', 
                       'jUL_j0_j1_j2', 
                       'jUL_j0_j1_j2_j3',  
                       'jUL_j0_j1_j2_j4', 
                       'jUL_j0_j1_j2_j3_j4']
    fit_opts_labels = ['UL', 
                       'UL, A_{0}', 
                       'UL, A_{0,1}',
                       'UL, A_{0,1,2}',
                       'UL, A_{0,1,2,3}',
                       'UL, A_{0,1,2,4}',
                       'UL, A_{0,1,2,3,4}' ]
    histos = {}
    for fit_opt in fit_opts:
        fin = ROOT.TFile('root/fit_'+tag+'_'+fit_opt+'_'+post_tag+'.root', 'READ')
        histos[fit_opt] = {
            'fit' : fin.Get('chi2_vs_mass').Clone('chi2_vs_mass_'+fit_opt),
        }
    #print(histos)
    c = ROOT.TCanvas("c", "canvas", 800, 600)
    c.cd()
    leg1 = ROOT.TLegend(0.35 ,0.65, 0.55, 0.88, "","brNDC")
    leg1.SetFillStyle(0)
    leg1.SetBorderSize(0)
    leg1.SetTextSize(0.03)
    leg1.SetFillColor(10)
    leg1.SetHeader(legend)
    count = 0
    mg = ROOT.TMultiGraph()
    for k,key in enumerate(fit_opts):
        h = histos[key]['fit']
        for i in range(0, h.GetN() ):
            h.SetPointX(i, (h.GetPointX(i)-offset)*1e+03)
        h.SetLineColor(count+1)
        h.SetLineWidth(2)
        h.SetMarkerColor(count+1)
        h.SetMarkerStyle(2)
        h.Fit("pol2", "Q0")
        #h.Fit("pol2", "Q")
        parabola = h.GetFunction("pol2")
        param0 = parabola.GetParameter(0); 
        param1 = parabola.GetParameter(1); 
        param2 = parabola.GetParameter(2); 
        sigmaM = 1./math.sqrt(param2);
        deltaM = -param1/param2*0.5;
        min = param0 + param1*deltaM + param2*deltaM*deltaM
        for i in range(0, h.GetN() ):
            h.SetPointY(i, h.GetPointY(i)-min)
        h.GetYaxis().SetRangeUser(xmin,xmax)
        mg.Add(h)
        leg1.AddEntry(h, fit_opts_labels[k]+', #sigma='+'{:03.1f}'.format(sigmaM)+' MeV, #delta='+'{:03.1f} MeV'.format(deltaM), 'LP')
        count += 1
    mg.GetYaxis().SetRangeUser(xmin,xmax)
    mg.GetYaxis().SetTitle('#Delta#chi^{2}')
    mg.GetXaxis().SetTitle('#DeltaM_{W} (MeV)')
    mg.Draw('apl')
    leg1.Draw()    
    line = ROOT.TF1("line", "1.0", mg.GetXaxis().GetXmin(), mg.GetXaxis().GetXmax())
    line.SetLineStyle(ROOT.kDashed)
    line.SetLineWidth(3)
    line.SetLineColor(ROOT.kBlack)
    line.Draw("same")
    #input()
    c.SaveAs('root/plot_chi2vsmass'+'_'+tag+'_'+post_tag+'.png')

def plot_pulls_1D(tag='testall_UL_10_4_A0_3_4_A1_3_4_A2_3_4_A3_3_4_A4_3_4_closure', post_tag='test', legend='', xmin=0.0, xmax=10):
    histos = {}
    fin = ROOT.TFile('root/fit_'+tag+'_'+post_tag+'.root', 'READ')
    mass_dict = {
        '0'  : '79.900',
        '18' : '79.972',
        '25' : '80.000',
        '32' : '80.028',
        '49' : '80.100',
        }
    for mass_opt in mass_dict.keys():
        histos[mass_opt] = {
            'pre' : fin.Get('hw_prefit'+mass_opt).Clone('hw_prefit'+mass_opt),
            'post' : fin.Get('hw_postfit'+mass_opt).Clone('hw_postfit'+mass_opt),
        }
    #print(histos)
    c = ROOT.TCanvas("c", "canvas", 800, 600)
    c.cd()
    leg1 = ROOT.TLegend(0.15, 0.65, 0.55, 0.88, "","brNDC")
    leg1.SetFillStyle(0)
    leg1.SetBorderSize(0)
    leg1.SetTextSize(0.03)
    leg1.SetFillColor(10)
    leg1.SetHeader(legend)
    count = 0    
    for key in mass_dict.keys():
        h_pre_1d  = histos[key]['pre'].ProjectionX("h_pre_1d"+key);
        h_post_1d = histos[key]['post'].ProjectionX("h_post_1d"+key);
        h_pull = h_pre_1d.Clone("hpull"+key);
        histos[key]['pull'] = h_pull
        for ix in range(1, h_pre_1d.GetNbinsX()+1):
            pull = (h_pre_1d.GetBinContent(ix) - h_post_1d.GetBinContent(ix))/ROOT.TMath.Sqrt( h_pre_1d.GetBinContent(ix) )
            h_pull.SetBinContent(ix, pull )
            h_pull.SetBinError(ix, 0.0)
        h_pull.SetLineColor(count+1)
        h_pull.SetLineWidth(2)
        leg1.AddEntry(h_pull, 'm_{W}='+mass_dict[key]+' GeV' , 'LP')
        count += 1
    c.cd()
    count = 0    
    for key in mass_dict.keys():
        if count==0:
            histos[key]['pull'].Draw()
            histos[key]['pull'].SetStats(0)
            histos[key]['pull'].GetYaxis().SetRangeUser(xmin,xmax)
            histos[key]['pull'].GetYaxis().SetTitle('(pre-post)/#sigma_{stat.}')
            histos[key]['pull'].GetXaxis().SetTitle('p_{T} (GeV)')
        else:
            histos[key]['pull'].Draw("SAME")
        count += 1
    leg1.Draw()
    #input()
    c.SaveAs('root/pull1D'+'_'+tag+'_'+post_tag+'.png')


tags      = [#'dev0_UL_10_4_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure',
             #'dev0_UL_12_4_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure',
             #'dev0_UL_10_4_A0_3_3_A1_3_3_A2_3_3_A3_4_3_A4_3_3_closure',
             #'dev0_UL_10_4_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_4_3_closure',
    #'dev0_UL_8_4_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_4_3_closure',
    #'dev0_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_4_3_closure',
    #'dev0_UL_10_4_A0_3_3_A1_3_3_A2_3_3_A3_3_4_A4_4_3_closure',
    #'dev0_UL_10_4_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_4_closure',
    'addmass3_UL_10_4_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure'
         ]
#post_tags = ['rebin11', 'rebin12', 'rebin13', 'rebin14', 'rebin16', 'rebin21', 'rebin31' ]
post_tags = ['rebin11' ]
legends   = ['100M events, 36/60']
#legends   = ['100M events, 36/240', '100M events, 36/120', '100M events, 36/80', '100M events, 36/60', '100M events, 36/40', '100M events, 18/240', '100M events, 12/240' ]
#legends   = ['100M events, 36/120', '100M events, 36/60', '100M events, 36/40', '100M events, 36/30', '100M events, 36/20', '100M events, 18/120', '100M events, 12/120' ]

for t in tags:
    continue
    for i,pt in enumerate(post_tags):        
        plot_fitopt(var='x_inty', tag=t, post_tag=pt, legend=legends[i], offset=0.003, xmin= 0.0, xmax=1.2, view_err=False)
        plot_fitopt(var='x_inty', tag=t, post_tag=pt, legend=legends[i], offset=0.003, xmin=-0.2, xmax=0.2, view_err=True)
        plot_fitopt(var='corrxy',   tag=t, post_tag=pt, legend=legends[i], offset=0.1, xmin= 0.0, xmax=1.2, view_err=False)
        plot_fitopt(var='corrxy',   tag=t, post_tag=pt, legend=legends[i], offset=0.1, xmin=-0.2, xmax=0.2, view_err=True)
        for Ai in ['A0', 'A1', 'A2', 'A3', 'A4']:
            plot_fitopt(var=Ai,   tag=t, post_tag=pt, legend=legends[i], offset=0.1, xmin=-0.5, xmax=1.5, view_err=False)
            plot_fitopt(var=Ai,   tag=t, post_tag=pt, legend=legends[i], offset=0.1, xmin=-0.5, xmax=0.5, view_err=True)

#plot_fitopt(var='A0', tag='dev0_UL_10_4_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure', post_tag='rebin11', legend='100M events, 36/120', offset=0.003, xmin=-0.3, xmax=0.3, view_err=True)
#plot_fitopt(var='A4', tag='polAifix_UL_10_4_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure')

tags      = [#'addmass00_UL_10_4_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure',
             #'addmass0_UL_10_4_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure',
             #'addmass1_UL_10_4_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure',
             'addmass2_UL_10_4_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure',
             'addmass3_UL_10_4_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure',
             'addmass4_UL_10_4_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure',
             #'addmass5_UL_10_4_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure',
             ]
post_tags = [#'rebin11'
    '1G_SCALE_ALL_BUT_A4',
    '4G_SCALE_ALL_BUT_A4',
    '10G_SCALE_ALL_BUT_A4'
             ] 
legends   = [#'1M  --> 100M events, 36/60',
             #'10M --> 100M events, 36/60',
             #'100M--> 100M events, 36/60',
             '1G  --> 100M events, 36/60',
             '4G  --> 100M events, 36/60',
             '10G --> 100M events, 36/60',
             #'20G --> 100M events, 36/60',
             ] 


'''
tags      = ['addmass0_UL_8_4_grid',
             'addmass1_UL_8_4_grid',
             'addmass2_UL_8_4_grid',
             'addmass3_UL_8_4_grid',
             'addmass4_UL_8_4_grid',
             'addmass5_UL_8_4_grid',
             'addmass6_UL_8_4_grid',
             'addmass7_UL_8_4_grid',
             'addmass8_UL_8_4_grid',
             ]
#post_tags = ['1M', '2M', '4M', '8M', '16M', '32M', '64M', '128M', '256M'] 
post_tags = ['SCALE_ALL_BUT_A4','SCALE_ALL_BUT_A4','SCALE_ALL_BUT_A4','SCALE_ALL_BUT_A4','SCALE_ALL_BUT_A4','SCALE_ALL_BUT_A4','SCALE_ALL_BUT_A4','SCALE_ALL_BUT_A4','SCALE_ALL_BUT_A4',]
legends   = ['1M  --> 100M events, 36/60',
             '2M  --> 100M events, 36/60',
             '4M  --> 100M events, 36/60',
             '8M  --> 100M events, 36/60',
             '16M --> 100M events, 36/60',
             '32M --> 100M events, 36/60',
             '64M --> 100M events, 36/60',
             '128M --> 100M events, 36/60',
             '256M --> 100M events, 36/60',
             ] 
'''
for j,t in enumerate(tags):
    continue
    #for i,pt in enumerate(post_tags):        
    pt = post_tags[j]
    plot_chi2(tag=t, post_tag=pt, legend=legends[j], offset=80.0, xmin=0.0, xmax=5)


plot_pulls_1D(tag='addmass2_UL_10_4_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure', post_tag='jUL_j0_j1_j2_j3_j4_DEBUG', legend='Mass', xmin=-2, xmax=2)
plot_pulls_1D(tag='addmass3_UL_10_4_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure', post_tag='jUL_j0_j1_j2_j3_j4_DEBUG', legend='Mass', xmin=-2, xmax=2)
plot_pulls_1D(tag='addmass4_UL_10_4_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure', post_tag='jUL_j0_j1_j2_j3_j4_DEBUG', legend='Mass', xmin=-2, xmax=2)
