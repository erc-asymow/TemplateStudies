import ROOT

def plot_fitopt(var='corrxy', tag='testall_UL_10_4_A0_3_4_A1_3_4_A2_3_4_A3_3_4_A4_3_4_closure'):
    fit_opts        = ['jUL', 'jUL_j0', 'jUL_j0_j1', 'jUL_j0_j1_j2', 'jUL_j0_j1_j2_j3', 'jUL_j0_j1_j2_j3_j4']
    fit_opts_labels = ['UL', 'UL, A_{0}', 'UL, A_{0}, A_{1}', 'UL, A_{0}, A_{1}, A_{2}', 'UL, A_{0}, A_{1}, A_{2}, A_{3}', 'UL, A_{0}, A_{1}, A_{2}, A_{3}, A_{4}' ]
    histos = {}
    for fit_opt in fit_opts:
        fin = ROOT.TFile('root/fit_'+tag+'_'+fit_opt+'.root', 'READ')
        histos[fit_opt] = {
            'fit' : fin.Get('fit_'+var).Clone('fit_'+var+'_'+fit_opt),
            'MC'  : fin.Get('fitMC_'+var).Clone('fit_'+var+'_'+fit_opt),            
        }
    #print(histos)
    c = ROOT.TCanvas("c", "canvas", 1200, 600)
    c.cd()
    leg1 = ROOT.TLegend(0.65,0.65,0.85,0.88, "","brNDC")
    leg1.SetFillStyle(0)
    leg1.SetBorderSize(0)
    leg1.SetTextSize(0.04)
    leg1.SetFillColor(10)
    count = 0
    mg = ROOT.TMultiGraph()
    for k,key in enumerate(fit_opts):
        h = histos[key]['fit']
        hMC = histos[key]['MC']
        hMC.SetMarkerColor(1)
        hMC.SetMarkerStyle(3)        
        for i in range(0, h.GetN() ):
            h.SetPointX(i, h.GetPointX(i)+0.1*count)
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
    mg.GetYaxis().SetRangeUser(0.0,1.2)
    mg.GetYaxis().SetTitle(var)
    mg.GetXaxis().SetTitle('POI index')
    mg.Draw('ap')
    leg1.Draw()    
    raw_input()
    c.SaveAs('root/plot_'+var+'_'+tag+'.png')

#plot_fitopt(var='x_inty', tag='polAifix_UL_10_4_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure')
plot_fitopt(var='A4', tag='polAifix_UL_10_4_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure')
#plot_fitopt(var='A4', tag='testall_y2p5_UL_10_4_A0_3_4_A1_3_4_A2_3_4_A3_3_4_A4_3_4_closure')
#plot_fitopt(var='corrxy')
