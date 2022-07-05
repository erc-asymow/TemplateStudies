import ROOT

def plot_fitopt(var='corrxy', tag='testall_UL_10_4_A0_3_4_A1_3_4_A2_3_4_A3_3_4_A4_3_4_closure'):
    #fit_opts = ['jUL', 'jUL_j0', 'jUL_j0_j1', 'jUL_j0_j1_j2', 'jUL_j0_j1_j2_j3', 'jUL_j0_j1_j2_j4', 'jUL_j0_j1_j2_j3_j4', 'jUL_j0_j4', 'jUL_j0_j3']
    #fit_opts = ['jUL', 'jUL_j0', 'jUL_j0_j3']
    #fit_opts = ['jUL', 'jUL_j0', 'jUL_j0_j4']
    fit_opts = ['jUL', 'jUL_j0', 'jUL_j0_j1', 'jUL_j0_j1_j2', 'jUL_j0_j1_j2_j3', 'jUL_j0_j1_j2_j3_j4']
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
    leg1 = ROOT.TLegend(0.15,0.75,0.35,0.88, "","brNDC")
    leg1.SetFillStyle(0)
    leg1.SetBorderSize(0)
    leg1.SetTextSize(0.04)
    leg1.SetFillColor(10)
    count = 0
    mg = ROOT.TMultiGraph()
    for key in fit_opts:
        h = histos[key]['fit']
        for i in range(0, h.GetN() ):
            h.SetPointX(i, h.GetPointX(i)+0.1*count)
        h.SetMarkerColor(count+1)
        h.SetMarkerStyle(2)
        h.SetLineColor(count+1)
        h.SetLineWidth(2)
        mg.Add(h)
        leg1.AddEntry(h, key, 'LP')
        count += 1
    mg.GetYaxis().SetRangeUser(0.0,1.0)
    mg.Draw('ap')
    leg1.Draw()
    raw_input()

#plot_fitopt(var='corrxy', tag='testall_UL_10_4_A0_2_3_A1_2_3_A2_2_3_A3_3_3_A4_3_4_closure')
plot_fitopt(var='A4', tag='testall_y2p5_UL_10_4_A0_3_4_A1_3_4_A2_3_4_A3_3_4_A4_3_4_closure')
#plot_fitopt(var='corrxy')
