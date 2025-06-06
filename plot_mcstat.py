import ROOT

import os.path
from sys import argv
#argv.append( '-b-' )
#ROOT.gROOT.SetBatch(True)
#argv.remove( '-b-' )
import math
import argparse
import numpy as np

parser = argparse.ArgumentParser(description='run')

parser.add_argument('--none', action='store_true'  , help = 'none')
parser.add_argument('--dryrun', action='store_true'  , help = 'dry run')
parser.add_argument('-f', '--file', default='1_200_0p015_corr'  , help = '')

args = parser.parse_args()


def plot_tstat(fname = "1_200_0p015_decorr"):
    c = ROOT.TCanvas("c", "canvas", 600, 600)
    c.SetLogy(1)
    leg1 = ROOT.TLegend(0.40,0.70, 0.85, 0.85, "","brNDC")
    leg1.SetFillStyle(0)
    leg1.SetBorderSize(0)
    leg1.SetTextSize(0.04)
    leg1.SetFillColor(10)
    xmin = 0.
    xmax = 10.
    fin = ROOT.TFile("root/mcstat_"+fname+".root", "READ")
    tree = fin.Get("treeFC")
    h = ROOT.TH1D("h", "; t_{#mu} ", 20, xmin, xmax)
    h.GetYaxis().SetTitleOffset(0.75)
    h.GetYaxis().SetTitleSize(0.06)
    h.GetXaxis().SetTitleOffset(0.73)
    h.GetXaxis().SetTitleSize(0.06)
    tree.Draw("tstat>>h")
    h.SetStats(0)
    h.Scale( 1./tree.GetEntries() )
    h.Scale( 1./h.GetXaxis().GetBinWidth(1) )
    h.SetMaximum(1.0)
    h.SetMinimum(0.001)
    h.SetLineWidth(2)
    c.cd()
    h.Draw("HISTE")
    leg1.AddEntry(h, 'Sampling distribution', 'LE')
    chi2 = ROOT.TF1("chi2", "1./TMath::Sqrt(2*x*TMath::Pi())*TMath::Exp(-0.5*x)", xmin, xmax)
    chi2.SetNpx(10000)
    chi2.Draw("SAME")
    chi2.SetLineStyle(ROOT.kDashed)
    leg1.AddEntry(chi2, '#chi^{2}_{1}', 'L')
    leg1.Draw()
    #c.SaveAs('tstat.png')
    #c.SaveAs('tstat.pdf')
    input()

def plot_muerr(var = "mu",  fname = "1_200_0p015_decorr"):
    c = ROOT.TCanvas("c", "canvas", 600, 600)
    c.SetLogy(1)
    leg1 = ROOT.TLegend(0.45,0.75, 0.85, 0.90, "","brNDC")
    leg1.SetFillStyle(0)
    leg1.SetBorderSize(0)
    leg1.SetTextSize(0.04)
    leg1.SetFillColor(10)
    xmin = -2 if var=="mu" else +0.05
    xmax = +4 if var=="mu" else +0.35
    fin = ROOT.TFile("root/mcstat_"+fname+".root", "READ")
    tree = fin.Get("tree")
    h = ROOT.TH1D("h", "; #hat{#mu} - 1 " if var=="mu" else "; #hat{#sigma} ", 40, xmin, xmax)
    h.GetYaxis().SetTitleOffset(0.75)
    h.GetYaxis().SetTitleSize(0.06)
    h.GetXaxis().SetTitleOffset(0.70)
    h.GetXaxis().SetTitleSize(0.06)
    h.SetLineWidth(2)
    tree.Draw(var+"_poisdata_BBfull>>h")
    #tree.Draw("errFCUp_poisdata_BBfull>>h")
    h.SetStats(0)
    #h.Scale( 1./tree.GetEntries() )
    #h.SetMaximum(1.0)
    #h.SetMinimum(0.001)
    c.cd()
    h.Draw("HISTE")
    leg1.AddEntry(h, 'Sampling distribution', 'LE')
    chi2 = ROOT.TF1("chi2", "[0]/TMath::Sqrt(2*TMath::Pi())/[1]*TMath::Exp(-0.5*(x-[2])*(x-[2])/[1]/[1])", xmin, xmax)        
    chi2.SetNpx(10000)
    chi2.SetParameter(0, h.GetEntries())
    chi2.SetParameter(1, h.GetRMS())
    chi2.SetParameter(2, h.GetMean())
    h.Fit(chi2, 'RQ')
    chi2.Draw("SAME")
    chi2.SetLineStyle(ROOT.kDashed)
    leg1.AddEntry(chi2, 'Gauss fit', 'L')
    leg1.Draw()
    c.SaveAs('hatmu.png' if var=="mu" else 'haterr.png')
    c.SaveAs('hatmu.pdf' if var=="mu" else 'haterr.pdf')
    input()

def plot_model(fname = "", eps = 0.03, N = 2e+06, k = 1, nbins = 200):
    ran = ROOT.TRandom3();
    c = ROOT.TCanvas("c", "canvas", 600, 600)
    #c.SetLogy(1)
    leg1 = ROOT.TLegend(0.70,0.75, 0.85, 0.85, "","brNDC")
    leg1.SetFillStyle(0)
    leg1.SetBorderSize(0)
    leg1.SetTextSize(0.04)
    leg1.SetFillColor(10)
    leg1.SetNColumns(2) 
    xmin = 0.
    xmax = 1.
    h1 = ROOT.TH1D("h1", "; bin index #font[12]{i}; ", nbins,0, nbins)
    h1.GetYaxis().SetNdivisions(5)
    h1.GetXaxis().SetNdivisions(5)
    h2 = ROOT.TH1D("h2", " ", nbins,0, nbins)
    h1t = ROOT.TH1D("h1t", " ", nbins,0, nbins)
    h2t = ROOT.TH1D("h2t", " ", nbins,0, nbins)
    for ibin in range(200):
        h1.SetBinContent(ibin+1, 1.0)
        h2.SetBinContent(ibin+1, 1.0+eps if ibin+1<nbins/2 else 1.0)
    totint = h1.Integral() + h2.Integral();
    h1.Scale( N/totint )
    h2.Scale( N/totint )

    for ibin in range(200):
        h1t.SetBinContent(ibin+1, ran.Poisson(h1.GetBinContent(ibin+1)*k)/k )
        h2t.SetBinContent(ibin+1, ran.Poisson(h2.GetBinContent(ibin+1)*k)/k )
        
    h1.GetYaxis().SetTitleOffset(0.75)
    h1.GetYaxis().SetTitleSize(0.06)
    h1.GetXaxis().SetTitleOffset(0.73)
    h1.GetXaxis().SetTitleSize(0.06)
    h1.SetStats(0)
    h1.SetMaximum( h1.GetMinimum()*1.08 )
    h1.SetMinimum( h1.GetMinimum()*0.94)
    h1.SetLineWidth(3)
    h1.SetLineStyle(ROOT.kDashed)
    h2.SetLineStyle(ROOT.kDashed)
    h1.SetLineColor(ROOT.kBlue)
    h2.SetLineWidth(3)
    h2.SetLineColor(ROOT.kBlack)
    h1t.SetLineWidth(1)
    h1t.SetLineColor(ROOT.kBlue)
    h2t.SetLineWidth(1)
    h2t.SetLineColor(ROOT.kBlack)
    c.cd()
    h1.Draw("HIST")
    h2.Draw("HISTSAME")
    h1t.Draw("HISTSAME")
    h2t.Draw("HISTSAME")
    leg1.AddEntry(h1t, 'T_{1i}', 'LE')
    leg1.AddEntry(h1, '#nu_{1i}', 'L')
    leg1.AddEntry(h2t, 'T_{2i}', 'LE')
    leg1.AddEntry(h2, '#nu_{2i}', 'L')
    leg1.Draw()
    #c.SaveAs('model.png')
    c.SaveAs(fname+'.pdf')
    #input()

    
def plot_simplemodel(fname = "simplemodel", var = 0.01):

    ran = ROOT.TRandom3();
    c = ROOT.TCanvas("c", "canvas", 600, 600)
    #c.SetLogy(1)
    leg1 = ROOT.TLegend(0.15,0.65, 0.55, 0.85, "","brNDC")
    leg1.SetFillStyle(0)
    leg1.SetBorderSize(0)
    leg1.SetTextSize(0.04)
    leg1.SetFillColor(10)
    #leg1.SetNColumns(2)

    nx = 50

    rho2s = np.zeros( nx, dtype = float)
    res1   = np.zeros( nx, dtype = float)
    res2   = np.zeros( nx, dtype = float)
    res3   = np.zeros( nx, dtype = float)

    for i in range(nx):
        print(i)
        rho2 = float(i)/nx
        rho2s[i] = rho2
        sin = math.sqrt(rho2)
        cos = math.sqrt(1-rho2)
        Si1 = 0.
        Si2 = 0.
        Si3 = 0.
        ntoys = 10000
        Sitrue = (cos*cos)
        for itoy in range(ntoys):
            alpha1 = ran.Gaus(0., var)
            alpha2 = ran.Gaus(0., var)
            Siitoy1 = (1 - (sin+alpha1)*(sin+alpha1)/( (sin+alpha1)*(sin+alpha1) + (cos+alpha2)*(cos+alpha2) ))
            Si1 += Siitoy1
            alpha1 *= 2
            alpha2 *= 2
            Siitoy2 = (1 - (sin+alpha1)*(sin+alpha1)/( (sin+alpha1)*(sin+alpha1) + (cos+alpha2)*(cos+alpha2) ))
            Si2 += Siitoy2
            alpha1 *= 2
            alpha2 *= 2
            Siitoy3 = (1 - (sin+alpha1)*(sin+alpha1)/( (sin+alpha1)*(sin+alpha1) + (cos+alpha2)*(cos+alpha2) ))
            Si3 += Siitoy3            
            #if i==0:
            #    print( Siitoy )
            #    print("....", Sitrue)
        Si1 /= ntoys
        Si2 /= ntoys
        Si3 /= ntoys
        res1[i] = (Si1/Sitrue)
        res2[i] = (Si2/Sitrue)
        res3[i] = (Si3/Sitrue)

    #print(rho2s)
    #print(res)
    gr1 = ROOT.TGraph( nx , rho2s, res1)
    gr2 = ROOT.TGraph( nx , rho2s, res2)
    gr3 = ROOT.TGraph( nx , rho2s, res3)
    mg = ROOT.TMultiGraph()
    mg.GetYaxis().SetRangeUser(0.98,1.10)
    mg.GetXaxis().SetRangeUser(0.,1.)
    mg.GetYaxis().SetTitle( '<#hat{S}>/S' )
    mg.GetXaxis().SetTitle( '#rho^{2}_{#mu}' )
    mg.GetYaxis().SetTitleOffset(0.71)
    mg.GetXaxis().SetTitleOffset(0.71)
    mg.GetYaxis().SetTitleSize(0.06)
    mg.GetXaxis().SetTitleSize(0.06)
    gr1.SetMarkerStyle(ROOT.kFullCircle)
    gr1.SetMarkerSize(1.)
    gr1.SetMarkerColor(ROOT.kRed)
    gr2.SetMarkerStyle(ROOT.kFullCircle)
    gr2.SetMarkerSize(1.)
    gr2.SetMarkerColor(ROOT.kBlue)
    gr3.SetMarkerStyle(ROOT.kFullCircle)
    gr3.SetMarkerSize(1.)
    gr3.SetMarkerColor(ROOT.kOrange)
    mg.Add(gr1)
    mg.Add(gr2)
    mg.Add(gr3)
    mg.Draw( 'APL' )
    leg1.SetHeader("#sqrt{<#alpha^{2}>}")
    leg1.AddEntry(gr1, "%.2f" % (var), "LP")
    leg1.AddEntry(gr2, "%.2f" % (var*2), "LP")
    leg1.AddEntry(gr3, "%.2f" % (var*4), "LP")
    line = ROOT.TF1('line','1',0, 1.)
    line.SetLineWidth(3)
    line.SetLineStyle(ROOT.kDashed)
    line.SetLineColor(ROOT.kBlack)
    line.Draw("same")
    leg1.Draw()
    c.Update()        
    #c.SaveAs(fname+'.pdf')
    input()

    
if __name__ == '__main__':
    #plot_tstat( fname=args.file )
    #plot_model( fname="model_k10", k=10)
    #plot_muerr( var = "mu", fname=args.file)
    plot_simplemodel("simplemodel", 0.03)
