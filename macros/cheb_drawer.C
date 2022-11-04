{

  int n = 4;

  double xx[20], yy[20];
  
  auto cheb = [](double* xs, double* ps){
    double x = xs[0];
    double scale = ps[0];
    double offset = ps[1];
    double n = int(ps[2]);
    double m = int(ps[3]);
    double den = 0.;
    double num = 0.;
    for(unsigned int i = 0; i <= n ; i++){
      int sign = i%2==0 ? +1 :-1;
      double xj = (TMath::Cos((n-i)*TMath::Pi()/n) + offset)*scale;
      if(x==xj) return 1.0;
      double val = sign/(x-xj);
      if(i==0 || i==n) val *= 0.5;
      den += val;
      if(i==m) num = val;
    }
    return num/den;
  };

  double offset = 0.0;
  double scale  = 1.0;
  double xmin = (-1.0 + offset)*scale;
  double xmax = (+1.0 + offset)*scale;
  
  //TF1* f = new TF1("f", cheb, xmin, xmax, 4);
  //f->SetParameter(0, scale);
  //f->SetParameter(1, offset);
  //f->SetMinimum(-1.1);
  //f->SetMaximum(1.1);
  //f->SetParameter(2, n);
  //f->SetNpx(10000);

  TLegend* leg1 = new TLegend(0.15 ,0.75, 0.75, 0.88, "","brNDC");
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.03);
  leg1->SetFillColor(10);
  leg1->SetHeader(Form("Base functions for n=%d", n));
  leg1->SetNColumns(4);
 
  TCanvas* c = new TCanvas("c", "canvas", 800, 800);
  vector<TF1*> funcs = {};
  for(int m=0; m<=n; m++){
    TF1* f = new TF1(Form("f%m"), cheb, xmin, xmax, 4);    
    f->SetTitle("");
    f->SetParameter(0, scale);
    f->SetParameter(1, offset);
    f->SetMinimum(-0.4);
    f->SetMaximum(1.4);
    f->SetParameter(2, n);
    f->SetParameter(3, m);
    f->SetNpx(10000);    
    //TF1* fm = (TF1*)f->Clone(Form("f%m"));    
    //fm->SetParameter(3, m);
    double x = (TMath::Cos((n-m)*TMath::Pi()/n) + offset)*scale;
    double xs [1] = {x};
    double ps[4]  = {scale,offset,double(n),double(m)};
    cout << "m = " << m << " --> " << f->Eval( x ) << " -- lambda = " << cheb(xs,ps) << endl;    
    f->SetLineColor(m!=9 ? (m+1) : 30);
    leg1->AddEntry(f, Form("m=%d",m), "L");
    //if(m==0) fm->Draw();
    //else fm->Draw("SAME");
    xx[m] = x;
    yy[m] = cheb(xs,ps);
    funcs.emplace_back(f);    
  }

  funcs[0]->Draw();
  TGraph* gr = new TGraph(n+1, xx, yy );
  gr->SetMarkerStyle(3);
  for(auto fm : funcs ) fm->Draw("same");    
  gr->Draw("PSAME");  
  leg1->AddEntry(gr, "Chebyshev points", "P");
  leg1->Draw();
  c->SaveAs(Form("chebyshev_interpolants_n%d.png", n));
 
}
