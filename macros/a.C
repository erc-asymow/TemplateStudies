{

  auto cheb = [](double x, double scale, double offset, unsigned int n, unsigned int m){
    double den = 0.;
    double num = 0.;
    for(unsigned int i = 0; i <= n ; i++){
      int sign = i%2==0 ? +1 :-1;
      double xj = (TMath::Cos((n-i)*TMath::Pi()/n) + offset)*scale;
      if(x==xj) return 1.0;// protect from nan      
      double val = sign/(x-xj);
      if(i==0 || i==n) val *= 0.5;
      den += val;
      if(i==m) num = val;
    }
    //std::cout << x << "==>" <<  num << "," << den << std::endl;                                             
    return num/den;
  };

  int deg = 12;
  unsigned int mid_deg = deg/2;
  double max_y = 3.0;
  for(unsigned int l = 0; l<=mid_deg; l++){
    for(int iy = 0; iy<10; iy++){
      double y = -max_y + iy*(2*max_y)/10.;
      double cheb_a = cheb( y, max_y, 0.0, deg, l) + cheb(y, max_y, 0.0, deg, deg-l) ;
      double cheb_b = cheb(-y, max_y, 0.0, deg, l) + cheb(-y, max_y, 0.0, deg, deg-l) ;
      cout << "y=" << y  << ", l=" << l << " -> " << cheb_a - cheb_b << endl;
    }
  }

}
