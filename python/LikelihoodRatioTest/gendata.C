{

  // another way to generate data
  TRandom3 rnd;

  double a = 2.0;
  TH1F *hist = new TH1F("hist","",5,0,1);
  for (int i=0;i<nevts;i++){
   	double r = rnd.Exp(1./a);
	hist->Fill(r);
  }
  hist->Draw();

  double sum =0; 
  for (int b=1;b<6;b++){
    //std::cout << "Bin content " << b << " = " << hist->GetBinContent(b) << std::endl; 
    sum+= hist->GetBinContent(b);
  }
  std::cout << " total = " << sum <<std::endl;
  
}
