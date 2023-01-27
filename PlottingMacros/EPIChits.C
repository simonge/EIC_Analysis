void EPIChits(){

  gStyle->SetStatW(0.3);
  gStyle->SetStatX(0.9);
  gStyle->SetStatColor(0);
  gStyle->SetStatBorderSize(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleYOffset(0.7);
  gStyle->SetTitleSize(0.3);

  gStyle->SetPadRightMargin(0.15);
  gStyle->SetPalette(kBird);

  //  TFile* iFile = new TFile("/scratch/EIC/Analysis/qr_out18x275_more_hists.root","READ");
  TFile* iFileQR   = new TFile("/scratch/EIC/ReconOut/QRHists2.root","READ");
  TFile* iFileBrem = new TFile("/scratch/EIC/ReconOut/BremHists2.root","READ");

  TString outName ="/home/simon/Documents/EIC/Plots/Benchmark/EPIC_Hit_Rate_18GeV_pixel.png";

  double Nevents   = 10000000;
  double meanQR    = 0.00374;
  double meanBrems = 11.65;
  double bunchFreq = 10e-9;

  TCanvas* can = new TCanvas("can","can",2500,1500);
  can->Divide(4,2);
  
  TH2F* tag1RawQR     = (TH2F*)iFileQR  ->Get("TaggerAcceptance/Tagger1XY");  
  TH2F* tag1RawBrem   = (TH2F*)iFileBrem->Get("TaggerAcceptance/Tagger1XY");  
  TH2F* tag2RawQR     = (TH2F*)iFileQR  ->Get("TaggerAcceptance/Tagger2XY");
  TH2F* tag2RawBrem   = (TH2F*)iFileBrem->Get("TaggerAcceptance/Tagger2XY");
  tag1RawBrem->Rebin2D();
  tag1RawQR->Rebin2D();
  tag2RawBrem->Rebin2D();
  tag2RawQR->Rebin2D();  

  double pixelArea        = 0.055*0.055;
  double binArea_1        = tag1RawQR->GetXaxis()->GetBinWidth(1)*tag1RawQR->GetYaxis()->GetBinWidth(1);
  double binArea_2        = tag2RawQR->GetXaxis()->GetBinWidth(1)*tag2RawQR->GetYaxis()->GetBinWidth(1);
  double areaScalePixel_1 = pixelArea/binArea_1; //pixel
  double areaScalemm2_1   = 1/binArea_1; //mm2
  double areaScalePixel_2 = pixelArea/binArea_2; //pixel
  double areaScalemm2_2   = 1/binArea_2; //mm2

  cout << tag1RawQR->GetXaxis()->GetBinWidth(1) << " " << tag1RawQR->GetYaxis()->GetBinWidth(1) << endl;
  cout << tag2RawQR->GetXaxis()->GetBinWidth(1) << " " << tag2RawQR->GetYaxis()->GetBinWidth(1) << endl;
  cout << pixelArea << " " << binArea_1 << " " << binArea_2 << endl;

  TH2F* tag1RateBrem = (TH2F*)tag1RawBrem->Clone("tag1rateBrems");
  TH2F* tag1RateQR   = (TH2F*)tag1RawQR->Clone("tag1rateQR");
  TH2F* tag2RateBrem = (TH2F*)tag2RawBrem->Clone("tag2rateBrems");
  TH2F* tag2RateQR   = (TH2F*)tag2RawQR->Clone("tag2rateQR");


  tag1RateBrem->Scale(meanBrems/(Nevents*bunchFreq)); //rate per bin
  tag1RateQR  ->Scale(meanQR/(Nevents*bunchFreq));
  tag2RateBrem->Scale(meanBrems/(Nevents*bunchFreq)); //rate per bin
  tag2RateQR  ->Scale(meanQR/(Nevents*bunchFreq));

  tag1RateBrem->Scale(areaScalePixel_1); //rate per bin
  tag1RateQR  ->Scale(areaScalePixel_1);
  tag2RateBrem->Scale(areaScalePixel_2); //rate per bin
  tag2RateQR  ->Scale(areaScalePixel_2);

//   tag1RateBrem->Scale(areaScalemm2_1); //rate per bin
//   tag1RateQR  ->Scale(areaScalemm2_1);
//   tag2RateBrem->Scale(areaScalemm2_2); //rate per bin
//   tag2RateQR  ->Scale(areaScalemm2_2);

  can->cd(1);
  gPad->SetLogz(1);
  tag1RateBrem->SetTitle("Tagger 1 Bremsstrahlung hit rate [Hz/pixel]; x [mm]; y [mm]");
  tag1RateBrem->Draw("colz");

  can->cd(2);
  gPad->SetLogz(1);
  tag1RateQR->SetTitle("Tagger 1 Quasi-real hit rate [Hz/pixel]; x [mm]; y [mm]");
  tag1RateQR->Draw("colz");

  can->cd(3);
  gPad->SetLogz(1);
  TH2F* tag1RateTotal = (TH2F*)tag1RateBrem->Clone("tag1RateTotal");
  tag1RateTotal->Add(tag1RateQR);
  tag1RateTotal->SetTitle("Tagger 1 Total hit rate [Hz/pixel]; x [mm]; y [mm]");
  tag1RateTotal->Draw("colz");

  can->cd(4);
  gPad->SetLogz(1);
  TH2F* tag1Frac = (TH2F*)tag1RawQR->Clone("tag1QRFrac");
  tag1Frac->Divide(tag1RateQR,tag1RateBrem);
  tag1Frac->SetTitle("Tagger 1 Quasi-real Event Fraction; x [mm]; y [mm]");
  tag1Frac->SetMaximum(1);
  tag1Frac->SetMinimum(0.00001);
  tag1Frac->Draw("colz");

  can->cd(5);
  gPad->SetLogz(1);
  tag2RateBrem->SetTitle("Tagger 2 Bremsstrahlung hit rate [Hz/pixel]; x [mm]; y [mm]");
  tag2RateBrem->Draw("colz");

  can->cd(6);
  gPad->SetLogz(1);
  tag2RateQR->SetTitle("Tagger 2 Quasi-real hit rate [Hz/pixel]; x [mm]; y [mm]");
  tag2RateQR->Draw("colz");

  can->cd(7);
  gPad->SetLogz(1);
  TH2F* tag2RateTotal = ((TH2F*)tag2RateBrem->Clone("tag2RateTotal"));
  tag2RateTotal->Add(tag2RateQR);
  tag2RateTotal->SetTitle("Tagger 2 Total hit rate [Hz/pixel]; x [mm]; y [mm]");
  tag2RateTotal->Draw("colz");

  can->cd(8);
  gPad->SetLogz(1);
  TH2F* tag2Frac = (TH2F*)tag2RawQR->Clone("tag2QRFrac");
  tag2Frac->Divide(tag2RateQR,tag2RateBrem);
  tag2Frac->SetTitle("Tagger 2 Quasi-real Event Fraction; x [mm]; y [mm]");
  tag2Frac->SetMaximum(1);
  tag2Frac->SetMinimum(0.00001);
  tag2Frac->Draw("colz");

  can->SaveAs(outName);

}
