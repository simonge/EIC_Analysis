void EPICAcceptance(){

  gStyle->SetStatW(0.3);
  gStyle->SetStatX(0.9);
  gStyle->SetStatColor(0);
  gStyle->SetStatBorderSize(1);
  gStyle->SetOptStat(10);
  gStyle->SetOptTitle(1);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleYOffset(0.7);
  gStyle->SetTitleSize(0.3);

  gStyle->SetPadRightMargin(0.15);
  gStyle->SetPalette(kBird);

  //  TFile* iFile = new TFile("/scratch/EIC/Analysis/qr_out18x275_more_hists.root","READ");
  TFile* iFile = new TFile("/scratch/EIC/Results/tempPlots.root","READ");

  TString outName ="/home/simon/Documents/EIC/Plots/Benchmark/EPICAcceptance_Brem_18GeV.png";
  
  TCanvas* can = new TCanvas("can","can",2000,1500);
  can->Divide(3,2);
  
  can->cd(1);
  gPad->SetLogz(1);
  TH2F* raw = (TH2F*)iFile->Get("TaggerAcceptance/EQraw");  
  raw->SetTitle("Generated Events; Electron Energy [GeV]; log_{10}(Q^{2})");
  raw->Draw("colz");
  
  can->cd(2);
  gPad->SetLogz(1);
  TH2F* tag1 = (TH2F*)iFile->Get("TaggerAcceptance/EQtag1");
  tag1->SetTitle("Tagger 1 Hits; Electron Energy [GeV]; log_{10}(Q^{2})");
  tag1->Draw("colz");
  
  can->cd(3);
  gPad->SetLogz(1);
  TH2F* tag2 = (TH2F*)iFile->Get("TaggerAcceptance/EQtag2");
  tag2->SetTitle("Tagger 2 Hits; Electron Energy [GeV]; log_{10}(Q^{2})");
  tag2->Draw("colz");

  can->cd(4);
  TH2F* totalA = (TH2F*)iFile->Get("TaggerAcceptance/EQaccept");
  totalA->SetTitle("Total Acceptance; Electron Energy [GeV]; log_{10}(Q^{2})");
  totalA->Draw("colz");

  can->cd(5);
  TH2F* tag1A = (TH2F*)iFile->Get("TaggerAcceptance/EQtag1accept");
  tag1A->SetTitle("Tagger 1 Acceptance; Electron Energy [GeV]; log_{10}(Q^{2})");
  tag1A->Draw("colz");

  can->cd(6);
  TH2F* tag2A = (TH2F*)iFile->Get("TaggerAcceptance/EQtag2accept");
  tag2A->SetTitle("Tagger 2 Acceptance; Electron Energy [GeV]; log_{10}(Q^{2})");
  tag2A->Draw("colz");

  can->SaveAs(outName);

}
