#include "TString.h"
#include "TCanvas.h"
#include "ROOT/RDataFrame.hxx"
#include <vector>
#include "TH1.h"
//#include "ROOT.h"


using namespace std;

std::vector<double> pixSize = {0,55,110,220,440,880,1760};

std::vector<TString> fileNames = {"/scratch/EIC/Analysis/temp.root"};


void Acceptances(){

  TString outName      = "AcceptancesPipe.pdf";
  TString outNamepng   = "AcceptancesPipe.png";
  TString outName2     = "AcceptancesBPipe.pdf";
  TString outNamepng2  = "AcceptancesBPipe.png";
  TString outName3     = "AcceptancesCPipe.pdf";
  TString outNamepng3  = "AcceptancesCPipe.png";

  gStyle->SetStatW(0.3);
  gStyle->SetStatColor(0);
  gStyle->SetStatBorderSize(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  gStyle->SetTitleStyle(0);
  //  gStyle->SetTitleSize(0.2);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.13);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetLabelSize(0.035,"X");
  gStyle->SetLabelSize(0.035,"Y");
  gStyle->SetLabelSize(0.035,"Z");
  gStyle->SetTitleSize(0.06,"X");
  gStyle->SetTitleSize(0.06,"Y");
  gStyle->SetTitleSize(0.06,"Z");

  gStyle->SetTitleOffset(1.0,"Y");

  TCanvas* can = new TCanvas("can","can",2200,1400);
  can->Divide(2,2);


  ROOT::RDataFrame df("temp",fileNames[0]);

  int    eBins     = 80;
  int    qBins     = 80;

  double eMin = 0;
  double eMax = 18;
  double qMin = -9;
  double qMax = 0;

  auto generated = df.Histo2D({"EQGen", "Generated Events;Electron Energy [GeV];Log(Q^{2})", eBins, eMin, eMax, qBins, qMin, qMax }, "eE","logQ2_3");
  auto thetaCut  = df.Filter("iFilter").Histo2D({"EQTheta", "Events - Theta < 10 mrad;Electron Energy [GeV];Log(Q^{2})", eBins, eMin, eMax, qBins, qMin, qMax }, "eE","logQ2_3");
  auto acceptedCut  = df.Filter("iFilter&&(Tag1_4||Tag2_4)").Histo2D({"EQaccepted", "Detector Hits;Electron Energy [GeV];Log(Q^{2})", eBins, eMin, eMax, qBins, qMin, qMax }, "eE","logQ2_3");

  auto acceptedCut1  = df.Filter("iFilter&&Tag1_4").Histo2D({"EQaccepted", "Detector Hits Tagger 1;Electron Energy [GeV];Log(Q^{2})", eBins, eMin, eMax, qBins, qMin, qMax }, "eE","logQ2_3");
  auto acceptedCut2  = df.Filter("iFilter&&Tag2_4").Histo2D({"EQaccepted", "Detector Hits Tagger 2;Electron Energy [GeV];Log(Q^{2})", eBins, eMin, eMax, qBins, qMin, qMax }, "eE","logQ2_3");
  

  auto acceptedPhiRes  = df.Filter("iFilter&&(Tag1_4||Tag2_4)&&(TMath::Pi()-scatteredElectron.Theta())>0.001").Histo2D({"EQacceptedPhi", "Phi Resolvable;Electron Energy [GeV];Log(Q^{2})", eBins, eMin, eMax, qBins, qMin, qMax }, "eE","logQ2_3");

  can->cd(1);
  gPad->SetLogz();
  generated->Draw("colz");
  can->cd(2);
  gPad->SetLogz();
  thetaCut->Draw("colz");
  cout << thetaCut->GetEntries() << endl;
  can->cd(3);
  gPad->SetLogz();
  acceptedCut->Draw("colz");
  cout << acceptedCut->GetEntries() << endl;

  TH2D* acceptance = (TH2D*)thetaCut->Clone("TaggerAcceptance");
  acceptance->Divide(acceptedCut.GetPtr(),thetaCut.GetPtr());
  acceptance->SetTitle("Acceptance");

  can->cd(4);
  acceptance->Draw("colz");
  
  can->SaveAs(outName);
  can->SaveAs(outNamepng);


  TCanvas* can2 = new TCanvas("can2","can2",2200,800);
  can2->Divide(3,1);
  
  TH2D* acceptance1 = (TH2D*)thetaCut->Clone("TaggerAcceptance1");
  acceptance1->Divide(acceptedCut1.GetPtr(),thetaCut.GetPtr());
  acceptance1->SetTitle("Tagger 1 Acceptance");

  TH2D* acceptance2 = (TH2D*)thetaCut->Clone("TaggerAcceptance2");
  acceptance2->Divide(acceptedCut2.GetPtr(),thetaCut.GetPtr());
  acceptance2->SetTitle("Tagger 2 Acceptance");

  TH2D* acceptanceSum = (TH2D*)thetaCut->Clone("TaggerAcceptanceSum");
  acceptanceSum->Add(acceptance1,acceptance2);
  acceptanceSum->SetTitle("Acceptance including double counting");
  
  can2->cd(1);
  acceptance1->Draw("colz");
  can2->cd(2);
  acceptance2->Draw("colz");
  can2->cd(3);
  acceptanceSum->Draw("colz");
  
  can2->SaveAs(outName2);
  can2->SaveAs(outNamepng2);

  gStyle->SetOptStat(0001);

  TCanvas* can3 = new TCanvas("can3","can3",2400,800);
  can3->Divide(4,1);
  
  can3->cd(1);
  gPad->SetLogz();
  generated->SetStats();
  generated->Draw("colz");

  can3->cd(2);
  gPad->SetLogz();
  thetaCut->SetStats();
  thetaCut->Draw("colz");

  can3->cd(3);
  gPad->SetLogz();
  acceptedCut->SetStats();
  acceptedCut->Draw("colz");

  can3->cd(4);
  gPad->SetLogz();
  acceptedPhiRes->SetStats();
  acceptedPhiRes->Draw("colz");

  
  can3->SaveAs(outName3);
  can3->SaveAs(outNamepng3);


}
