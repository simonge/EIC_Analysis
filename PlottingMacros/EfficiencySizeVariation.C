#include "TString.h"
#include "TCanvas.h"
#include "ROOT/RDataFrame.hxx"
#include <vector>
#include "TH1.h"
//#include "ROOT.h"


using namespace std;

std::vector<double> pixSize = {0,55,110,220,440,880,1760};

std::vector<TString> fileNames = {"/scratch/EIC/Analysis/tempClusterQR.root"};
//std::vector<TString> fileNames = {"/scratch/EIC/Analysis/clusterFrontWindow.root"};

//std::vector<TString> fileNames = {"/scratch/EIC/Analysis/clusterLarge.root"};

void EfficiencySizeVariation(){

  //TString outNamepng   = "plots/EfficiencySizeVariationFrontWindow.png";
  TString outNamepng   = "plots/EfficiencySizeVariation.png";
  //  TString outNamepng2   = "plots/EfficiencySizeVariation2.png";
  TString outNamepng3   = "plots/EfficiencySizeVariation30.png";

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
  can->Divide(3,1);

  ROOT::RDataFrame df("temp",fileNames[0]);

  int    eBins     = 80;
  int    qBins     = 80;

  double eMin = 0;
  double eMax = 18;
  double qMin = -9;
  double qMax = 0;
  
  auto df1 = df.Filter("iFilter");

  auto generated = df.Histo1D({"EQGen", "Generated Events;Electron Energy [GeV]", eBins, eMin, eMax }, "eE");
  auto forward   = df1.Histo1D({"EQForward", "Forward Generated Events;Electron Energy [GeV]", eBins, eMin, eMax }, "eE");
  auto accepted  = df1.Filter("NTracks>=1").Histo1D({"EQAccepted", "AcceptedEvents Events;Electron Energy [GeV]", eBins, eMin, eMax }, "eE");
  auto Chip1  = df1.Filter("Any(fit_minY>-256*1&&fit_maxY<256*1)").Histo1D({"EQAccepted1", "AcceptedEvents Events;Electron Energy [GeV]", eBins, eMin, eMax }, "eE");
  auto Chip2  = df1.Filter("Any(fit_minY>-256*2&&fit_maxY<256*2)").Histo1D({"EQAccepted2", "AcceptedEvents Events;Electron Energy [GeV]", eBins, eMin, eMax }, "eE");
  auto Chip3  = df1.Filter("Any(fit_minY>-256*3&&fit_maxY<256*3)").Histo1D({"EQAccepted3", "AcceptedEvents Events;Electron Energy [GeV]", eBins, eMin, eMax }, "eE");
  auto Chip4  = df1.Filter("Any(fit_minY>-256*4&&fit_maxY<256*4)").Histo1D({"EQAccepted4", "AcceptedEvents Events;Electron Energy [GeV]", eBins, eMin, eMax }, "eE");
  auto Chip5  = df1.Filter("Any(fit_minY>-256*5&&fit_maxY<256*5)").Histo1D({"EQAccepted5", "AcceptedEvents Events;Electron Energy [GeV]", eBins, eMin, eMax }, "eE");
  auto Chip6  = df1.Filter("Any(fit_minY>-256*6&&fit_maxY<256*6)").Histo1D({"EQAccepted6", "AcceptedEvents Events;Electron Energy [GeV]", eBins, eMin, eMax }, "eE");

  can->cd(1);
  //  generated->Draw("hist");
  int GenCounts = generated->GetEntries();
  int ForwardCounts = forward->GetEntries();
  int acceptedCounts = accepted->GetEntries();
  int Chip1Counts = Chip1->GetEntries();
  int Chip2Counts = Chip2->GetEntries();
  int Chip3Counts = Chip3->GetEntries();
  int Chip4Counts = Chip4->GetEntries();
  int Chip5Counts = Chip5->GetEntries();
  int Chip6Counts = Chip6->GetEntries();
  
  TH1D* forw = (TH1D*)generated->Clone("For");
  forw->Divide(forward.GetPtr(),generated.GetPtr());
  forw->SetMinimum(0);
  forw->SetMaximum(1);
  forw->SetTitle("Energy Dependent Efficiency");
  forw->Draw("hist");
  TH1D* acc = (TH1D*)generated->Clone("Acc");
  acc->Divide(accepted.GetPtr(),generated.GetPtr());
  acc->SetLineColor(kViolet);
  acc->Draw("hist same");
  TH1D* C1 = (TH1D*)generated->Clone("Acc");
  C1->Divide(Chip1.GetPtr(),generated.GetPtr());
  C1->SetLineColor(kRed);
  C1->Draw("hist same");;
  TH1D* C2 = (TH1D*)generated->Clone("Acc");
  C2->Divide(Chip2.GetPtr(),generated.GetPtr());
  C2->SetLineColor(kBlue);
  C2->Draw("hist same");;
  TH1D* C3 = (TH1D*)generated->Clone("Acc");
  C3->Divide(Chip3.GetPtr(),generated.GetPtr());
  C3->SetLineColor(kGreen);
  C3->Draw("hist same");;
  TH1D* C4 = (TH1D*)generated->Clone("Acc");
  C4->Divide(Chip4.GetPtr(),generated.GetPtr());
  C4->SetLineColor(kYellow);
  C4->Draw("hist same");;
  TH1D* C5 = (TH1D*)generated->Clone("Acc");
  C5->Divide(Chip5.GetPtr(),generated.GetPtr());
  C5->SetLineColor(kOrange);
  C5->Draw("hist same");;
  TH1D* C6 = (TH1D*)generated->Clone("Acc");
  C6->Divide(Chip6.GetPtr(),generated.GetPtr());
  C6->SetLineColor(kGray);
  C6->Draw("hist same");
  
  
  auto legend = new TLegend(0.2,0.8,0.45,0.88);
  legend->AddEntry(forw,"Forward (theta<10 mrad)","l");
  legend->AddEntry(acc,"Simulation Accepted","l");
//   legend->AddEntry(C6,"16.9 cm height","l");
//   legend->AddEntry(C5,"14.1 cm height","l");
//   legend->AddEntry(C4,"11.3 cm height","l");
//   legend->AddEntry(C3,"8.4 cm height","l");
//   legend->AddEntry(C2,"5.6 cm height","l");
//   legend->AddEntry(C1,"2.8 cm height","l");
  legend->Draw();
  
  TH1D* fracEff1 = (TH1D*)C1->Clone("FracEff1");
  fracEff1->Divide(fracEff1,acc);
  TH1D* fracEff2 = (TH1D*)C2->Clone("FracEff2");
  fracEff2->Divide(fracEff2,acc);
  TH1D* fracEff3 = (TH1D*)C3->Clone("FracEff3");
  fracEff3->Divide(fracEff3,acc);
  TH1D* fracEff4 = (TH1D*)C4->Clone("FracEff4");
  fracEff4->Divide(fracEff4,acc);
  TH1D* fracEff5 = (TH1D*)C5->Clone("FracEff5");
  fracEff5->Divide(fracEff5,acc);
  TH1D* fracEff6 = (TH1D*)C6->Clone("FracEff6");
  fracEff6->Divide(fracEff6,acc);

  can->cd(2);
  fracEff1->SetLineColor(kRed);
  fracEff1->Draw("hist");
  fracEff2->SetLineColor(kBlue);
  fracEff2->Draw("same hist");
  fracEff3->SetLineColor(kGreen);
  fracEff3->Draw("same hist");
  fracEff4->SetLineColor(kYellow);
  fracEff4->Draw("same hist");
  fracEff5->SetLineColor(kOrange);
  fracEff5->Draw("same hist");
  fracEff6->SetLineColor(kGray);
  fracEff6->Draw("same hist");

  ROOT::VecOps::RVec<double> xpoints = {-2,-1,0,1,2,3,4,5,6};
  ROOT::VecOps::RVec<double> ypoints = {(double)GenCounts,(double)ForwardCounts,(double)acceptedCounts,(double)Chip1Counts,(double)Chip2Counts,(double)Chip3Counts,(double)Chip4Counts,(double)Chip5Counts,(double)Chip6Counts};
  ypoints/=ForwardCounts;
  //  ypoints/=GenCounts;

  TGraph eff(xpoints.size(), &xpoints[0], &ypoints[0]);
  can->cd(3);
  eff.SetMinimum(0);
  eff.SetMarkerStyle(2);
  eff.DrawClone("AP");    


  
  can->SaveAs(outNamepng);


  TCanvas* can2 = new TCanvas("can2","can2",1400,1400);
  forw->Draw("hist");
  acc->Draw("hist same");
//   C1->Draw("hist same");
//   C2->Draw("hist same");
//   C3->Draw("hist same");
//   C4->Draw("hist same");
//   C5->Draw("hist same"); 
//   C6->Draw("hist same");
  legend->Draw();
  
  can2->SaveAs(outNamepng3);
  


}
