#include "TString.h"
#include "TCanvas.h"
#include "ROOT/RDataFrame.hxx"
#include <vector>
#include "TH1.h"
//#include "ROOT.h"


using namespace std;

std::vector<double> pixSize = {0,55,110,220,440,880,1760};

std::vector<TString> fileNames = {"/scratch/EIC/Analysis/tempClusterQR.root"};

//std::vector<TString> fileNames = {"/scratch/EIC/Analysis/clusterLarge.root"};

void AcceptancesSizeVariation(){

  TString outNamepng   = "AcceptancesSizeVariation.png";
  TString outNamepng2   = "AcceptancesSizeVariation2.png";

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
  can->Divide(5,2);
  TCanvas* can2 = new TCanvas("can2","can2",2200,1400);
  can2->Divide(5,2);

  ROOT::RDataFrame df("temp",fileNames[0]);

  int    eBins     = 80;
  int    qBins     = 80;

  double eMin = 0;
  double eMax = 18;
  double qMin = -9;
  double qMax = 0;
  
  auto df1 = df.Filter("iFilter");

  auto generated = df.Histo2D({"EQGen", "Generated Events;Electron Energy [GeV];Log(Q^{2})", eBins, eMin, eMax, qBins, qMin, qMax }, "eE","logQ2_3");
  auto forward = df1.Histo2D({"EQForward", "Forward Generated Events;Electron Energy [GeV];Log(Q^{2})", eBins, eMin, eMax, qBins, qMin, qMax }, "eE","logQ2_3");
  auto accepted  = df1.Filter("NTracks>=1").Histo2D({"EQAccepted", "AcceptedEvents Events;Electron Energy [GeV];Log(Q^{2})", eBins, eMin, eMax, qBins, qMin, qMax }, "eE","logQ2_3");
  auto Chip1  = df1.Filter("Any(fit_minY>-256&&fit_maxY<256)").Histo2D({"EQAccepted1", "AcceptedEvents Events;Electron Energy [GeV];Log(Q^{2})", eBins, eMin, eMax, qBins, qMin, qMax }, "eE","logQ2_3");
  auto Chip2  = df1.Filter("Any(fit_minY>-512&&fit_maxY<512)").Histo2D({"EQAccepted2", "AcceptedEvents Events;Electron Energy [GeV];Log(Q^{2})", eBins, eMin, eMax, qBins, qMin, qMax }, "eE","logQ2_3");
  auto Chip3  = df1.Filter("Any(fit_minY>-256*3&&fit_maxY<256*3)").Histo2D({"EQAccepted3", "AcceptedEvents Events;Electron Energy [GeV];Log(Q^{2})", eBins, eMin, eMax, qBins, qMin, qMax }, "eE","logQ2_3");
  auto Chip4  = df1.Filter("Any(fit_minY>-256*4&&fit_maxY<256*4)").Histo2D({"EQAccepted4", "AcceptedEvents Events;Electron Energy [GeV];Log(Q^{2})", eBins, eMin, eMax, qBins, qMin, qMax }, "eE","logQ2_3");
  auto Chip5  = df1.Filter("Any(fit_minY>-256*5&&fit_maxY<256*5)").Histo2D({"EQAccepted5", "AcceptedEvents Events;Electron Energy [GeV];Log(Q^{2})", eBins, eMin, eMax, qBins, qMin, qMax }, "eE","logQ2_3");
  auto Chip6  = df1.Filter("Any(fit_minY>-256*6&&fit_maxY<256*6)").Histo2D({"EQAccepted6", "AcceptedEvents Events;Electron Energy [GeV];Log(Q^{2})", eBins, eMin, eMax, qBins, qMin, qMax }, "eE","logQ2_3");

  can->cd(1);
  gPad->SetLogz();
  generated->Draw("colz");
  int GenCounts = generated->GetEntries();
  can->cd(2);
  gPad->SetLogz();
  forward->Draw("colz");
  int ForwardCounts = forward->GetEntries();
  can->cd(3);
  gPad->SetLogz();
  accepted->Draw("colz");
  int acceptedCounts = accepted->GetEntries();
  can->cd(4);
  gPad->SetLogz();
  Chip1->Draw("colz");
  int Chip1Counts = Chip1->GetEntries();
  can->cd(4);
  gPad->SetLogz();
  Chip2->Draw("colz");
  int Chip2Counts = Chip2->GetEntries();
  can->cd(5);
  gPad->SetLogz();
  Chip3->Draw("colz");
  int Chip3Counts = Chip3->GetEntries();
  can->cd(6);
  gPad->SetLogz();
  Chip3->Draw("colz");
  int Chip4Counts = Chip4->GetEntries();
  can->cd(7);
  gPad->SetLogz();
  Chip4->Draw("colz");
  int Chip5Counts = Chip5->GetEntries();
  can->cd(8);
  gPad->SetLogz();
  Chip6->Draw("colz");
  int Chip6Counts = Chip6->GetEntries();
  
  ROOT::VecOps::RVec<double> xpoints = {-2,-1,0,1,2,3,4,5,6};
  ROOT::VecOps::RVec<double> ypoints = {(double)GenCounts,(double)ForwardCounts,(double)acceptedCounts,(double)Chip1Counts,(double)Chip2Counts,(double)Chip3Counts,(double)Chip4Counts,(double)Chip5Counts,(double)Chip6Counts};
  ypoints/=ForwardCounts;

  TGraph eff(xpoints.size(), &xpoints[0], &ypoints[0]);
  can->cd(9);
  eff.SetMinimum(0);
  eff.SetMarkerStyle(2);
  eff.Draw("AP");    
  
  can->SaveAs(outNamepng);


  TH2D* acc = (TH2D*)generated->Clone("TaggerAcceptance0");
  acc->Divide(accepted.GetPtr(),accepted.GetPtr());
  acc->SetTitle("Acc");
  acc->Draw("colz"); 
  double ForwardInt = acc->Integral();

  can2->cd(1);
  TH2D* acceptance = (TH2D*)generated->Clone("TaggerAcceptanceA");
  acceptance->Divide(accepted.GetPtr(),forward.GetPtr());
  acceptance->SetTitle("Acceptance");
  acceptance->Draw("colz");  
  double AcceptedInt = acceptance->Integral();
  can2->cd(2);
  TH2D* acceptance1 = (TH2D*)generated->Clone("TaggerAcceptance1");
  acceptance1->Divide(Chip1.GetPtr(),forward.GetPtr());
  acceptance1->SetTitle("Acceptance1");
  acceptance1->Draw("colz");
  double Chip1Int = acceptance1->Integral();
  can2->cd(3);
  TH2D* acceptance2 = (TH2D*)generated->Clone("TaggerAcceptance2");
  acceptance2->Divide(Chip2.GetPtr(),forward.GetPtr());
  acceptance2->SetTitle("Acceptance2");
  acceptance2->Draw("colz");
  double Chip2Int = acceptance2->Integral();
  can2->cd(4);
  TH2D* acceptance3 = (TH2D*)generated->Clone("TaggerAcceptance3");
  acceptance3->Divide(Chip3.GetPtr(),forward.GetPtr());
  acceptance3->SetTitle("Acceptance3");
  acceptance3->Draw("colz");
  double Chip3Int = acceptance3->Integral();
  can2->cd(5);
  TH2D* acceptance4 = (TH2D*)generated->Clone("TaggerAcceptance4");
  acceptance4->Divide(Chip4.GetPtr(),forward.GetPtr());
  acceptance4->SetTitle("Acceptance4");
  acceptance4->Draw("colz");
  double Chip4Int = acceptance4->Integral();
  can2->cd(6);
  TH2D* acceptance5 = (TH2D*)generated->Clone("TaggerAcceptance5");
  acceptance5->Divide(Chip5.GetPtr(),forward.GetPtr());
  acceptance5->SetTitle("Acceptance5");
  acceptance5->Draw("colz");
  double Chip5Int = acceptance5->Integral();
  can2->cd(7);
  TH2D* acceptance6 = (TH2D*)generated->Clone("TaggerAcceptance6");
  acceptance6->Divide(Chip6.GetPtr(),forward.GetPtr());
  acceptance6->SetTitle("Acceptance6");
  double Chip6Int = acceptance6->Integral();
  acceptance6->Draw("colz");

  ROOT::VecOps::RVec<double> xpoints2 = {-1,0,1,2,3,4,5,6};
  ROOT::VecOps::RVec<double> ypoints2 = {ForwardInt,AcceptedInt,Chip1Int,Chip2Int,Chip3Int,Chip4Int,Chip5Int,Chip6Int};
  ypoints2/=ForwardInt;

  TGraph eff2(xpoints2.size(), &xpoints2[0], &ypoints2[0]);
  can2->cd(9);
  eff2.SetMinimum(0);
  eff2.SetMarkerStyle(2);
  eff2.Draw("AP");    


  can2->SaveAs(outNamepng2);

}
