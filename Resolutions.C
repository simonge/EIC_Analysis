#include "TString.h"
#include "TCanvas.h"
#include "ROOT/RDataFrame.hxx"
#include <vector>
#include "TH1.h"
//#include "ROOT.h"


using namespace std;

std::vector<double> pixSize = {0,55,110,220,440,880,1760};

std::string tag = "Decor-NormTE-4Layers-Bigger-4STEP5";

std::vector<TString> fileNames = {"/scratch/EIC/Results/ML-Out/"+tag+"_real_ETP.root","/home/simon/Analysis/EIC_Analysis/Reg_cell-"+tag+".root"};
//std::vector<TString> fileNames = {"/scratch/EIC/Results/ML-Out/4layer_real_ETP.root","/home/simon/Analysis/EIC_Analysis/Reg_point-RealHits4layer.root"};
//std::vector<TString> fileNames = {"/scratch/EIC/Results/ML-Out/4layer_real_ETP.root","/home/simon/Analysis/EIC_Analysis/Reg_cell-RealHits4layer.root"};
//std::vector<TString> fileNames = {"/scratch/EIC/Results/ML-Out/test_cell_ETP.root"};


void Resolutions(){

//   TString outName  = "PTestResolutionsA.pdf";
//   TString outName2 = "PTestResolutionsB.pdf";
//   TString outName  = "ResolutionsA-New.pdf";
//   TString outName2 = "ResolutionsB-New.pdf";
//   TString outNamepng  = "ResolutionsA-JustPipe.png";
//   TString outNamepng2 = "ResolutionsB-JustPipe.png";
  TString outNamepng  = "ETPResolutionsA-"+tag+".png";
  TString outNamepng2 = "ETPResolutionsB-"+tag+".png";

  gStyle->SetStatW(0.3);
  gStyle->SetStatColor(0);
  gStyle->SetStatBorderSize(1);
  gStyle->SetOptStat(0);
  //  gStyle->SetOptStat(1100);
  gStyle->SetOptTitle(1);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.13);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetLabelSize(0.05,"X");
  gStyle->SetLabelSize(0.05,"Y");
  gStyle->SetLabelSize(0.042,"Z");
  gStyle->SetTitleSize(0.06,"X");
  gStyle->SetTitleSize(0.06,"Y");
  gStyle->SetTitleSize(0.06,"Z");

  gStyle->SetTitleOffset(1.0,"Y");

  TCanvas* can = new TCanvas("can","can",2200,1400);
  can->Divide(3,2);

  TCanvas* can2 = new TCanvas("can2","can2",2200,1400);
  can2->Divide(3,2);

//   ROOT::RDataFrame dfcheet("dataset/TestTree",fileNames[0]);

//   //Cheating
//   auto dfcheet_2 = dfcheet.Define("CalPhi","atan2(DNN_CPU.sin_phiV_,DNN_CPU.cos_phiV_)")
//     .Define("Phi","atan2(sin_phiV_,cos_phiV_)")
//     .Define("PhiRes","(Phi-CalPhi)")
//     .Define("ThetaGen","thetaE*1000")
//     .Define("ThetaRecon","DNN_CPU.thetaE*1000")
//     .Define("ThetaRes","(ThetaGen-ThetaRecon)")
//     .Define("ePred","DNN_CPU.eE")
//     .Define("ERes","100*(eE-DNN_CPU.eE)/eE");
  ROOT::RDataFrame dfcheet("predict",fileNames[1]);

  auto dfcheet_2 = dfcheet.Define("CalPhi","pPred")
    .Define("Phi","phiV")
    .Define("PhiRes","(Phi-CalPhi)")
    .Define("ThetaGen","thetaE*1000")
    .Define("ThetaRecon","tPred*1000")
    .Define("ThetaRes","(ThetaGen-ThetaRecon)")
    .Define("ERes","100*(eE-ePred)/eE");

  ROOT::RDataFrame df55("predict",fileNames[1]);

  //Not Cheating
  auto df55_2 = df55.Define("CalPhi","pPred")
    .Define("Phi","phiV")
    .Define("PhiRes","(Phi-CalPhi)")
    .Define("ThetaGen","thetaE*1000")
    .Define("ThetaRecon","tPred*1000")
    .Define("ThetaRes","(ThetaGen-ThetaRecon)")
    .Define("ERes","100*(eE-ePred)/eE");

  int    eBins     = 80;
  int    phiBins   = 80;
  int    thetaBins = 80;
  int    resBins   = 50;

  double eMin = 6;
  double eMax = 18;
  double phiMin = -TMath::Pi();
  double phiMax = TMath::Pi();
  double thetaMin = 0.0;
  double thetaMax = 0.011*1000;

  double PResRange = 1.5;
  double TResRange = 0.002*1000;
  double EResRange = 0.05*100;

  auto Phi1D    = dfcheet_2.Filter("ThetaGen>1").Histo1D({"PRes", ";#phi gen-recon [rad]", resBins, -PResRange, PResRange }, "PhiRes");
  auto Theta1D  = dfcheet_2.Histo1D({"TRes", ";#theta gen-recon [mrad]", resBins, -TResRange, TResRange }, "ThetaRes");
  auto Energy1D = dfcheet_2.Histo1D({"ERes", ";Energy (gen-recon)/gen [%]", resBins, -EResRange, EResRange }, "ERes");

  auto Phi1D55    = df55_2.Filter("ThetaGen>1").Histo1D({"PRes55", ";#phi gen-recon [rad]", resBins, -PResRange, PResRange }, "PhiRes");
  auto Theta1D55  = df55_2.Histo1D({"TRes55", ";#theta gen-recon [mrad]", resBins, -TResRange, TResRange }, "ThetaRes");
  auto Energy1D55 = df55_2.Histo1D({"ERes55", ";Energy (gen-recon)/gen [%]", resBins, -EResRange, EResRange }, "ERes");

  auto Phi2D    = dfcheet_2.Filter("ThetaGen>1").Histo2D({"Phi",    "Phi Reconstruction (Theta > 1 mrad);#phi gen [rad];#phi recon [rad]",      phiBins,   phiMin,   phiMax,   phiBins,   phiMin,   phiMax,  }, "Phi",   "CalPhi");
  auto Theta2D  = dfcheet_2.Histo2D({"Theta",  "Theta Reconstruction;#theta gen [mrad] ;#theta recon [mrad]", thetaBins, thetaMin, thetaMax, thetaBins, thetaMin, thetaMax }, "ThetaGen","ThetaRecon");
  auto Energy2D = dfcheet_2.Histo2D({"Energy", "Energy Reconstruction;Energy gen [GeV] ;Energy recon [GeV]", eBins,     eMin,     eMax,     eBins,     eMin,     eMax     }, "eE",    "ePred");

  auto Phi2D2    = dfcheet_2.Filter("ThetaGen>1").Histo2D({"PRes", ";#phi gen-recon [rad];Electron energy [GeV] ",   resBins, -PResRange, PResRange, eBins, eMin, eMax },"PhiRes",  "eE");
  auto Theta2D2  = dfcheet_2.Histo2D({"TRes", ";#theta gen-recon [mrad];Electron energy [GeV] ", resBins, -TResRange, TResRange, eBins, eMin, eMax },"ThetaRes","eE");
  auto Energy2D2 = dfcheet_2.Histo2D({"ERes", ";Energy (gen-recon)/gen [%];Electron energy [GeV] ",       resBins, -EResRange, EResRange, eBins, eMin, eMax },"ERes",    "eE");

  auto Phi3D    = dfcheet_2.Histo3D({"PRes", ";Electron energy [GeV] ;#theta [mrad]", eBins, eMin, eMax, thetaBins, thetaMin, thetaMax, resBins, -PResRange, PResRange }, "eE","ThetaGen","PhiRes");
  auto Theta3D  = dfcheet_2.Histo3D({"TRes", ";Electron energy [GeV] ;#theta [mrad]", eBins, eMin, eMax, thetaBins, thetaMin, thetaMax, resBins, -TResRange, TResRange }, "eE","ThetaGen","ThetaRes");
  auto Energy3D = dfcheet_2.Histo3D({"ERes", ";Electron energy [GeV] ;#theta [mrad]", eBins, eMin, eMax, thetaBins, thetaMin, thetaMax, resBins, -EResRange, EResRange }, "eE","ThetaGen","ERes");

  (Phi3D.GetPtr())->FitSlicesZ();
  (Theta3D.GetPtr())->FitSlicesZ();
  (Energy3D.GetPtr())->FitSlicesZ();

  TH2D* phiSigma    = (TH2D*)(gDirectory->Get("PRes_2"));
  TH2D* thetaSigma  = (TH2D*)(gDirectory->Get("TRes_2"));
  TH2D* energySigma = (TH2D*)(gDirectory->Get("ERes_2"));
//   TH2D* phiSigma    = (TH2D*)(gDirectory->Get("PRes_chi2"));
//   TH2D* thetaSigma  = (TH2D*)(gDirectory->Get("TRes_chi2"));
//   TH2D* energySigma = (TH2D*)(gDirectory->Get("ERes_chi2"));

  can->cd(1);
  energySigma->SetTitle("Energy Resolutions [%]");
  energySigma->SetStats(0);
  energySigma->SetMaximum(EResRange*2);
  energySigma->Draw("colz");
  can->cd(2);
  thetaSigma->SetTitle("Theta Resolutions [mrad]");
  thetaSigma->SetStats(0);
  thetaSigma->SetMaximum(TResRange*2);
  thetaSigma->Draw("colz");
  can->cd(3);  
  phiSigma->SetTitle("Phi Resolutions [rad]");
  phiSigma->SetStats(0);
  phiSigma->SetMaximum(0.3);
  phiSigma->Draw("colz");
   
  can->cd(4);
  gPad->SetLogz();
  Energy2D2.GetPtr()->Draw("colz");
  Energy2D2->SetStats(0);
  can->cd(5);
  gPad->SetLogz();
  Theta2D2.GetPtr()->Draw("colz");
  Theta2D2->SetStats(0);
  can->cd(6);
  gPad->SetLogz();
  Phi2D2.GetPtr()->Draw("colz");
  Phi2D2->SetStats(0);
  
//   can->SaveAs(outName2);
  can->SaveAs(outNamepng2);

  can->cd(1);
  gPad->SetLogz();
  Energy2D.GetPtr()->Draw("colz");
  Energy2D->SetStats(0);
  can->cd(2);
  gPad->SetLogz();
  Theta2D.GetPtr()->Draw("colz");
  Theta2D->SetStats(0);
  can->cd(3);
  gPad->SetLogz();
  Phi2D.GetPtr()->Draw("colz");
  Phi2D->SetStats(0);
  
  TH1D* ERes = (TH1D*)(Energy1D.GetPtr());
  TH1D* TRes = (TH1D*)(Theta1D.GetPtr());
  TH1D* PRes = (TH1D*)(Phi1D.GetPtr());

  TH1D* ERes55 = (TH1D*)(Energy1D55.GetPtr());
  TH1D* TRes55 = (TH1D*)(Theta1D55.GetPtr());
  TH1D* PRes55 = (TH1D*)(Phi1D55.GetPtr());

  can->cd(4);
  ERes->Draw();
  ERes->Fit("gaus");
//   ERes55->Draw("same");
//   ERes55->Fit("gaus");
  ERes->SetStats(0);
  can->cd(5);
  TRes->Draw();
  TRes->Fit("gaus");
//   TRes55->Draw("same");
//   TRes55->Fit("gaus");
  TRes->SetStats(0);
  can->cd(6);
  PRes->Draw();
  PRes->Fit("gaus");
//   PRes55->Draw("same");
//   PRes55->Fit("gaus");
  PRes->SetStats(0);
  
//   can->SaveAs(outName);
  can->SaveAs(outNamepng);


}
