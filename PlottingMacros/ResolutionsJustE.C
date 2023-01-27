#include "TString.h"
#include "TCanvas.h"
#include "ROOT/RDataFrame.hxx"
#include <vector>
#include "TH1.h"
//#include "ROOT.h"


using namespace std;

std::vector<double> pixSize = {0,55,110,220,440,880,1760};

std::string tag = "Decor-NormTE-4Layers-Bigger-4STEP5";
std::string tag2 = "cell-RealHitsDecor-NormE-4Layers-Bigger10-2STEP4FrontWindow";

std::vector<TString> fileNames = {"/scratch/EIC/Results/ML-Out/"+tag+"_real_E.root","/home/simon/Analysis/EIC_Analysis/Reg_"+tag2+".root"};
//std::vector<TString> fileNames = {"/scratch/EIC/Results/ML-Out/4layer_real_ETP.root","/home/simon/Analysis/EIC_Analysis/Reg_point-RealHits4layer.root"};
//std::vector<TString> fileNames = {"/scratch/EIC/Results/ML-Out/4layer_real_ETP.root","/home/simon/Analysis/EIC_Analysis/Reg_cell-RealHits4layer.root"};
//std::vector<TString> fileNames = {"/scratch/EIC/Results/ML-Out/test_cell_ETP.root"};


void ResolutionsJustE(){

//   TString outNamepng = "EResolutions-"+tag+".png";
  TString outNamepng = "EResolutions-"+tag2+".png";

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
  can->Divide(2,2);


//   ROOT::RDataFrame dfcheet("dataset/TestTree",fileNames[0]);

//   //Cheating
//   auto dfcheet_2 = dfcheet
//     .Define("ePred","DNN_CPU.eE")
//     .Define("ERes","100*(eE-ePred)/eE");

  ROOT::RDataFrame dfcheet("predict",fileNames[1]);

  auto dfcheet_2 = dfcheet
    .Define("ERes","100*(eE-ePred)/eE");



  int    eBins     = 200;
  int    resBins   = 100;

  double eMin = 6;
  double eMax = 18;

  double PResRange = 1.5;
  double TResRange = 0.002*1000;
  double EResRange = 0.05*100;

  auto Energy1D  = dfcheet_2.Histo1D({"ERes", ";Energy (gen-recon)/gen [%]", resBins, -EResRange, EResRange }, "ERes");

  auto Energy2D  = dfcheet_2.Histo2D({"Energy", "Energy Reconstruction;Energy gen [GeV] ;Energy recon [GeV]", eBins,     eMin,     eMax,     eBins,     eMin,     eMax     }, "eE",    "ePred");

  auto Energy2D2 = dfcheet_2.Histo2D({"ERes", ";Energy (gen-recon)/gen [%];Electron energy [GeV] ",       resBins, -EResRange, EResRange, eBins, eMin, eMax },"ERes",    "eE");


  (Energy2D2.GetPtr())->FitSlicesX();
  TH1D* energySigma = (TH1D*)(gDirectory->Get("ERes_2"));
   

  can->cd(1);
  energySigma->SetTitle("Energy Resolutions [%]");
  energySigma->SetStats(0);
  energySigma->SetMinimum(0);
  energySigma->SetMaximum(1);
  energySigma->Draw("colz");

  can->cd(2);
  gPad->SetLogz();
  Energy2D2.GetPtr()->Draw("colz");
  Energy2D2->SetStats(0);
  

  can->cd(3);
  gPad->SetLogz();
  Energy2D.GetPtr()->Draw("colz");
  Energy2D->SetStats(0);
  
  TH1D* ERes = (TH1D*)(Energy1D.GetPtr());
  can->cd(4);
  ERes->Draw();
  ERes->Fit("gaus");
  ERes->SetStats(0);

  can->SaveAs(outNamepng);

  TH1D* energySigma2 = (TH1D*)(energySigma->Clone("Clone2"));
  energySigma->Rebin(eBins);
  energySigma->Scale(1/(double)eBins);
  energySigma->Draw("hist");
  std::cout << energySigma->GetBinContent(1) << std::endl;


}
