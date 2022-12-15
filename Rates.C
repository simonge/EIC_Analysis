#include "TString.h"
#include "TCanvas.h"
#include "ROOT/RDataFrame.hxx"
#include <vector>
#include "TH1.h"
//#include "ROOT.h"


using namespace std;

std::vector<double> pixSize = {0,55,110,220,440,880,1760};

std::vector<TString> fileNames = {"/scratch/EIC/Analysis/temp.root","/scratch/EIC/Analysis/tempBrems.root"};


void Rates(){

  TString outName      = "Rates.pdf";
  TString outNamepng   = "Rates.png";

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
  gStyle->SetLabelSize(0.042,"X");
  gStyle->SetLabelSize(0.042,"Y");
  gStyle->SetLabelSize(0.042,"Z");
  gStyle->SetTitleSize(0.06,"X");
  gStyle->SetTitleSize(0.06,"Y");
  gStyle->SetTitleSize(0.06,"Z");

  gStyle->SetTitleOffset(1.0,"Y");

  gStyle->SetPalette(kBird);

  double Nevents   = 10000000;
  double meanQR    = 0.00374;
  double meanBrems = 11.65;
  double bunchFreq = 10e-9;

  double pixelArea        = 0.055*0.055;

  TCanvas* can = new TCanvas("can","can",1800,1400);
  can->Divide(2,2);


  ROOT::RDataFrame df("temp",fileNames[0]);

  auto df2 = df.Filter("iFilter")
    .Define("Tag1X","xID[moduleID==1&&layerID==0]*0.055")
    .Define("Tag1Y","yID[moduleID==1&&layerID==0]*0.055")
    .Define("Tag2X","xID[moduleID==2&&layerID==0]*0.055")
    .Define("Tag2Y","yID[moduleID==2&&layerID==0]*0.055");
  
  ROOT::RDataFrame dfb("temp",fileNames[1]);

  auto dfb2 = dfb.Filter("iFilter")
    .Define("Tag1X","xID[moduleID==1&&layerID==0]*0.055")
    .Define("Tag1Y","yID[moduleID==1&&layerID==0]*0.055")
    .Define("Tag2X","xID[moduleID==2&&layerID==0]*0.055")
    .Define("Tag2Y","yID[moduleID==2&&layerID==0]*0.055");

  int    tagBins  = 150;

  double t1x = 150;
  double t2x = 120;
  double t1y = 100;
  double t2y = 75;

  auto Tag1_QR = df2.Histo2D({"Tag1_QR", "Tagger 1 QR Hit Distribution [Hz/ 55#mum pixel];x [mm];y [mm]", tagBins, -t1x, t1x,tagBins, -t1y, t1y }, "Tag1X","Tag1Y");
  auto Tag2_QR = df2.Histo2D({"Tag2_QR", "Tagger 2 QR Hit Distribution [Hz/ 55#mum pixel];x [mm];y [mm]", tagBins, -t2x, t2x,tagBins, -t2y, t2y }, "Tag2X","Tag2Y");
  auto Tag1_Brem = dfb2.Histo2D({"Tag1_Brem", "Tagger 1 Brem Hit Distribution [Hz/ 55#mum pixel];x [mm];y [mm]", tagBins, -t1x, t1x,tagBins, -t1y, t1y }, "Tag1X","Tag1Y");
  auto Tag2_Brem = dfb2.Histo2D({"Tag2_Brem", "Tagger 2 Brem Hit Distribution [Hz/ 55#mum pixel];x [mm];y [mm]", tagBins, -t2x, t2x,tagBins, -t2y, t2y }, "Tag2X","Tag2Y");
 
  double binArea_1        = Tag1_QR->GetXaxis()->GetBinWidth(1)*Tag1_QR->GetYaxis()->GetBinWidth(1);
  double binArea_2        = Tag2_QR->GetXaxis()->GetBinWidth(1)*Tag2_QR->GetYaxis()->GetBinWidth(1);

  double areaScalePixel_1 = pixelArea/binArea_1; //pixel
  double areaScalemm2_1   = 1/binArea_1; //mm2
  double areaScalePixel_2 = pixelArea/binArea_2; //pixel
  double areaScalemm2_2   = 1/binArea_2; //mm2

  Tag1_QR->Scale(meanQR/(Nevents*bunchFreq));
  Tag2_QR->Scale(meanQR/(Nevents*bunchFreq));
  Tag1_Brem->Scale(meanBrems/(Nevents*bunchFreq));
  Tag2_Brem->Scale(meanBrems/(Nevents*bunchFreq));

  Tag1_Brem->Scale(areaScalePixel_1); //rate per bin
  Tag1_QR  ->Scale(areaScalePixel_1);
  Tag2_Brem->Scale(areaScalePixel_2); //rate per bin
  Tag2_QR  ->Scale(areaScalePixel_2);

  can->cd(1);
  gPad->SetLogz();
  Tag1_QR->Draw("colz");
  can->cd(2);
  gPad->SetLogz();
  Tag2_QR->Draw("colz");
  can->cd(3);
  gPad->SetLogz();
  Tag1_Brem->Draw("colz");
  can->cd(4);
  gPad->SetLogz();
  Tag2_Brem->Draw("colz");

  
  can->SaveAs(outName);
  can->SaveAs(outNamepng);

  

}
