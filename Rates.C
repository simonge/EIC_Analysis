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
  TString outQ2Name    = "Q2Rates.pdf";
  TString outQ2Namepng = "Q2Rates.png";

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

  double Nevents    = 10000000;
  double meanQR     = 0.00374;
  double meanQRtrig = 0.5;
  double meanBrems  = 11.65;
  double bunchFreq  = 10e-9;

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

  int    Q2bins  = 400;
  double Q2min   = -15;
  double Q2max   = 2;

  TCanvas* can2 = new TCanvas("can","can",1800,900);
  can2->Divide(2,1);

  // Q2 rate plots
  auto QR_Q2    = df2.Filter("Tag1_4||Tag2_4").Histo1D({"QR_Q2",";logQ2;Event Rate per trigger",Q2bins,Q2min,Q2max},"logQ2");
  auto Brem_Q2  = dfb2.Filter("Tag1_4||Tag2_4").Histo1D({"Brem_Q2",";logQ2;Event Rate per trigger",Q2bins,Q2min,Q2max},"logQ2");

  QR_Q2  ->Scale(meanQRtrig/Nevents);
  Brem_Q2->Scale(meanBrems/Nevents);

  auto sum_Q2 = (TH1*)QR_Q2->Clone("Total_Q2");
  sum_Q2->Add(QR_Q2.GetPtr(),Brem_Q2.GetPtr());

  auto ratio_Q2 = (TH1*)QR_Q2->Clone("Total_Q2");
  ratio_Q2->Divide(QR_Q2.GetPtr(),sum_Q2);

  QR_Q2   ->SetLineColor(kRed);
  Brem_Q2 ->SetLineColor(kBlue);
  sum_Q2  ->SetLineColor(kGreen);

  ratio_Q2->GetYaxis()->SetTitle("Signal Fration");
  
  can2->cd(1);
  gPad->SetLogy();
  sum_Q2 ->Draw("hist");
  QR_Q2  ->Draw("hist same");
  Brem_Q2->Draw("hist same");

  can2->cd(2);
  ratio_Q2->Draw();

  can2->SaveAs(outQ2Name);
  can2->SaveAs(outQ2Namepng);

}
