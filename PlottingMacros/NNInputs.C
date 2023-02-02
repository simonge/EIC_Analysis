#include "TString.h"
#include "TCanvas.h"
#include "ROOT/RDataFrame.hxx"
#include <vector>
#include "TH1.h"
//#include "ROOT.h"


using namespace std;

std::vector<double> pixSize = {0,55,110,220,440,880,1760};

std::vector<TString> fileNames = {"/scratch/EIC/Analysis/temp.root"};


void NNInputs(){

  TString outNamepng   = "plots/InputPos.png";


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

  // gStyle->SetPalette(kBird);

  TCanvas* can = new TCanvas("can","can",800,1400);
  can->Divide(1,2);

  ROOT::RDataFrame df("temp",fileNames[0]);

  auto df2 = df.Filter("(Tag1_4||Tag2_4)&&iFilter")
    .Define("Pos0","real_cut[0]")
    .Define("PosY","Pos0.y()")
    .Define("PosZ","Pos0.z()/1000")
    .Define("Vec0","real_vector[0]")
    .Define("VecX","Vec0.x()*100")
    .Define("VecY","Vec0.y()*100");
  

  int    bins  = 400;

  double zPosMin = -14;
  double zPosMax = -11;
  double yPosMin = -20;
  double yPosMax =  20;


  double xVecMin = -4.5;
  double xVecMax = -1.8;
  double yVecMin = -0.1;
  double yVecMax =  0.1;

  auto Pos = df2.Histo2D({"Pos_Inputs", "Input positions to neural net;y [mm];z [m]", bins, yPosMin, yPosMax,bins, zPosMin, zPosMax }, "PosY","PosZ");
  auto Vec = df2.Histo2D({"Vec_Inputs", "Input vectors to neural net;x [unit vec percent];y [unit vec percent]", bins, xVecMin, xVecMax,bins, yVecMin, yVecMax }, "VecX","VecY");
 
 
  can->cd(1);
  gPad->SetLogz();
  Pos->Draw("colz");
  can->cd(2);
  gPad->SetLogz();
  Vec->Draw("colz");

  
  can->SaveAs(outNamepng);


}
