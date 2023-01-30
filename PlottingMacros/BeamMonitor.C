#include "TString.h"
#include "TCanvas.h"
#include "ROOT/RDataFrame.hxx"
#include <vector>
#include "TH1.h"
//#include "ROOT.h"


using namespace std;

std::vector<double> pixSize = {0,55,110,220,440,880,1760};

std::vector<TString> fileNames = {"/scratch/EIC/Jarda/maps_basic_v1.root"};
//std::vector<TString> fileNames = {"/scratch/EIC/Results/ML-Out/test_cell_ETP.root"};


void BeamMonitor(){

  TString outNamepng = "plots/BeamPos.png";

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

  ROOT::RDataFrame dfcheet("s2_tracks",fileNames[0]);

  //Cheating
  auto dfcheet_2 = dfcheet.Define("MeanY","ROOT::VecOps::Mean(pos_y)");


  int     bins     = 200;
  double  range    = 20;
  
  auto YPos    = dfcheet_2.Histo1D({"MeanT", ";#meanY [rad]", bins, -range, range }, "MeanY");

  YPos.GetPtr()->Draw("colz");
  
  can->SaveAs(outNamepng);


}
