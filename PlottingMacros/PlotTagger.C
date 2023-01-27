R__LOAD_LIBRARY(libfmt.so)
#include "fmt/core.h"

#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/VectorUtil.h"
//#include <Vector4D.h>

#include "ROOT/RDataFrame.hxx"
#include "TCanvas.h"
#include "TChain.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TVector3.h"

#include <tuple>
#include <vector>
#include <cmath>
 
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"

#pragma link C++ class std::vector<TVector3>;

  void PlotTagger(TString inName          = "/scratch/EIC/Analysis/temp.root",
		  TString outName         = "/scratch/EIC/Results/tempPlots.root"){
  
    ROOT::EnableImplicitMT();

  using namespace ROOT::Math;
  using namespace std;

  TFile* oFile = new TFile(outName,"RECREATE");


  // Input Data Chain
  TChain* t = new TChain("temp");
  t->Add(inName);

  ROOT::RDataFrame d0(*t);
  auto nGen = d0.Count();
  cout << nGen.GetValue() << endl;

  auto dCut = d0.Filter("iFilter&&Any(vector_filter)");
  //  auto dCut = d0.Filter("Any(vector_filter)");
  auto nCut = dCut.Count();
  cout << nCut.GetValue() << endl;

  double maxE  = 18;
  double maxQ2 = 0;
  double minQ2 = -9;

  std::vector<ROOT::RDF::RResultPtr< ::TH2D>> EQHists;
    
EQHists.push_back(d0.Histo2D({"gen_distributionEQ", "Generated Events;Electron Energy [GeV];log_{10}(Q^{2}) [GeV]", 100, 0.0, maxE,100,minQ2,maxQ2}, "eE","logQ2"));

 EQHists.push_back(dCut.Histo2D({"cut_distributionEQ", "Generated Events with hit and angle cut;Electron Energy [GeV];log_{10}(Q^{2}) [GeV]", 100, 0.0, maxE,100,minQ2,maxQ2}, "eE","logQ2"));

 EQHists.push_back(dCut.Filter("Tag1_1").Histo2D({"tag1_1_distributionEQ", "Generated Events with hit and angle cut;Electron Energy [GeV];log_{10}(Q^{2}) [GeV]", 100, 0.0, maxE,100,minQ2,maxQ2}, "eE","logQ2"));
 EQHists.push_back(dCut.Filter("Tag1_2").Histo2D({"tag1_2_distributionEQ", "Generated Events with hit and angle cut;Electron Energy [GeV];log_{10}(Q^{2}) [GeV]", 100, 0.0, maxE,100,minQ2,maxQ2}, "eE","logQ2"));
 EQHists.push_back(dCut.Filter("Tag1_3").Histo2D({"tag1_3_distributionEQ", "Generated Events with hit and angle cut;Electron Energy [GeV];log_{10}(Q^{2}) [GeV]", 100, 0.0, maxE,100,minQ2,maxQ2}, "eE","logQ2"));
 EQHists.push_back(dCut.Filter("Tag1_4").Histo2D({"tag1_4_distributionEQ", "Generated Events with hit and angle cut;Electron Energy [GeV];log_{10}(Q^{2}) [GeV]", 100, 0.0, maxE,100,minQ2,maxQ2}, "eE","logQ2"));
 
 EQHists.push_back(dCut.Filter("Tag2_1").Histo2D({"tag2_1_distributionEQ", "Generated Events with hit and angle cut;Electron Energy [GeV];log_{10}(Q^{2}) [GeV]", 100, 0.0, maxE,100,minQ2,maxQ2}, "eE","logQ2"));
 EQHists.push_back(dCut.Filter("Tag2_2").Histo2D({"tag2_2_distributionEQ", "Generated Events with hit and angle cut;Electron Energy [GeV];log_{10}(Q^{2}) [GeV]", 100, 0.0, maxE,100,minQ2,maxQ2}, "eE","logQ2"));
 EQHists.push_back(dCut.Filter("Tag2_3").Histo2D({"tag2_3_distributionEQ", "Generated Events with hit and angle cut;Electron Energy [GeV];log_{10}(Q^{2}) [GeV]", 100, 0.0, maxE,100,minQ2,maxQ2}, "eE","logQ2"));
 EQHists.push_back(dCut.Filter("Tag2_4").Histo2D({"tag2_4_distributionEQ", "Generated Events with hit and angle cut;Electron Energy [GeV];log_{10}(Q^{2}) [GeV]", 100, 0.0, maxE,100,minQ2,maxQ2}, "eE","logQ2"));
    
  std::vector<ROOT::RDF::RResultPtr< ::TH2D>> rawHits;

  float modSizeX[2]  = {150,120};
  float modSizeY[2]  = {100,75};
  int   positionbins = 200;

  //-------------------------------------------------------
  // Raw Hit distributions
  //-------------------------------------------------------

  // Module1    
  TString histName  = "Module 1, Layer 1";
  TString histTitle = "Module1Layer1";
  TString cutString = "layerID==0&&moduleID==1";

  rawHits.push_back(dCut.Define("Cut","layerID==0&&moduleID==1")
		    .Define("hitX","xID[Cut]*0.055")
		    .Define("hitY","yID[Cut]*0.055")
		    .Histo2D({histTitle, histName+";x hit position [mm];y hit position [mm]", positionbins, -modSizeX[0], modSizeX[0],positionbins,-modSizeY[0],modSizeY[0]}, "hitX","hitY"));

  histName  = "Module 1, Layer 2";
  histTitle = "Module1Layer2";
  rawHits.push_back(dCut.Define("Cut","layerID==1&&moduleID==1")
		    .Define("hitX","xID[Cut]*0.055")
		    .Define("hitY","yID[Cut]*0.055")
		    .Histo2D({histTitle, histName+";x hit position [mm];y hit position [mm]", positionbins, -modSizeX[0], modSizeX[0],positionbins,-modSizeY[0],modSizeY[0]}, "hitX","hitY"));

  histName  = "Module 1, Layer 3";
  histTitle = "Module1Layer3";
  rawHits.push_back(dCut.Define("Cut","layerID==2&&moduleID==1")
		    .Define("hitX","xID[Cut]*0.055")
		    .Define("hitY","yID[Cut]*0.055")
		    .Histo2D({histTitle, histName+";x hit position [mm];y hit position [mm]", positionbins, -modSizeX[0], modSizeX[0],positionbins,-modSizeY[0],modSizeY[0]}, "hitX","hitY"));

  histName  = "Module 1, Layer 4";
  histTitle = "Module1Layer4";
  rawHits.push_back(dCut.Define("Cut","layerID==3&&moduleID==1")
		    .Define("hitX","xID[Cut]*0.055")
		    .Define("hitY","yID[Cut]*0.055")
		    .Histo2D({histTitle, histName+";x hit position [mm];y hit position [mm]", positionbins, -modSizeX[0], modSizeX[0],positionbins,-modSizeY[0],modSizeY[0]}, "hitX","hitY"));

  //Module2
  histName  = "Module 2, Layer 1";
  histTitle = "Module2Layer1";
  rawHits.push_back(dCut.Define("Cut","layerID==0&&moduleID==2")
		    .Define("hitX","xID[Cut]*0.055")
		    .Define("hitY","yID[Cut]*0.055")
		    .Histo2D({histTitle, histName+";x hit position [mm];y hit position [mm]", positionbins, -modSizeX[1], modSizeX[1],positionbins,-modSizeY[1],modSizeY[1]}, "hitX","hitY"));

  histName  = "Module 2, Layer 2";
  histTitle = "Module2Layer2";
  rawHits.push_back(dCut.Define("Cut","layerID==1&&moduleID==2")
		    .Define("hitX","xID[Cut]*0.055")
		    .Define("hitY","yID[Cut]*0.055")
		    .Histo2D({histTitle, histName+";x hit position [mm];y hit position [mm]", positionbins, -modSizeX[1], modSizeX[1],positionbins,-modSizeY[1],modSizeY[1]}, "hitX","hitY"));

  histName  = "Module 2, Layer 3";
  histTitle = "Module2Layer3";
  rawHits.push_back(dCut.Define("Cut","layerID==2&&moduleID==2")
		    .Define("hitX","xID[Cut]*0.055")
		    .Define("hitY","yID[Cut]*0.055")
		    .Histo2D({histTitle, histName+";x hit position [mm];y hit position [mm]", positionbins, -modSizeX[1], modSizeX[1],positionbins,-modSizeY[1],modSizeY[1]}, "hitX","hitY"));

  histName  = "Module 2, Layer 4";
  histTitle = "Module2Layer4";
  rawHits.push_back(dCut.Define("Cut","layerID==3&&moduleID==2")
		    .Define("hitX","xID[Cut]*0.055")
		    .Define("hitY","yID[Cut]*0.055")
		    .Histo2D({histTitle, histName+";x hit position [mm];y hit position [mm]", positionbins, -modSizeX[1], modSizeX[1],positionbins,-modSizeY[1],modSizeY[1]}, "hitX","hitY"));

  


  for(auto hist:rawHits){
    hist->Write();
  }

  for(auto hist:EQHists){
    hist->Write();
  }
  oFile->Close();

}

//--------------------------------------------------------------------------------------------
