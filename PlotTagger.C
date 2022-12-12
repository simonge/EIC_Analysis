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
  
    //ROOT::EnableImplicitMT();

  using namespace ROOT::Math;
  using namespace std;

  TFile* oFile = new TFile(outName,"RECREATE");


  // Input Data Chain
  TChain* t = new TChain("temp");
  t->Add(inName);

  ROOT::RDataFrame d0(*t);


  auto dCut = d0.Filter("Any(vector_filter)");

  double maxE  = 18;
  double maxQ2 = 0;
  double minQ2 = -9;
    
  auto gen_distributionEQ = d0.Histo2D({"gen_distributionE", "Generated Events;Electron Energy [GeV];log_{10}(Q^{2}) [GeV]", 100, 0.0, maxE,100,minQ2,maxQ2}, "eE","logQ2");

  auto cut_distributionEQ = dCut.Histo2D({"gen_distributionE", "Generated Events;Electron Energy [GeV];log_{10}(Q^{2}) [GeV]", 100, 0.0, maxE,100,minQ2,maxQ2}, "eE","logQ2");
    
  std::vector<ROOT::RDF::RResultPtr< ::TH2D>> rawHits;

  int   nModules  = 2;
  int   nLayers   = 4;
  float pixelSize = 0.055;
  float modSizeX[2]  = {150,120};
  float modSizeY[2]  = {100,75};
    
  TString cut;

  for(int i=0; i<nModules; i++){
    auto dModule = dCut.Define("modHits",[&i](ROOT::VecOps::RVec<int> module){
			    return module==i+1;
			  },{"moduleID"});
    

    for(int j=0; j<nLayers; j++){
      
      auto dLayer = dModule.Define("layerHits",[&j](ROOT::VecOps::RVec<int> layer){
			     return layer==j;
			   },{"layerID"})
	.Define("hits","cell_position[layerHits&&modHits]")
	.Define("hitX","xID[layerHits&&modHits]*0.055")
	.Define("hitY","yID[layerHits&&modHits]*0.055");


      TString histName;
      histName.Form("Module %d, Layer %d",i+1,j+1);
      
      rawHits.push_back(dLayer.Histo2D({histName, histName+"x hit position [mm];y hit position [mm]", 100, -modSizeX[i], modSizeX[i],100,-modSizeY[i],modSizeY[i]}, "hitX","hitY"));
			
    }
  }
  
  for(auto hist:rawHits){
    hist->Write();
  }

  gen_distributionEQ->Write();
  cut_distributionEQ->Write();
  oFile->Close();

}

//--------------------------------------------------------------------------------------------
