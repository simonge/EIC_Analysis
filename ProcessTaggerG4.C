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

#include <algorithm>
#include <iterator>
#include <tuple>
#include <vector>
#include <cmath>

#include "DD4hep/Detector.h"
#include "DDRec/CellIDPositionConverter.h"

#include "dd4pod/Geant4ParticleCollection.h"
#include "dd4pod/TrackerHitCollection.h"
#include "dd4pod/CalorimeterHitCollection.h"
 
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
  
//-----------------------------------------------------------------------------------------
// Grab Component functor
//-----------------------------------------------------------------------------------------
  struct getSubID{
    getSubID(std::string cname, dd4hep::Detector& det, std::string rname = "TaggerTrackerHits") : componentName(cname), detector(det), readoutName(rname){}
    
    ROOT::VecOps::RVec<int> operator()(const std::vector<edm4hep::SimTrackerHitData>& evt) {
      auto decoder = detector.readout(readoutName).idSpec().decoder();
      auto indexID = decoder->index(componentName);
      ROOT::VecOps::RVec<int> result;
      for(const auto& h: evt) {
	result.push_back(decoder->get(h.cellID,indexID));      
      }
      return result;    
    };
    
    void SetComponent(std::string cname){
      componentName = cname;
    }
    void SetReadout(std::string rname){
      readoutName = rname;
    }

    private: 
    std::string componentName;
    dd4hep::Detector& detector;
    std::string readoutName;
  };

//-----------------------------------------------------------------------------------------
// Grab Particle functor
//-----------------------------------------------------------------------------------------
struct getParticle{
  getParticle(int genStat, int pdg) : generatorStatus(genStat), PDG(pdg){}
  
  ROOT::Math::PxPyPzMVector operator()(const std::vector<edm4hep::MCParticleData>& evt) {
    for(const auto& h: evt) {
      if(h.PDG!=PDG || h.generatorStatus!=generatorStatus) continue;
      return ROOT::Math::PxPyPzMVector(h.momentum.x,h.momentum.y,h.momentum.z,h.mass);      
    }
    return ROOT::Math::PxPyPzMVector();
  };
  
  void SetGeneratorStatus(int genStat){
    generatorStatus = genStat;
  }
  void SetPDG(int pdg){
    PDG = pdg;
  }
  
private: 
  int generatorStatus;
  int PDG;
};

// Particle definitions and frame names.
struct partDetails{
  std::string fName;
  int pdg;
  int genStatus;
};

int beamID = 3;
int simID  = 1;

std::vector<partDetails> parts = {{"beamElectron",11,beamID},{"beamProton",2212,beamID},{"scatteredElectron",11,simID}};

//-----------------------------------------------------------------------------------------
// Main Analysis Call
//-----------------------------------------------------------------------------------------
  void ProcessTaggerG4(TString inName          = "/scratch/EIC/G4out/derek/x_5_100_16.root",
		       TString outName         = "/scratch/EIC/Analysis/temp.root",
		       std::string compactName = "/home/simon/geant4/eic/ip6/eic_ip6.xml"){
  
  ROOT::EnableImplicitMT();

  using namespace ROOT::Math;
  using namespace std;

  // Input Data Chain
  TChain* t = new TChain("events");
  t->Add(inName);

  ROOT::RDataFrame d0(*t, {"TaggerTrackerHits", "MCParticles"});
  //  d0.Range(0,1000); // Limit events to analyse 

  // -------------------------
  // Get the DD4hep instance
  // Load the compact XML file
  // Initialize the position converter tool
  dd4hep::Detector& detector = dd4hep::Detector::getInstance();
  detector.fromCompact(compactName);
  dd4hep::rec::CellIDPositionConverter cellid_converter(detector);
  // -------------------------

  //--------------------------------------------------------------------------------------------
  //Lambda Functions
  //--------------------------------------------------------------------------------------------
  
  // -------------------------
  // Beam Vector 
  // -------------------------
  auto beamVertex = [&](const std::vector<edm4hep::MCParticleData>& evt) {    

    for (const auto& h : evt) {
      // The actual hit position:
      if(h.generatorStatus!=simID) continue;
      return Cartesian3D(h.vertex.x,h.vertex.y,h.vertex.z);
      break;
    }
    return Cartesian3D();
  };


  // -------------------------
  // Real Hit Vector 
  // -------------------------
  auto real_vector = [&](const std::vector<edm4hep::SimTrackerHitData>& hits) {
    
    for (const auto& h : hits) {
      TVector3 result(h.momentum.x,h.momentum.y,h.momentum.z);
      result *= 1/result.Mag();
      return result;
    }
    return TVector3();
  };

  auto ids = detector.readout("TaggerTrackerHits").idSpec().fields();
  std::vector<std::string> ID_Vec;
  std::vector<std::string> Part_Vec;


//       //Do some calculations
//  auto d1 = d0.Define("Hit Filter",[](auto conts){
   auto d1 = d0.Define("nHits", "TaggerTrackerHits.size()")
     .Define("nParticles", "MCParticles.size()")
     .Define("real_position_x",   "TaggerTrackerHits.position.x")
     .Define("real_position_y",   "TaggerTrackerHits.position.y")
     .Define("real_position_z",   "TaggerTrackerHits.position.z")
     .Define("real_vector",     real_vector,           {"TaggerTrackerHits"});

   for(auto part: parts){
     std::string colName = part.fName;
     Part_Vec.push_back(colName);
     d1 = d1.Define(colName,getParticle(part.genStatus,part.pdg),{"MCParticles"});
   }

   for(auto id: ids){
     std::string colName = id.first+"ID";
     ID_Vec.push_back(colName);
     d1 = d1.Define(colName,getSubID(id.first,detector),{"TaggerTrackerHits"});
   }
   

   d1 = d1.Define("vertex", beamVertex , {"MCParticles"})
     .Define("eE", "scatteredElectron.energy()")
     .Define("pseudorapidity", "scatteredElectron.eta()")
     .Define("scatteredV", "beamElectron-scatteredElectron")
     .Define("qE", "scatteredV.energy()")
     .Define("Q2", "-scatteredV.M2()")
     .Define("logQ2", "log10(Q2)")

     .Define("x11","xID[layerID==0&&systemID==195]")
     .Define("y11","yID[layerID==0&&systemID==195]")
     .Define("Npix11","x11.size()")
     .Define("x12","xID[layerID==1&&systemID==195]")
     .Define("y12","yID[layerID==1&&systemID==195]")
     .Define("Npix12","x12.size()")
     .Define("x13","xID[layerID==2&&systemID==195]")
     .Define("y13","yID[layerID==2&&systemID==195]")
    .Define("Npix13","x13.size()");


   d1 = d1.Define("CapHit",getSubID("system",detector,"Cap_Track"));


   ROOT::RDF::RSnapshotOptions opts;
   opts.fMode = "UPDATE";

   d1.Snapshot("ids",outName,ID_Vec);
   d1.Snapshot("parts",outName,Part_Vec,opts);
   d1.Snapshot("temp",outName,{"vertex","nParticles","nHits"},opts);


  

//    auto d2 = d1.Filter("(Npix11!=0 && Npix12!=0)");// && Npix13!=0)");
//    d1.Snapshot("input",outName,{"nhitsT1","nhits2","vertex","Q2","logQ2","ex","ey","ez","eTheta","ePhi","eE","qx","qy","qz","qTheta","qPhi","qE","pseudorapidity"});
//    //d2.Snapshot("detector1",outname,{"x11","y11","Npix11","x12","y12","Npix12","x13","y13","Npix13","Q2","logQ2","eE","qTheta","qPhi","eTheta","ePhi","pseudorapidity"},opts);
//    //   d1.Filter("B2BPass!=0")
//    //  .Snapshot("detectors",outname,{"nhitsT1","layerID","vertex","Q2","logQ2","eE","qTheta","qPhi","eTheta","ePhi","pseudorapidity","real_position_x","real_position_y","real_position_z","real_vector"},opts);
//    d2.Snapshot("detector1",outName,{"vertex","x11","y11","Npix11","x12","y12","Npix12","x13","y13","Npix13","Q2","logQ2","eE","qTheta","qPhi","eTheta","ePhi","pseudorapidity","real_position_x","real_position_y","real_position_z","real_vector"},opts);


}

//--------------------------------------------------------------------------------------------
