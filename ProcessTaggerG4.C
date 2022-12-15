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
#include "TLinearFitter.h"

#include <tuple>
#include <vector>
#include <cmath>

#include "DD4hep/Detector.h"
#include "DDRec/CellIDPositionConverter.h"
 
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"

#pragma link C++ class std::vector<TVector3>;
  
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

//-----------------------------------------------------------------------------------------
// Linear fit functor
//-----------------------------------------------------------------------------------------
struct fitPoints{
  fitPoints(int nLayers,int mod): maxLayer(nLayers), module(mod) { }
  
  std::pair<TVector3,double> operator()(const std::vector<TVector3>& positions, const ROOT::VecOps::RVec<float>& energies, const ROOT::VecOps::RVec<int>& moduleID, const ROOT::VecOps::RVec<int>& layerID) {
    
    TVector3 outVec;
    TVector3 outPos;
    double   outChi2 = 999999;
//     TLinearFitter* lf = new TLinearFitter(1);
//     lf->SetFormula( "pol1");
    TLinearFitter* lf = new TLinearFitter(2);
    lf->SetFormula( "1++x[0]++x[1]");
    
    ROOT::VecOps::RVec<double> x;
    ROOT::VecOps::RVec<double> y;
    ROOT::VecOps::RVec<double> z;
    ROOT::VecOps::RVec<double> e;

    for(uint i=0; i<positions.size(); i++){
      if(moduleID[i]!=module) continue;
      if(layerID[i]>maxLayer) continue;
      
//       x.push_back(positions[i].x());
//       y.push_back(positions[i].y());
      x.push_back(positions[i].x());
      y.push_back(positions[i].y());
      z.push_back(positions[i].z());
      e.push_back(energies[i]);

    }    
   
    if(e.size()>=2){
         
      double totalWeight = Sum(e);
      double xMean = Sum(x*e)/totalWeight;
      double yMean = Sum(y*e)/totalWeight;
      double zMean = Sum(z*e)/totalWeight;
      outPos.SetXYZ(xMean,yMean,zMean);

      x -= xMean;
      y -= yMean;
      z -= zMean;
      
      ROOT::VecOps::RVec<double> v;
      for(uint i=0; i<x.size(); i++){
	v.push_back(x[i]);
	v.push_back(y[i]);
      }



      lf->AssignData(e.size(), 2, &v[0], &z[0], &e[0]);    
      lf->Eval();
      TVectorD params;
//       TVectorD errors;
      lf->GetParameters(params);
//       params.Print();
//       lf->GetErrors(errors);
//       errors.Print();
      outChi2=lf->GetChisquare();
//       cout << outChi2 << endl;
      outVec.SetXYZ(params[0],params[1],params[2]);
//       cout << lf->GetNpoints() <<endl;
      //outVec.Print();
      //(outVec.Unit()).Print();
    }
    

    std::pair<TVector3,double> outPair = {outVec.Unit(),outChi2};

    return outPair;

  };
  
  void SetMaxLayers(int n){
    maxLayer = n;
  }
  
  
private: 
  int module;
  int maxLayer;
};
  

// Particle definitions and frame names.
struct partDetails{
  std::string fName;
  int pdg;
  int genStatus;
};

int beamID = 3;
int simID  = 1;
//int simID  = 1;

std::vector<partDetails> parts = {{"beamElectron",11,beamID},{"beamProton",2212,beamID},{"scatteredElectron",11,simID}};

//-----------------------------------------------------------------------------------------
// Main Analysis Call
//-----------------------------------------------------------------------------------------
//   void ProcessTaggerG4(TString inName          = "/scratch/EIC/G4out/derek/x_5_100_16.root",
// 		       TString outName         = "/scratch/EIC/Analysis/temp.root",
// 		       std::string compactName = "/home/simon/geant4/eic/ip6/eic_ip6.xml"){
//   void ProcessTaggerG4(TString inName          = "/scratch/EIC/G4out/qr_18x275_beam_out_*.edm4hep.root",
// 		       TString outName         = "/scratch/EIC/Analysis/temp.root",
// 		       std::string compactName = "/home/simon/EIC/epic/epic_ip6.xml"){
  void ProcessTaggerG4(TString inName          = "/scratch/EIC/G4out/lgen_18x275_beam_out_*.edm4hep.root",
		       TString outName         = "/scratch/EIC/Analysis/tempBrems.root",
		       std::string compactName = "/home/simon/EIC/epic/epic_ip6.xml"){
//   void ProcessTaggerG4(TString inName          = "/scratch/EIC/G4out/qr_18x275_beam_ReallyNoSolenoid_*.edm4hep.root",
// 		       TString outName         = "/scratch/EIC/Analysis/ReallyNoSolenoid.root",
// 		       std::string compactName = "/home/simon/EIC/epic/epic_ip6.xml"){
  
    ROOT::EnableImplicitMT();

  using namespace ROOT::Math;
  using namespace std;

  // Input Data Chain
  TChain* t = new TChain("events");
  t->Add(inName);

  ROOT::RDataFrame d0(*t, {"TaggerTrackerHits", "MCParticles"});
  //d0.Range(0,1000); // Limit events to analyse 

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
      if(h.generatorStatus!=simID) continue;
      return h.vertex;
      break;
    }
    return edm4hep::Vector3d();
  };

  // -------------------------
  // Beam Time 
  // -------------------------
  auto beamTime = [&](const std::vector<edm4hep::MCParticleData>& evt) {    

    for (const auto& h : evt) {
      if(h.generatorStatus!=simID) continue;
      return h.time;
    }
    return (float)0.0;
  };


  // -------------------------
  // Real Hit Vector 
  // -------------------------
  auto real_vector = [&](const std::vector<edm4hep::SimTrackerHitData>& hits) {
    
    std::vector<TVector3> vectors;

    for (const auto& h : hits) {
      TVector3 result(h.momentum.x,h.momentum.y,h.momentum.z);
      result *= 1/result.Mag();
      vectors.push_back(result);
    }
    return vectors;
  };

  // -------------------------
  // Real Hit Positions
  // -------------------------
  auto real_position = [&](const std::vector<edm4hep::SimTrackerHitData>& hits) {
    
    std::vector<TVector3> positions;

    for (const auto& h : hits) {
      TVector3 result(h.position.x,h.position.y,h.position.z);
      positions.push_back(result);
    }
    return positions;
  };

  // -------------------------
  // Cell Hit Positions
  // -------------------------
  auto cell_position = [&](const std::vector<edm4hep::SimTrackerHitData>& hits) {
    
    std::vector<TVector3> positions;

    for (const auto& h : hits) {
      auto pos1 = cellid_converter.position(h.cellID);
      TVector3 result(pos1.x()*10,pos1.y()*10,pos1.z()*10);
      positions.push_back(result);
    }
    return positions;
  };

  // -------------------------
  // Vector XZ origin cut
  // -------------------------
//   auto vector_cut = [&](const std::vector<TVector3>& positions, const std::vector<TVector3>& vectors) {
    
//     if(positions.size()<=0) return TVector3();

//     auto vec = vectors[0];
//     auto pos = positions[0];

//     float con = pos.x()/vec.x();

//     auto   xPoint = pos-con*vec;
//     return xPoint;

//   };

  auto vector_cut = [&](const std::vector<TVector3>& positions, const std::vector<TVector3>& vectors) {
    
    std::vector<TVector3> points;

    for(uint i=0; i<positions.size(); i++){
      
      float con = positions[i].x()/vectors[i].x();

      points.push_back(positions[i]-con*vectors[i]);

    }
    return points;

  };

//   auto vector_cut2 = [&](const std::vector<TVector3>& positions, const TVector3& vectors) {
    
//     std::vector<TVector3> points;
      
//     float con = positions[0].x()/vectors.x();

//     points.push_back(positions[0]-con*vectors);

//     return points;

//   };
  
  auto vector_filter = [&](const std::vector<TVector3>& positions) {
    
    std::vector<bool> good;

    for(auto pos: positions){
      if(pos.z()<-8000&&pos.z()>-15000&&abs(pos.y())<100) good.push_back(true);
      else good.push_back(false);

    }
    return good;

  };
  

  auto ids = detector.readout("TaggerTrackerHits").idSpec().fields();
  std::vector<std::string> ID_Vec;
  std::vector<std::string> Part_Vec;

//       //Do some calculations
//  auto d1 = d0.Define("Hit Filter",[](auto conts){
   auto d1 = d0.Define("nHits", "TaggerTrackerHits.size()")
     .Define("real_position",     real_position,           {"TaggerTrackerHits"})
     .Define("cell_position",     cell_position,           {"TaggerTrackerHits"})
     .Define("cell_vector2",      "cell_position[0]-cell_position")
     .Define("cell_cut",          vector_cut,              {"cell_position","cell_vector2"})
     .Define("real_time",         "TaggerTrackerHits.time")
     .Define("real_EDep",         "TaggerTrackerHits.EDep")
     .Define("real_vector",       real_vector,             {"TaggerTrackerHits"})
     .Define("vector_cut",        vector_cut,              {"real_position","real_vector"})
     .Define("vector_filter",     vector_filter,           {"vector_cut"});

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
   

//    d1 = d1.Define("CapHit",getSubID("system",detector,"Cap_Track"))
//      .Define("CapPass","CapHit.size()");
   d1 = d1.Define("iFilter","(TMath::Pi()-scatteredElectron.Theta())<0.0105")
     .Define("vertex", beamVertex , {"MCParticles"})
     .Define("nParticles",        "MCParticles.size()")
     .Define("time", beamTime , {"MCParticles"})
     .Define("eE", "scatteredElectron.energy()")
     .Define("thetaE", "acos((scatteredElectron.Vect().Unit()).Dot(beamElectron.Vect().Unit()))")
     //     .Define("phiE", "acos((((beamElectron.Vect()).Cross(TVector3(0,1,0))).Unit()).Dot((scatteredElectron.Vect()).Cross(TVector3(0,1,0)).Unit()))")
     .Define("pseudorapidity", "scatteredElectron.eta()")
     .Define("scatteredV", "beamElectron-scatteredElectron")
     .Define("thetaV", "scatteredV.theta()")
     .Define("phiV", "scatteredV.phi()")
     .Define("qE", "scatteredV.energy()")
     .Define("Q2", "-scatteredV.M2()")
     .Define("logQ2", "log10(Q2)")
     .Define("Tag1_1","Any(layerID==0&&moduleID==1)")
     .Define("Tag1_2","Tag1_1&&Any(layerID==1&&moduleID==1)")
     .Define("Tag1_3","Tag1_2&&Any(layerID==2&&moduleID==1)")
     .Define("Tag1_4","Tag1_3&&Any(layerID==3&&moduleID==1)")
     .Define("Tag2_1","Any(layerID==0&&moduleID==2)")
     .Define("Tag2_2","Tag2_1&&Any(layerID==1&&moduleID==2)")
     .Define("Tag2_3","Tag2_2&&Any(layerID==2&&moduleID==2)")
     .Define("Tag2_4","Tag2_3&&Any(layerID==3&&moduleID==2)")
     .Define("x11","xID[layerID==0&&moduleID==1]")
     .Define("y11","yID[layerID==0&&moduleID==1]")
     .Define("Npix11","x11.size()")
     .Define("x12","xID[layerID==1&&moduleID==1]")
     .Define("y12","yID[layerID==1&&moduleID==1]")
     .Define("Npix12","x12.size()")
     .Define("x13","xID[layerID==2&&moduleID==1]")
     .Define("y13","yID[layerID==2&&moduleID==1]")
     .Define("Npix13","x13.size()")
     .Define("x14","xID[layerID==3&&moduleID==1]")
     .Define("y14","yID[layerID==3&&moduleID==1]")
     .Define("Npix14","x14.size()")
     .Define("x21","xID[layerID==0&&moduleID==2]")
     .Define("y21","yID[layerID==0&&moduleID==2]")
     .Define("Npix21","x21.size()")
     .Define("x22","xID[layerID==1&&moduleID==2]")
     .Define("y22","yID[layerID==1&&moduleID==2]")
     .Define("Npix22","x22.size()")
     .Define("x23","xID[layerID==2&&moduleID==2]")
     .Define("y23","yID[layerID==2&&moduleID==2]")
     .Define("Npix23","x23.size()")
     .Define("x24","xID[layerID==3&&moduleID==2]")
     .Define("y24","yID[layerID==3&&moduleID==2]")
     .Define("Npix24","x24.size()");

//    d1 = d1.Define("fit_temp", fitPoints(3,1) , {"vector_cut","real_EDep","moduleID","layerID"})
//      .Define("fit_vector","fit_temp.first")
//      .Define("fit_chi2","fit_temp.second");
   
   ROOT::RDF::RSnapshotOptions opts;
   opts.fMode = "UPDATE";

   d1.Snapshot("ids",outName,ID_Vec);
   d1.Snapshot("parts",outName,Part_Vec,opts);

   

   std::vector<std::string> Out_Vec = {"vertex","nParticles","nHits","real_position","cell_position","real_vector","iFilter","vector_cut","vector_filter","thetaV","thetaE","phiV","qE","eE","logQ2","Tag1_1","Tag1_2","Tag1_3","Tag1_4","Tag2_1","Tag2_2","Tag2_3","Tag2_4","x11","y11","x12","y12","x13","y13","x14","y14","x21","y21","x22","y22","x23","y23","x24","y24"};//,"fit_vector","fit_chi2"};

   for(auto a:ID_Vec) Out_Vec.push_back(a);
   for(auto a:Part_Vec) Out_Vec.push_back(a);

   d1.Snapshot("temp",outName,Out_Vec,opts);


   //Hits distribution/layer
   //Acceptance 
  

//    auto d2 = d1.Filter("(Npix11!=0 && Npix12!=0)");// && Npix13!=0)");
//    d1.Snapshot("input",outName,{"nhitsT1","nhits2","vertex","Q2","logQ2","ex","ey","ez","eTheta","ePhi","eE","qx","qy","qz","qTheta","qPhi","qE","pseudorapidity"});
//    //d2.Snapshot("detector1",outname,{"x11","y11","Npix11","x12","y12","Npix12","x13","y13","Npix13","Q2","logQ2","eE","qTheta","qPhi","eTheta","ePhi","pseudorapidity"},opts);
//    //   d1.Filter("B2BPass!=0")
//    //  .Snapshot("detectors",outname,{"nhitsT1","layerID","vertex","Q2","logQ2","eE","qTheta","qPhi","eTheta","ePhi","pseudorapidity","real_position_x","real_position_y","real_position_z","real_vector"},opts);
//    d2.Snapshot("detector1",outName,{"vertex","x11","y11","Npix11","x12","y12","Npix12","x13","y13","Npix13","Q2","logQ2","eE","qTheta","qPhi","eTheta","ePhi","pseudorapidity","real_position_x","real_position_y","real_position_z","real_vector"},opts);


}

//--------------------------------------------------------------------------------------------
