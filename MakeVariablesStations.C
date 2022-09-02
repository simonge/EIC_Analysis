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


//--------------------------------------------------------------------------------------------

  void MakeVariablesStations(const char* fname = "/scratch/EIC/G4out/CentralTagFrontWindow/qr_18x275_beam_out_0.root"){
  //  void MakeVariablesFPUpdate(const char* fname = "/scratch/EIC/G4out/FP_Tagger_Test.root"){
  //  void MakeVariablesFP(const char* fname = "/scratch/EIC/G4out/lgen_18x275_beam_out_*.root"){
//   void MakeVariablesNew(const char* fname = "/scratch/EIC/G4out/qr_18x275_out_more.root"){
  
  ROOT::EnableImplicitMT();
  TString outname = "/scratch/EIC/Analysis/CentralTagFrontWindow_events.root";
  //  TString outname = "/scratch/EIC/Analysis/qr_out18x275_events.root";

  using namespace ROOT::Math;
  using namespace std;

  //Chain
  TChain* t = new TChain("events");
  t->Add(fname);

  ROOT::RDataFrame d0(*t, {"TaggerTrackerHits", "MCParticles"});
  //d0.Range(0,10000);

  // -------------------------
  // Get the DD4hep instance
  // Load the compact XML file
  // Initialize the position converter tool
  dd4hep::Detector& detector = dd4hep::Detector::getInstance();
  detector.fromCompact("/home/simon/geant4/eic/ip6/eic_ip6.xml");
  dd4hep::rec::CellIDPositionConverter cellid_converter(detector);
  // -------------------------


  // -------------------------
  //How to identify the detector and layer from the cellID
  auto decoder = detector.readout("TaggerTrackerHits").idSpec().decoder();
  fmt::print("{}\n", decoder->fieldDescription());
  auto x_index         = decoder->index("x");
  auto y_index         = decoder->index("y");
  auto layer_index     = decoder->index("layer");
  auto detector_index  = decoder->index("system");
  // -------------------------


  //--------------------------------------------------------------------------------------------
  //Lambda Functions
  //--------------------------------------------------------------------------------------------
  
  // Simple lambda to define nhits branch
  // -------------------------
  auto nhits  = [](const std::vector<edm4hep::SimTrackerHitData>& evt)     { return (int)evt.size(); };
  auto nhits2 = [](const std::vector<edm4hep::MCParticleData>& evt) { return (int)evt.size(); }; 
  
  // -------------------------
  // ID of detector hit
  // -------------------------
  auto detectorID = [&](const std::vector<edm4hep::SimTrackerHitData>& evt) {         
    ROOT::VecOps::RVec<int> result;    
    for(const auto& h: evt) {       
      result.push_back(decoder->get(h.cellID,detector_index));            
    }     
    return result;     
  };

  // -------------------------
  // ID of calorimeter detector hit
  // -------------------------
  auto detectorIDC = [&](const std::vector<dd4pod::CalorimeterHitData>& evt) {         
    ROOT::VecOps::RVec<int> result;    
    for(const auto& h: evt) result.push_back(decoder->get(h.cellID,detector_index));            
    return result;     
  };

  // -------------------------
  // ID of layer hit
  // -------------------------
  auto layer = [&](const std::vector<edm4hep::SimTrackerHitData>& evt) {        
    ROOT::VecOps::RVec<int> result;    
    for(const auto& h: evt) {      
      result.push_back(decoder->get(h.cellID,layer_index));
    }     
    return result;    
  };

  // -------------------------
  // X-Pixel Hit
  // -------------------------
  auto xpix = [&](const std::vector<edm4hep::SimTrackerHitData>& evt) {         
    ROOT::VecOps::RVec<int> result;     
    for(const auto& h: evt) {       
      result.push_back(decoder->get(h.cellID,x_index));                    
    }     
    return result;     
  };

  // -------------------------  
  // Y-Pixel Hit
  // -------------------------
  auto ypix = [&](const std::vector<edm4hep::SimTrackerHitData>& evt) {    
    ROOT::VecOps::RVec<int> result;
    for(const auto& h: evt) {
      result.push_back(decoder->get(h.cellID,y_index));      
    }
    return result;    
  };

  // -------------------------
  // Electron Vector 
  // -------------------------
  auto electronVector = [&](const std::vector<edm4hep::MCParticleData>& evt) {    
    double px=0;
    double py=0;
    double pz=0;
    double mass=0;

    for (const auto& h : evt) {
      // The actual hit position:
      //if(h.ID!=7) continue; //Derek
      if(h.PDG!=11) continue;
      if(h.generatorStatus!=1) continue;
      
      px   = h.momentum.x;
      py   = h.momentum.y;
      pz   = h.momentum.z;
      mass = h.mass;

      //cout << px << " " <<py << " " <<pz << " " << mass << endl;

      break;
    }

    PxPyPzMVector result(px,py,pz,mass);
    return result;
  };

  // -------------------------
  // Beam Vector 
  // -------------------------
  auto beamVertex = [&](const std::vector<edm4hep::MCParticleData>& evt) {    
    double x=0;
    double y=0;
    double z=0;

    for (const auto& h : evt) {
      // The actual hit position:
      //if(h.ID!=7) continue; //Derek
      
      x   = h.vertex.x;
      y   = h.vertex.y;
      z   = h.vertex.z;
      
      
      break;
    }

    Cartesian3D result(x,y,z);
    return result;
  };

  // -------------------------
  // Beam Vector 
  // -------------------------
  auto beamVector = [&](const std::vector<edm4hep::MCParticleData>& evt) {    
    double px=0;
    double py=0;
    double pz=0;
    double mass=0;

    Double_t eMass = 0.000511;
    for (const auto& h : evt) {
      // The actual hit position:
      //if(h.ID!=7) continue; //Derek
      if(h.PDG!=11) continue;
      if(h.generatorStatus!=3) continue;
   
      px   = h.momentum.x;
      py   = h.momentum.y;
      pz   = h.momentum.z;
      mass = h.mass;

      //cout << px << " " <<py << " " <<pz << endl;

      break;
    }

    PxPyPzMVector result(px,py,pz,mass);
    return result;
  };

  // -------------------------
  // Calculate Q^2 
  // -------------------------
  auto scatterVector = [&](const PxPyPzMVector vec) {    
    Double_t eMass = 0.000511;
    //Double_t beamE = 10;
    Double_t beamE = 18;
    PxPyPzEVector beamVec(0,0,-sqrt(eMass*eMass+beamE*beamE),beamE);
    auto tVec = beamVec - vec;  
    return tVec;
  };

  // -------------------------
  // Project Real Hit Position 
  // -------------------------
  // auto project_position = [&](const std::vector<edm4hep::SimTrackerHitData>& hits) {
//     dd4pod::VectorXYZ result;
//     for (const auto& h : hits) {
//       // The actual hit position:
//       auto pos0 = (h.position);
//       auto vec0 = h.momentum.scale(1/h.momentum.mag());
//       auto pos1 = pos0.subtract(vec0.scale((pos0.z+14865)/vec0.z));// cellid_converter.position(h.cellID);
//       result = pos1;
//       break;
//     }
//     return result;
//   };

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


  // -------------------------
  // Recon Hit Position 
  // -------------------------
  auto recon_position = [&](const std::vector<edm4hep::SimTrackerHitData>& hits) {
    std::vector<dd4pod::VectorXYZ> result;
    for (const auto& h : hits) {
      auto pos1 = cellid_converter.position(h.cellID);
      //cout << pos1 << endl;
      result.push_back(dd4pod::VectorXYZ(pos1.x(),pos1.y(),pos1.z()));
    }
    return result;
  };


  // -------------------------
  // Recon Hit Vector 
  // -------------------------
  auto recon_vector = [&](const std::vector<dd4pod::VectorXYZ>& hits) {
    dd4pod::VectorXYZ result;
    int i=0;
    for (const auto& h : hits) {
      i++;
      if(i==1) result = h;
      if(i>=1){
	result = result.subtract(h);
	break;
      }
    }
    result = result.scale(1/result.mag());
    //cout << result << endl;
    return result;
  };

  // -------------------------
  // Recon Hit Position 
  // -------------------------
  auto recon_projection = [&](const std::vector<dd4pod::VectorXYZ>& hit, const dd4pod::VectorXYZ vector ) {
    dd4pod::VectorXYZ result;
    for (const auto& h : hit) {
      result = h.subtract(vector.scale((h.z+14865)/vector.z));// cellid_converter.position(h.cellID);
      break;
    }
    return result;
  };


//    auto Segmentation = [&](const std::vector<dd4pod::CalorimeterHitData>& evt ) {    
   
//      std::vector<int> result;
     
//      for(const auto& h: evt) {
       
//        int detector = decoder->get(h.cellID,index);
       
//        result.push_back(detector);
       
//      }
//      return result;
//    };


//       //Do some calculations
//  auto d1 = d0.Define("Hit Filter",[](auto conts){
   auto d1 = d0.Define("nhitsT1", nhits, {"TaggerTrackerHits"})
     //     .Define("nhitsT2", nhits, {"Tag_2_Track"})
//      .Filter("nhits>3")
//                 //.Filter([=](const std::vector<edm4hep::SimTrackerHitData>& hits) {
//                 //      for (auto h : hits) {
//                 //        auto pos = ROOT::Math::XYZVector(h.position.x,h.position.y,h.position.z);
//                 //        if ((pos.r() > 100.0) && (std::abs(pos.phi()-M_PI/2.0)< M_PI/6)) {
//                 //          return true;
//                 //        }
//                 //      }
//                 //      return false;
//                 //    },
//                 //    {"GEMTrackerHits"})
     .Define("real_position_x",   "TaggerTrackerHits.position.x")
     .Define("real_position_y",   "TaggerTrackerHits.position.y")
     .Define("real_position_z",   "TaggerTrackerHits.position.z")
    //     .Define("real_position",   project_position,      {"TaggerTrackerHits"})
     //     .Define("TaggerTrackerHits_Prime",  "TaggerTrackerHits[TaggerTrackerHits.PDG==11 && TaggerTrackerHits.generatorStatus==3]")
     .Define("real_vector",     real_vector,           {"TaggerTrackerHits"})
     .Define("recon_position",  recon_position,        {"TaggerTrackerHits"})
     .Define("recon_vector",    recon_vector,          {"recon_position"})
     .Define("recon_projection",recon_projection,      {"recon_position","recon_vector"})
//      .Define("real_vector",     "Tag_1_Track.momentum.scale(1/Tag_1_Track.momentum.mag())")
     //     .Define("fit_position",    fit_position,      {"Tag_1_Track"})
     .Define("detectorTag",  detectorID, {"TaggerTrackerHits"})
     .Define("layerTag",     layer,      {"TaggerTrackerHits"})
     .Define("XpixTag",      xpix,       {"TaggerTrackerHits"})
     .Define("YpixTag",      ypix,       {"TaggerTrackerHits"})
    //      .Define("layerTag2",   layer,      {"Tag_2_Track"})
//      .Define("XpixTag2",    xpix,       {"Tag_2_Track"})
//      .Define("YpixTag2",    ypix,       {"Tag_2_Track"})
     .Define("nhits2", nhits2,      {"MCParticles"})
     .Define("vertex", beamVertex , {"MCParticles"})
     .Define("beamE",  beamVector , {"MCParticles"})
     .Define("scatteredE", electronVector , {"MCParticles"})
     .Define("ex", "scatteredE.px()")
     .Define("ey", "scatteredE.py()")
     .Define("ez", "scatteredE.pz()")
     //     .Define("eTheta", "scatteredE.Theta()")
     .Define("eTheta", "acos(scatteredE.Vect().Unit().Dot(-beamE.Vect().Unit()))")
     //     .Define("ePhi",   "scatteredE.Phi()")
     .Define("ePhi",   "scatteredE.BoostToCM(beamE).Phi()")
     .Define("eE", "scatteredE.energy()")
     .Define("pseudorapidity", "scatteredE.eta()")
     .Define("scatteredV", "beamE-scatteredE")
//      .Define("scatteredV", scatterVector , {"scatteredE"})
     .Define("qx", "scatteredV.px()")
     .Define("qy", "scatteredV.py()")
     .Define("qz", "scatteredV.pz()")
     .Define("qTheta", "scatteredV.Theta()")
     .Define("qPhi",   "scatteredV.Phi()")
     .Define("qE", "scatteredV.energy()")
     .Define("Q2", "-scatteredV.M2()")
     .Define("logQ2", "log10(Q2)")

     .Define("AnglePix11","atan(TaggerTrackerHits.momentum.x/TaggerTrackerHits.momentum.z)[layerTag==0]")
     .Define("AnglePix12","atan(TaggerTrackerHits.momentum.x/TaggerTrackerHits.momentum.z)[layerTag==1]")
   //   .Define("AnglePix21","atan(Tag_2_Track.momentum.x/Tag_2_Track.momentum.z)[layerTag2==0]")
//      .Define("AnglePix22","atan(Tag_2_Track.momentum.x/Tag_2_Track.momentum.z)[layerTag2==1]")

     .Define("Xpix11","XpixTag[layerTag==0&&detectorTag==195]")
     .Define("Ypix11","YpixTag[layerTag==0&&detectorTag==195]")
     .Define("Npix11","Xpix11.size()")
     .Define("Xpix12","XpixTag[layerTag==1&&detectorTag==195]")
     .Define("Ypix12","YpixTag[layerTag==1&&detectorTag==195]")
     .Define("Npix12","Xpix12.size()")
     .Define("Xpix13","XpixTag[layerTag==2&&detectorTag==195]")
     .Define("Ypix13","YpixTag[layerTag==2&&detectorTag==195]")
    .Define("Npix13","Xpix13.size()");





  
   ROOT::RDF::RSnapshotOptions opts;
   opts.fMode = "UPDATE";

   auto d2 = d1.Filter("(Npix11!=0 && Npix12!=0)");// && Npix13!=0)");
   d1.Snapshot("input",outname,{"nhitsT1","nhits2","vertex","Q2","logQ2","ex","ey","ez","eTheta","ePhi","eE","qx","qy","qz","qTheta","qPhi","qE","pseudorapidity"});
   //d2.Snapshot("detector1",outname,{"Xpix11","Ypix11","Npix11","Xpix12","Ypix12","Npix12","Xpix13","Ypix13","Npix13","Q2","logQ2","eE","qTheta","qPhi","eTheta","ePhi","pseudorapidity"},opts);
   //   d1.Filter("B2BPass!=0")
   //  .Snapshot("detectors",outname,{"nhitsT1","layerTag","vertex","Q2","logQ2","eE","qTheta","qPhi","eTheta","ePhi","pseudorapidity","real_position_x","real_position_y","real_position_z","real_vector"},opts);
   d2.Snapshot("detector1",outname,{"vertex","Xpix11","Ypix11","Npix11","Xpix12","Ypix12","Npix12","Xpix13","Ypix13","Npix13","Q2","logQ2","eE","qTheta","qPhi","eTheta","ePhi","pseudorapidity","real_position_x","real_position_y","real_position_z","real_vector"},opts);


}

//--------------------------------------------------------------------------------------------
