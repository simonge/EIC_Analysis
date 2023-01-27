R__LOAD_LIBRARY(libfmt.so)
#include "fmt/core.h"

#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/VectorUtil.h"
#include "TDecompSVD.h"
#include "TMatrixD.h"
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
#include "Math/Vector3D.h" 
#include "TLinearFitter.h"

#include <tuple>
#include <vector>
#include <cmath>

#include "DD4hep/Detector.h"
#include "DDRec/CellIDPositionConverter.h"
 
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
using namespace ROOT::Math;

struct track{  
  XYZVector pos;
  XYZVector vec;
  double    chi2;
  ROOT::VecOps::RVec<int> clusters;
};

using RVecT = ROOT::VecOps::RVec<track>;

  
   //int beamID = 4;
   int beamID = 3;
   int simID  = 1;
   //int simID  = 1;   
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
  fitPoints(ulong nLayers,ulong nModules,ROOT::VecOps::RVec<double> weight): maxLayer(nLayers), maxModules(nModules) {SetLayerWeights(weight); }
  
  RVecT operator()(const std::vector<XYZVector>& positions, const ROOT::VecOps::RVec<int>& moduleID, const ROOT::VecOps::RVec<int>& layerID) {
    
    RVecT outTracks;
    if(positions.size()<maxLayer) return outTracks;
//     if(positions.size()>maxLayer*3) return outTracks;

    ROOT::VecOps::RVec<int>  indeces  (positions.size());
    for(ulong i = 0; i<indeces.size(); i++)
      indeces[i] = i;

    
    auto rPositions = ROOT::VecOps::RVec<XYZVector>(positions);
    
    for(ulong module=1; module<=maxModules; module++){
      auto modHits = rPositions[moduleID==module];
      if(modHits.size()<maxLayer) continue;
      if(modHits.size()>50) continue;
      auto modHitLay  = layerID[moduleID==module];
      auto modIndeces = indeces[moduleID==module];
      
      bool tryFit = true;
      ROOT::VecOps::RVec<ROOT::VecOps::RVec<XYZVector>> layerHits;
      ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>>       layerIndeces;
      for(ulong layer=0; layer<maxLayer; layer++){
	auto layHit = modHits[modHitLay==layer];
	auto layIndeces = modIndeces[modHitLay==layer];
	if(!layHit.size()){
	  tryFit=false;
	  break;
	}
	layerHits.push_back(layHit);
	layerIndeces.push_back(layIndeces);
      }      
      if(!tryFit) continue;

      ROOT::VecOps::RVec<double> x(maxLayer,0);
      ROOT::VecOps::RVec<double> y(maxLayer,0);
      ROOT::VecOps::RVec<double> z(maxLayer,0);
      ROOT::VecOps::RVec<int> index(maxLayer,0);

      //Please try and neaten me up
      TMatrixD ma(3,1);
      TMatrixD mb(3,1);
      TMatrixD mc(3,1);
      TMatrixD md(3,1);
                       
      //Need to somehow generalise number of layers... recursive is ugly?
      int i=0;
      for(auto hit0: layerHits[0]){
	x[0] = hit0.x()*layerWeights[0];
	y[0] = hit0.y()*layerWeights[0];
	z[0] = hit0.z()*layerWeights[0];
	index[0] = layerIndeces[0][i++];
	int j=0;
	for(auto hit1: layerHits[1]){
	  x[1] = hit1.x()*layerWeights[1];
	  y[1] = hit1.y()*layerWeights[1];
	  z[1] = hit1.z()*layerWeights[1];
	  index[1] = layerIndeces[1][j++];
	  int k=0;
	  for(auto hit2: layerHits[2]){
	    x[2] = hit2.x()*layerWeights[2];
	    y[2] = hit2.y()*layerWeights[2];
	    z[2] = hit2.z()*layerWeights[2];
	    index[2] = layerIndeces[2][k++];
	    int l=0;
	    for(auto hit3: layerHits[3]){
	      x[3] = hit3.x()*layerWeights[3];
	      y[3] = hit3.y()*layerWeights[3];
	      z[3] = hit3.z()*layerWeights[3];
	      index[3] = layerIndeces[3][l++];

	      XYZVector outPos = XYZVector(Mean(x),Mean(y),Mean(z))/meanWeight;


	      XYZPoint((hit0-outPos)*layerWeights[0]).GetCoordinates(ma.GetMatrixArray());
	      XYZPoint((hit1-outPos)*layerWeights[1]).GetCoordinates(mb.GetMatrixArray());
	      XYZPoint((hit2-outPos)*layerWeights[2]).GetCoordinates(mc.GetMatrixArray());
	      XYZPoint((hit3-outPos)*layerWeights[3]).GetCoordinates(md.GetMatrixArray());


	      TMatrixD vecMatrix(3,4);
	      vecMatrix.SetSub(0,0,ma);
	      vecMatrix.SetSub(0,1,mb);
	      vecMatrix.SetSub(0,2,mc);
	      vecMatrix.SetSub(0,3,md);
	      
	      TDecompSVD decomp(vecMatrix.T());
	      
	      decomp.Decompose();

	      
	      auto decompVec = decomp.GetV().GetMatrixArray();
	      
	      auto varMatrix = vecMatrix*(decomp.GetV());
	      //varMatrix.Print();
	      auto subMat = varMatrix.GetSub(0,3,1,2);
// 	      subMat.Print();
	      auto subAsArray  = subMat.GetMatrixArray();
	      ROOT::VecOps::RVec<double> subAsVector(subAsArray,subAsArray+8);
	      double outChi2 = Sum(subAsVector*subAsVector)/8;
    

// 	      lf->AssignData(maxLayer, 2, &v[0], &z[0]);  
// 	      lf->Eval();
	      TVectorD params;
	      TVectorD errors;
// 	      lf->GetParameters(params);
// 	      params.Print();
// 	      lf->GetErrors(errors);
 	      //double outChi2=1;
// 	      double outChi2=lf->GetChisquare();
	      //	      XYZVector outVec = XYZVector(params[0],params[1],params[0]);
	      XYZVector outVec = XYZVector(decompVec[0],decompVec[3],decompVec[6]);
// 	      std::cout << outPos << std::endl;
// 	      std::cout << outVec.Unit() << std::endl;
// 	      std::cout << outChi2 << std::endl<< std::endl;	      
	      track outTrack = {outPos,outVec.Unit(),outChi2,index};
	      outTracks.push_back(outTrack);


	    }
	  }
	}
      }
    }
    
    return outTracks;

  };

  
  
  void SetMaxLayers(ulong n){
    maxLayer = n;
  }
  void SetLayerWeights(ROOT::VecOps::RVec<double> weights){
    layerWeights = weights;
    totalWeight  = Sum(layerWeights);
    meanWeight   = Mean(layerWeights);
  }
    
private: 
  ulong maxModules;
  ulong maxLayer;
  ROOT::VecOps::RVec<double> layerWeights = {1,0.8,0.6,0.4};
  double totalWeight = Sum(layerWeights);
  double meanWeight  = Mean(layerWeights);
};
  

// Particle definitions and frame names.
struct partDetails{
  std::string fName;
  int pdg;
  int genStatus;
};


//-----------------------------------------------------------------------------------------
// Main Analysis Call
//-----------------------------------------------------------------------------------------
//   void ProcessTaggerG4(TString inName          = "/scratch/EIC/G4out/derek/x_5_100_16.root",
// 		       TString outName         = "/scratch/EIC/Analysis/temp.root",
// 		       std::string compactName = "/home/simon/geant4/eic/ip6/eic_ip6.xml"){

// void ProcessTaggerG4(TString inName          = "/scratch/EIC/G4out/qr_18x275_beam_out_1.edm4hep.root",
// 		     TString outName         = "/scratch/EIC/Analysis/tempSmall.root",
// 		     std::string compactName = "/home/simon/EIC/epic/epic_ip6.xml"){

// void ProcessTaggerG4Clusters(TString inName          = "/scratch/EIC/G4out/qr_18x275_beam_out_*.edm4hep.root",
// 		     TString outName         = "/scratch/EIC/Analysis/tempClusterQR.root",
// 		     std::string compactName = "/home/simon/EIC/epic/epic_ip6.xml"){
 
// void ProcessTaggerG4(TString inName          = "/scratch/EIC/G4out/qr_18x275_beam_FrontWindow_*.edm4hep.root",
// 		     TString outName         = "/scratch/EIC/Analysis/tempFrontWindow.root",
// 		     std::string compactName = "/home/simon/EIC/epic/epic_ip6.xml"){

//   void ProcessTaggerG4(TString inName          = "/scratch/EIC/G4out/lgen_18x275_beam_out_*.edm4hep.root",
// 		       TString outName         = "/scratch/EIC/Analysis/tempBrems.root",
// 		       std::string compactName = "/home/simon/EIC/epic/epic_ip6.xml"){

//   void ProcessTaggerG4Clusters(TString inName          = "/scratch/EIC/G4out/lgen_18x275_beam_out_1.edm4hep.root",
// 		       TString outName         = "/scratch/EIC/Analysis/clusterSmall.root",
// 		       std::string compactName = "/home/simon/EIC/epic/epic_ip6.xml"){

  void ProcessTaggerG4Clusters(TString inName          = "/scratch/EIC/G4out/gen_uni_FrontWindow_*.edm4hep.root",
		       TString outName         = "/scratch/EIC/Analysis/clusterLarge.root",
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
  //d0.Range(0,10000); // Limit events to analyse 

  // -------------------------
  // Get the DD4hep instance
  // Load the compact XML file
  // Initialize the position converter tool
  dd4hep::Detector& detector = dd4hep::Detector::getInstance();
  detector.fromCompact(compactName);
  dd4hep::rec::CellIDPositionConverter cellid_converter(detector);
  // -------------------------

  //--------------------------------------------------------------------------------------------
  // Lambda Functions
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
    
    std::vector<XYZVector> vectors;

    for (const auto& h : hits) {
      XYZVector result(h.momentum.x,h.momentum.y,h.momentum.z);
      result *= 1/result.R();
      vectors.push_back(result);
    }
    return vectors;
  };

  // -------------------------
  // Real Hit Positions
  // -------------------------
  auto real_position = [&](const std::vector<edm4hep::SimTrackerHitData>& hits) {
    
    std::vector<XYZVector> positions;

    for (const auto& h : hits) {
      XYZVector result(h.position.x,h.position.y,h.position.z);
      positions.push_back(result);
    }
    return positions;
  };

  // -------------------------
  // Cell Hit Positions
  // -------------------------
  auto cell_position = [&](const std::vector<edm4hep::SimTrackerHitData>& hits) {
    
    std::vector<XYZVector> positions;

    for (const auto& h : hits) {
      auto pos1 = cellid_converter.position(h.cellID);
      XYZVector result(pos1.x()*10,pos1.y()*10,pos1.z()*10);
      positions.push_back(result);
    }
    return positions;
  };

  // -------------------------
  // Cell Hit Positions
  // -------------------------
  auto cell_vector = [&](const std::vector<XYZVector>& hits, const ROOT::VecOps::RVec<int>& layerID, const  ROOT::VecOps::RVec<int>& moduleID) {
  //  auto cell_vector = [&](const std::vector<XYZVector>& hits) {
    
    std::vector<XYZVector> vectors;
    XYZVector iPos;
    int module;
    
    for(uint i=0; i<hits.size(); i++){
      if(i==0){
	module=moduleID[i];
	iPos=hits[i];      
	if(layerID[i]!=0) break;
      }
      if(i>0 && module==moduleID[i] && layerID[i]>=1) vectors.push_back((hits[i]-iPos).Unit());
    }
    
    while(vectors.size()<hits.size()){
      vectors.push_back(XYZVector(1,0,0));
    }

    return vectors;
  };

  // -------------------------
  // Vector XZ origin cut
  // -------------------------
//   auto vector_cut = [&](const std::vector<XYZVector>& positions, const std::vector<XYZVector>& vectors) {
    
//     if(positions.size()<=0) return XYZVector();

//     auto vec = vectors[0];
//     auto pos = positions[0];

//     float con = pos.x()/vec.x();

//     auto   xPoint = pos-con*vec;
//     return xPoint;

//   };

  auto vector_cut = [&](const  std::vector<XYZVector>& positions, const  std::vector<XYZVector>& vectors) {
    
    std::vector<XYZVector> points;

    for(uint i=0; i<positions.size(); i++){
      
      float con = positions[i].x()/vectors[i].x();

      points.push_back(positions[i]-con*vectors[i]);

    }
    return points;

  };

//   auto vector_cut2 = [&](const std::vector<XYZVector>& positions, const XYZVector& vectors) {
    
//     std::vector<XYZVector> points;
      
//     float con = positions[0].x()/vectors.x();

//     points.push_back(positions[0]-con*vectors);

//     return points;

//   };
  
  auto vector_filter = [&](const std::vector<XYZVector>& positions) {
    
    ROOT::VecOps::RVec<bool> good;

    for(auto pos: positions){
      if(pos.z()<-8000&&pos.z()>-15000&&abs(pos.y())<100) good.push_back(true);
      else good.push_back(false);

    }
    return good;

  };
  
  //Labels hits with cluster IDs
  auto cluster_hits = [&](const ROOT::VecOps::RVec<int> module, const ROOT::VecOps::RVec<int> layer,const ROOT::VecOps::RVec<int> x,const ROOT::VecOps::RVec<int> y,const ROOT::VecOps::RVec<float> t,const ROOT::VecOps::RVec<float> e) {
    
    ROOT::VecOps::RVec<bool> avaliable(module.size(),1);
    ROOT::VecOps::RVec<int>  clusterNo(module.size(),0);
    ROOT::VecOps::RVec<int>  indeces  (module.size());
    for(ulong i = 0; i<indeces.size(); i++)
      indeces[i] = i;
   
    int seedNumber = 0;

    while(Any(avaliable)){

      //Cluster seed
      auto maxIndex = ArgMax(e*avaliable);

      avaliable[maxIndex] = 0;
      clusterNo[maxIndex] = seedNumber;

      ROOT::VecOps::RVec<ulong> indexList = {maxIndex};
      
      while(indexList.size()){

	auto index = indexList[0];
	auto filter = avaliable*(module==module[index])*(layer==layer[index])*(abs(x-x[index])<=1)*(abs(y-y[index])<=1)*(abs(t-t[index])<1);
	avaliable = avaliable*(!filter);
	clusterNo += filter*seedNumber;
	indexList = Concatenate(indexList,indeces[filter]);
	
	indexList.erase(indexList.begin());
      }
      seedNumber++;
    }    

//     std::cout << seedNumber << " " << Max(clusterNo) << " " << clusterNo.size() << std::endl<< std::endl;

    return clusterNo;

  };

  // -------------------------
  // Cluster Hit Positions
  // -------------------------
  auto cluster_position = [&](const ROOT::VecOps::RVec<int>  clusterNo, const ROOT::VecOps::RVec<float> e, const std::vector<XYZVector> positions) {
    
    std::vector<XYZVector> clusterPos;
    
    for(int index=0; index<=Max(clusterNo); index++){
      auto filter = (index==clusterNo);
      auto eSum            = Sum(e[filter]);
      XYZVector sumPositions;
      for(ulong i=0; i<filter.size(); i++){
	if(filter[i]==0) continue;
	sumPositions += positions[i]*e[i];
      }
      
      XYZVector outvec = sumPositions/eSum;
      clusterPos.push_back(outvec);
    }
    return clusterPos;
  };

  // -------------------------
  // Cluster detID
  // -------------------------
  auto cluster_detID = [&](const ROOT::VecOps::RVec<int>  clusterNo, const ROOT::VecOps::RVec<int> detID) {
    
    ROOT::VecOps::RVec<int> clusterID;
    
    for(int index=0; index<=Max(clusterNo); index++){
      auto filter = (index==clusterNo);
      auto ID     = (detID[filter])[0];
      clusterID.push_back(ID);
    }
    
    return clusterID;
  };

  // -------------------------
  // Cluster maxID
  // -------------------------
  auto cluster_maxID = [&](const ROOT::VecOps::RVec<int>  clusterNo, const ROOT::VecOps::RVec<int> detID) {
    
    ROOT::VecOps::RVec<int> clusterID;
    
    for(int index=0; index<=Max(clusterNo); index++){
      auto filter = (index==clusterNo);
      auto ID     = Max(detID[filter]);
      clusterID.push_back(ID);
    }    
    return clusterID;
  };

  // -------------------------
  // Cluster minID
  // -------------------------
  auto cluster_minID = [&](const ROOT::VecOps::RVec<int>  clusterNo, const ROOT::VecOps::RVec<int> detID) {
    
    ROOT::VecOps::RVec<int> clusterID;
    
    for(int index=0; index<=Max(clusterNo); index++){
      auto filter = (index==clusterNo);
      auto ID     = Min(detID[filter]);
      clusterID.push_back(ID);
    }    
    return clusterID;
  };

  // -------------------------
  // Cluster size
  // -------------------------
  auto cluster_size = [&](const ROOT::VecOps::RVec<int>  clusterNo) {
    
    ROOT::VecOps::RVec<int> clusterID;
    
    for(int index=0; index<=Max(clusterNo); index++){
      auto filter = (index==clusterNo);
      auto ID     = Sum(filter);
      clusterID.push_back(ID);
    }    
    return clusterID;
  };

  // -------------------------
  // track minID
  // -------------------------
  auto track_minID = [&](const std::vector<ROOT::VecOps::RVec<int>>  trackClusters, const ROOT::VecOps::RVec<int> clusterLim) {
    
    ROOT::VecOps::RVec<int> lim;
    
    for(auto clusters:trackClusters){
      auto ID     = Min(Take(clusterLim,clusters));
      lim.push_back(ID);
    }    
    return lim;
  };
  // -------------------------
  // track maxID
  // -------------------------
  auto track_maxID = [&](const std::vector<ROOT::VecOps::RVec<int>>  trackClusters, const ROOT::VecOps::RVec<int> clusterLim) {
    
    ROOT::VecOps::RVec<int> lim;
    
    for(auto clusters:trackClusters){
      auto ID     = Max(Take(clusterLim,clusters));
      lim.push_back(ID);
    }    
    return lim;
  };
  // -------------------------
  // track detID
  // -------------------------
  auto track_detID = [&](const std::vector<ROOT::VecOps::RVec<int>>  trackClusters, const ROOT::VecOps::RVec<int> clusterMod) {
    
    ROOT::VecOps::RVec<int> mod;
    
    for(auto clusters:trackClusters){
      auto ID     = Take(clusterMod,clusters)[0];
      mod.push_back(ID);
    }    
    return mod;
  };

  //----------------------------------------------------------------------------
  // Fill Dataframes
  //----------------------------------------------------------------------------
  

   //-----------------------------------------
   // Hits positions by real/cell positions
   //-----------------------------------------
   auto d1 = d0.Define("nHits", "TaggerTrackerHits.size()")
     .Define("real_position",     real_position,           {"TaggerTrackerHits"})
     .Define("cell_position",     cell_position,           {"TaggerTrackerHits"})
     .Define("real_time",         "TaggerTrackerHits.time")
     .Define("real_EDep",         "TaggerTrackerHits.EDep")
     .Define("real_vector",       real_vector,             {"TaggerTrackerHits"})
     .Define("real_cut",          vector_cut,              {"real_position","real_vector"})
     .Define("vector_filter",     vector_filter,           {"real_cut"});
   
   //-----------------------------------------
   // Chosen Particles to save
   //-----------------------------------------
   std::vector<partDetails> parts = {{"beamElectron",11,beamID},{"beamProton",2212,beamID},{"scatteredElectron",11,simID}};
   std::vector<std::string> Part_Vec;
   for(auto part: parts){
     std::string colName = part.fName;
     Part_Vec.push_back(colName);
     d1 = d1.Define(colName,getParticle(part.genStatus,part.pdg),{"MCParticles"});
   }

   //-----------------------------------------
   // Hit detector IDs
   //-----------------------------------------
   auto ids = detector.readout("TaggerTrackerHits").idSpec().fields();
   std::vector<std::string> ID_Vec;
   for(auto id: ids){
     std::string colName = id.first+"ID";
     ID_Vec.push_back(colName);
     d1 = d1.Define(colName,getSubID(id.first,detector),{"TaggerTrackerHits"});
   }
   
   //-----------------------------------------
   // Cluster Hits
   //-----------------------------------------

   d1 = d1.Define("hitClusterNo",cluster_hits,{"moduleID","layerID","xID","yID","real_time","real_EDep"})
     .Define("clusterPos",cluster_position,{"hitClusterNo","real_EDep","cell_position"})
     .Define("clusterMod",cluster_detID,{"hitClusterNo","moduleID"})
     .Define("clusterLay",cluster_detID,{"hitClusterNo","layerID"})
     .Define("clusterSize",cluster_size,{"hitClusterNo"})
     .Define("clusterMinX",cluster_minID,{"hitClusterNo","xID"})
     .Define("clusterMinY",cluster_minID,{"hitClusterNo","yID"})
     .Define("clusterMaxX",cluster_maxID,{"hitClusterNo","xID"})
     .Define("clusterMaxY",cluster_maxID,{"hitClusterNo","yID"})
     .Define("NHits","TaggerTrackerHits.size()")
     .Define("NClusters","Max(hitClusterNo)");

   std::vector<std::string> Cluster_Vec = {"NClusters","clusterPos","clusterSize","clusterMod","clusterLay"};

   d1 = d1.Define("fit_temp", fitPoints(4,2,{1,0.8,0.7,0.6}) , {"clusterPos","clusterMod","clusterLay"})
     .Define("fit_position", [](RVecT tracks){std::vector<XYZVector> outvec; for(auto track:tracks) outvec.push_back(track.pos);  return outvec;},{"fit_temp"})
     .Define("fit_vector",   [](RVecT tracks){std::vector<XYZVector> outvec; for(auto track:tracks) outvec.push_back(track.vec);  return outvec;},{"fit_temp"})
     .Define("fit_chi2",     [](RVecT tracks){std::vector<double>    outvec; for(auto track:tracks) outvec.push_back(track.chi2); return outvec;},{"fit_temp"})
     .Define("fit_clusters", [](RVecT tracks){std::vector<ROOT::VecOps::RVec<int>>    outvec; for(auto track:tracks) outvec.push_back(track.clusters); return outvec;},{"fit_temp"})
     .Define("fit_minX", track_minID ,{"fit_clusters","clusterMinX"})
     .Define("fit_maxX", track_maxID ,{"fit_clusters","clusterMaxX"})
     .Define("fit_minY", track_minID ,{"fit_clusters","clusterMinY"})
     .Define("fit_maxY", track_maxID ,{"fit_clusters","clusterMaxY"})
     .Define("fit_module", track_detID ,{"fit_clusters","clusterMod"})
     .Define("NTracks",    "fit_position.size()")
     .Define("fit_cut",      vector_cut,       {"fit_position","fit_vector"});

   std::vector<std::string> Track_Vec = {"fit_position","fit_vector","fit_chi2","fit_cut","fit_minX","fit_maxX","fit_minY","fit_maxY","fit_module","NTracks"};


//    d1.Snapshot("clusters",outName,{"real_EDep","moduleID","layerID","xID","yID","real_time","hitClusterNo","NClusters","NHits","clusterPos","fit_position","fit_vector","fit_chi2","fit_cut"});


   
//    d1 = d1.Define("CapHit",getSubID("system",detector,"Cap_Track"))
//      .Define("CapPass","CapHit.size()");
   d1 = d1.Define("cell_vector",       cell_vector,             {"cell_position","layerID","moduleID"})
     .Define("cell_cut",          vector_cut,              {"cell_position","cell_vector"})
     .Define("iFilter","(TMath::Pi()-scatteredElectron.Theta())<0.0105")
     .Define("vertex", beamVertex , {"MCParticles"})
     .Define("nParticles",        "MCParticles.size()")
     .Define("time", beamTime , {"MCParticles"})
     .Define("eE", "scatteredElectron.energy()")
     .Define("thetaE", "acos((scatteredElectron.Vect().Unit()).Dot(beamElectron.Vect().Unit()))")
     //     .Define("phiE", "acos((((beamElectron.Vect()).Cross(XYZVector(0,1,0))).Unit()).Dot((scatteredElectron.Vect()).Cross(XYZVector(0,1,0)).Unit()))")
     .Define("pseudorapidity", "scatteredElectron.eta()")
     .Define("scatteredV", "beamElectron-scatteredElectron")
     .Define("thetaV", "scatteredV.theta()")
     .Define("phiV", "scatteredV.phi()")
     .Define("qE", "scatteredV.energy()")
     .Define("Q2", "2*beamElectron.energy()*eE*(1-cos(TMath::Pi()-scatteredElectron.theta()))")
     .Define("logQ2", "log10(Q2)")
     .Define("Q2_2", "2*beamElectron.energy()*eE*(1-cos(thetaE))")
     .Define("logQ2_2", "log10(Q2_2)")
     .Define("Q2_3", "-scatteredV.M2()")
     .Define("logQ2_3", "log10(Q2_3)")
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

   // Specific variables for training network faster?

//    d1 = d1.Define("real_cut_x")
//      .Define("cosPhi","cos(phiV)")
//      .Define("cosPhi","sin(phiV)");

   //[&](const std::vector<XYZVector>& vec){}

   
    ROOT::RDF::RSnapshotOptions opts;
    opts.fMode = "UPDATE";
    
//     d1.Snapshot("ids",outName,ID_Vec);
//     d1.Snapshot("parts",outName,Part_Vec,opts);

   

    std::vector<std::string> Out_Vec = {"vertex","nParticles","nHits","real_position","cell_position","cell_vector","cell_cut","real_vector","iFilter","real_cut","vector_filter","thetaV","thetaE","phiV","qE","eE","logQ2","logQ2_2","logQ2_3","Tag1_1","Tag1_2","Tag1_3","Tag1_4","Tag2_1","Tag2_2","Tag2_3","Tag2_4","x11","y11","x12","y12","x13","y13","x14","y14","x21","y21","x22","y22","x23","y23","x24","y24"};//,"fit_vector","fit_chi2"};

    for(auto a:ID_Vec) Out_Vec.push_back(a);
    for(auto a:Part_Vec) Out_Vec.push_back(a);
    for(auto a:Cluster_Vec) Out_Vec.push_back(a);
    for(auto a:Track_Vec) Out_Vec.push_back(a);


    d1.Snapshot("temp",outName,Out_Vec);
    //    d1.Snapshot("temp",outName,Out_Vec,opts);
//    //Filter("(Tag1_4||Tag2_4)&&iFilter").

   //Hits distribution/layer
   //Acceptance 
  

//    auto d2 = d1.Filter("(Npix11!=0 && Npix12!=0)");// && Npix13!=0)");
//    d1.Snapshot("input",outName,{"nhitsT1","nhits2","vertex","Q2","logQ2","ex","ey","ez","eTheta","ePhi","eE","qx","qy","qz","qTheta","qPhi","qE","pseudorapidity"});
//    //d2.Snapshot("detector1",outname,{"x11","y11","Npix11","x12","y12","Npix12","x13","y13","Npix13","Q2","logQ2","eE","qTheta","qPhi","eTheta","ePhi","pseudorapidity"},opts);
//    //   d1.Filter("B2BPass!=0")
//    //  .Snapshot("detectors",outname,{"nhitsT1","layerID","vertex","Q2","logQ2","eE","qTheta","qPhi","eTheta","ePhi","pseudorapidity","cell_position_x","cell_position_y","cell_position_z","cell_vector"},opts);
//    d2.Snapshot("detector1",outName,{"vertex","x11","y11","Npix11","x12","y12","Npix12","x13","y13","Npix13","Q2","logQ2","eE","qTheta","qPhi","eTheta","ePhi","pseudorapidity","cell_position_x","cell_position_y","cell_position_z","cell_vector"},opts);


}

//--------------------------------------------------------------------------------------------
