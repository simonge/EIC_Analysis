#include <edm4hep/utils/vector_utils.h>
#include <edm4hep/MCParticle.h>
#include <edm4eic/ReconstructedParticle.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TFile.h>

// Define aliases for the data types 
using MCP = edm4hep::MCParticleData;
using RecoP = edm4eic::ReconstructedParticleData;


// Define function to vectorize the edm4hep::utils methods
template <typename T>
auto getEta = [](ROOT::VecOps::RVec<T> momenta) {
  return ROOT::VecOps::Map(momenta, [](const T& p) { return edm4hep::utils::eta(p.momentum); });
};

template <typename T>
auto getPhi = [](ROOT::VecOps::RVec<T> momenta) {
  return ROOT::VecOps::Map(momenta, [](const T& p) { return edm4hep::utils::angleAzimuthal(p.momentum); });
};

// Define the function to perform the efficiency analysis
void EfficiencyAnalysisRDF(TString infile="pythia8NCDIS_18x275_minQ2=10_beamEffects_xAngle=-0.025_hiDiv_5.0002.eicrecon.tree.edm4eic.root"){
   
  // Set up input file 
  ROOT::RDataFrame df("events", infile);

  // Define new dataframe node with additional columns
  auto df1 =  df.Define("statusFilter",  "MCParticles.generatorStatus == 1"    )
                .Define("absPDG",        "abs(MCParticles.PDG)"                )
                .Define("pdgFilter",     "absPDG == 11 || absPDG == 13 || absPDG == 211 || absPDG == 321 || absPDG == 2212")
                .Define("particleFilter","statusFilter && pdgFilter"           )
                .Define("filtMCParts",   "MCParticles[particleFilter]"         )
                .Define("assoFilter",    "Take(particleFilter,_ReconstructedChargedParticleAssociations_sim.index)") // Incase any of the associated particles happen to not be charged
                .Define("assoMCParts",   "Take(MCParticles,_ReconstructedChargedParticleAssociations_sim.index)[assoFilter]")
                .Define("assoRecParts",  "Take(ReconstructedChargedParticles,_ReconstructedChargedParticleAssociations_rec.index)[assoFilter]")
                .Define("filtMCEta",     getEta<MCP>   , {"filtMCParts"} )
                .Define("filtMCPhi",     getPhi<MCP>   , {"filtMCParts"} )
                .Define("accoMCEta",     getEta<MCP>   , {"assoMCParts"} )
                .Define("accoMCPhi",     getPhi<MCP>   , {"assoMCParts"} )
                .Define("assoRecEta",    getEta<RecoP> , {"assoRecParts"})
                .Define("assoRecPhi",    getPhi<RecoP> , {"assoRecParts"})
                .Define("deltaR",        "ROOT::VecOps::DeltaR(assoRecEta, accoMCEta, assoRecPhi, accoMCPhi)");

  // Define histograms
  auto partEta                = df1.Histo1D({"partEta","Eta of Thrown Charged Particles;Eta",100,-5.,5.},"filtMCEta");
  auto matchedPartEta         = df1.Histo1D({"matchedPartEta","Eta of Thrown Charged Particles That Have Matching Track",100,-5.,5.},"accoMCEta");
  auto matchedPartTrackDeltaR = df1.Histo1D({"matchedPartTrackDeltaR","Delta R Between Matching Thrown and Reconstructed Charged Particle",5000,0.,5.},"deltaR");

  // Write histograms to file
  TFile *ofile = TFile::Open("EfficiencyAnalysis_Out_RDF.root","RECREATE");

  // Booked Define and Histo1D lazy actions are only performed here
  partEta->Write();
  matchedPartEta->Write();
  matchedPartTrackDeltaR->Write();
      
  ofile->Close(); // Close output file
}
