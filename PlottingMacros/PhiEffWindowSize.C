#include "TString.h"
#include "TCanvas.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include <vector>
#include "TH1.h"
//#include "ROOT.h"


using namespace std;

//std::string tag = "Decor-NormTE-4Layers-Bigger-4STEP4";
//std::string tag = "RealHits4layerclusterFrontWindow";
std::string tag = "RealHits4layertempClusterQR";
std::string tag2 = "RealHits4layertempClusterBrems";
//std::string tag = "FitTrain_uniform";

std::vector<TString> fileNames = {"/scratch/EIC/Results/ML-Out/"+tag+"_ETP.root","/home/simon/Analysis/EIC_Analysis/NeuralNetwork/Reg_point-"+tag+".root","/home/simon/Analysis/EIC_Analysis/NeuralNetwork/Reg_point2-"+tag2+".root"};
//std::vector<TString> fileNames = {"/scratch/EIC/Results/ML-Out/"+tag+"_ETP.root","/home/simon/Analysis/EIC_Analysis/NeuralNetwork/Reg_point-"+tag+".root"};

// struct countRMS{
//   ROOT::RDF::RResultPtr<unsigned long long> counts;
//   ROOT::RDF::RResultPtr<double>  rms;
// };
struct countRMS{
  double counts;
  double rms;
};

struct GetCounts2{
  GetCounts2(THnD* hist): Hist(hist){}

  countRMS operator()(ULong64_t ebin,ULong64_t tbin,ULong64_t pbin){
    Hist->GetAxis(0)->SetRange(ebin+1,ebin+1);
    Hist->GetAxis(1)->SetRange(tbin+1,tbin+1);
    Hist->GetAxis(2)->SetRange(pbin+1,pbin+1);
    auto HistProj = Hist->Projection(3);
    
    auto count = HistProj->GetEntries();
    double rms = 0;
    if(count!=0){
    //   TFitResultPtr r = HistProj->Fit("gaus","S");
//       rms   = r->Parameter(2);
      rms   = HistProj->GetStdDev();
    }
    
    countRMS retVal = {count,rms};
    return retVal; 
  };

private:
  THnD* Hist;

};


void PhiEffWindowSize(){
  //ROOT::EnableImplicitMT();

  TString outNamepng0  = "plots/PSizeEfficiencies-"+tag+"00.png";
  TString outNamepng  = "plots/PSizeEfficiencies-"+tag+".png";
  TString outNamepng2  = "plots/BremEfficiencies-"+tag+".png";

  gStyle->SetStatW(0.3);
  gStyle->SetStatColor(0);
  gStyle->SetStatBorderSize(1);
  gStyle->SetOptStat(0);
  //  gStyle->SetOptStat(1111);
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

  ROOT::RDataFrame dfQR("predict",fileNames[1]);

  auto dfQR_2 = dfQR//.Filter("fit_minchi2<0.002")
    .Define("Phi","phiV*TMath::RadToDeg()")
    .Define("ThetaGen","thetaE*1000")
    .Define("CalPhiC","pPredC*TMath::RadToDeg()")
    .Define("PhiResC","(Phi-CalPhiC)")
    .Define("ThetaReconC","tPredC*1000")
    .Define("ThetaResC","(ThetaGen-ThetaReconC)")
    .Define("EResC","100*(eE-ePredC)/eE");

  ROOT::RDataFrame dfBrem("predict",fileNames[2]);

//   auto dfBrem_2 = dfBrem//.Filter("fit_minchi2<0.002")
//     .Define("Phi","phiV*TMath::RadToDeg()")
//     .Define("ThetaGen","thetaE*1000")
//     .Define("CalPhiC","pPredC*TMath::RadToDeg()")
//     .Define("PhiResC","(Phi-CalPhiC)")
//     .Define("ThetaReconC","tPredC*1000")
//     .Define("ThetaResC","(ThetaGen-ThetaReconC)")
//     .Define("EResC","100*(eE-ePredC)/eE");

  int    eBins     = 30;
  int    pBins     = 30;
  int    tBins     = 30;
  int    resBins   = 400;

  double eMin = 6;
  double eMax = 18;
  double pMin = -180;
  double pMax = 180;
  double tMin = 0.0;
  double tMax = 0.011*1000;
    
  double PResRange = TMath::Pi()*TMath::RadToDeg();
  double TResRange = 0.002*1000;
  double EResRange = 0.05*100;
  
  const std::vector< int >    bins = {eBins,tBins,pBins,resBins};
  const std::vector< double > min  = {eMin,tMin,pMin,-PResRange};
  const std::vector< double > max  = {eMax,tMax,pMax, PResRange};
  
  double tSize = (tMax-tMin)/tBins;
  double pSize = (pMax-pMin)/pBins;
  double eSize = (eMax-eMin)/eBins;
 
  auto Phi4D    = dfQR_2.HistoND({"PResC", ";Electron energy [GeV] ;#theta [mrad]",4, bins,min,max }, {"eE","ThetaGen","Phi","PhiResC"});
  auto Phi4D1   = dfQR_2.Filter("(ymin>(-256*1))&&(ymax<(256*1))").HistoND({"PResC", ";Electron energy [GeV] ;#theta [mrad]",4, bins,min,max }, {"eE","ThetaGen","Phi","PhiResC"});
  auto Phi4D2   = dfQR_2.Filter("(ymin>(-256*2))&&(ymax<(256*2))").HistoND({"PResC", ";Electron energy [GeV] ;#theta [mrad]",4, bins,min,max }, {"eE","ThetaGen","Phi","PhiResC"});
  auto Phi4D3   = dfQR_2.Filter("(ymin>(-256*3))&&(ymax<(256*3))").HistoND({"PResC", ";Electron energy [GeV] ;#theta [mrad]",4, bins,min,max }, {"eE","ThetaGen","Phi","PhiResC"});
  auto Phi4D4   = dfQR_2.Filter("(ymin>(-256*4))&&(ymax<(256*4))").HistoND({"PResC", ";Electron energy [GeV] ;#theta [mrad]",4, bins,min,max }, {"eE","ThetaGen","Phi","PhiResC"});
  auto Phi4D5   = dfQR_2.Filter("(ymin>(-256*5))&&(ymax<(256*5))").HistoND({"PResC", ";Electron energy [GeV] ;#theta [mrad]",4, bins,min,max }, {"eE","ThetaGen","Phi","PhiResC"});
  auto Phi4D6   = dfQR_2.Filter("(ymin>(-256*6))&&(ymax<(256*6))").HistoND({"PResC", ";Electron energy [GeV] ;#theta [mrad]",4, bins,min,max }, {"eE","ThetaGen","Phi","PhiResC"});
 
  auto RPhi4D    = dfQR_2.HistoND({"PResC", ";Electron energy [GeV] ;#theta [mrad]",4, bins,min,max }, {"eE","ThetaGen","Phi","PhiResC"});
  auto RPhi4D1   = dfQR_2.Filter("(ymin>(-256*1))&&(ymax<(256*1))").HistoND({"PResC", ";Electron energy [GeV] ;#theta [mrad]",4, bins,min,max }, {"ePredC","ThetaReconC","CalPhiC","PhiResC"});
  auto RPhi4D2   = dfQR_2.Filter("(ymin>(-256*2))&&(ymax<(256*2))").HistoND({"PResC", ";Electron energy [GeV] ;#theta [mrad]",4, bins,min,max }, {"ePredC","ThetaReconC","CalPhiC","PhiResC"});
  auto RPhi4D3   = dfQR_2.Filter("(ymin>(-256*3))&&(ymax<(256*3))").HistoND({"PResC", ";Electron energy [GeV] ;#theta [mrad]",4, bins,min,max }, {"ePredC","ThetaReconC","CalPhiC","PhiResC"});
  auto RPhi4D4   = dfQR_2.Filter("(ymin>(-256*4))&&(ymax<(256*4))").HistoND({"PResC", ";Electron energy [GeV] ;#theta [mrad]",4, bins,min,max }, {"ePredC","ThetaReconC","CalPhiC","PhiResC"});
  auto RPhi4D5   = dfQR_2.Filter("(ymin>(-256*5))&&(ymax<(256*5))").HistoND({"PResC", ";Electron energy [GeV] ;#theta [mrad]",4, bins,min,max }, {"ePredC","ThetaReconC","CalPhiC","PhiResC"});
  auto RPhi4D6   = dfQR_2.Filter("(ymin>(-256*6))&&(ymax<(256*6))").HistoND({"PResC", ";Electron energy [GeV] ;#theta [mrad]",4, bins,min,max }, {"ePredC","ThetaReconC","CalPhiC","PhiResC"});

  //Brems
  auto dfBrem_2 = dfBrem
    .Define("Phi","phiV*TMath::RadToDeg()")
    .Define("ThetaGen","thetaE*1000")
    .Define("CalPhiC","pPredC*TMath::RadToDeg()")
    .Define("PhiResC","(Phi-CalPhiC)")
    .Define("ThetaReconC","tPredC*1000")
    .Define("ThetaResC","(ThetaGen-ThetaReconC)")
    .Define("EResC","100*(eE-ePredC)/eE");

  auto BremPhi4D    = dfBrem_2.HistoND({"BremPResC", ";Electron energy [GeV] ;#theta [mrad]",4, bins,min,max }, {"ePredC","ThetaReconC","CalPhiC","PhiResC"});
  auto BremPhi4D1   = dfBrem_2.Filter("(ymin>(-256*1))&&(ymax<(256*1))").HistoND({"BremPResC1", ";Electron energy [GeV] ;#theta [mrad]",4, bins,min,max }, {"ePredC","ThetaReconC","CalPhiC","PhiResC"});
  auto BremPhi4D2   = dfBrem_2.Filter("(ymin>(-256*2))&&(ymax<(256*2))").HistoND({"BremPResC2", ";Electron energy [GeV] ;#theta [mrad]",4, bins,min,max }, {"ePredC","ThetaReconC","CalPhiC","PhiResC"});
  auto BremPhi4D3   = dfBrem_2.Filter("(ymin>(-256*3))&&(ymax<(256*3))").HistoND({"BremPResC3", ";Electron energy [GeV] ;#theta [mrad]",4, bins,min,max }, {"ePredC","ThetaReconC","CalPhiC","PhiResC"});
  auto BremPhi4D4   = dfBrem_2.Filter("(ymin>(-256*4))&&(ymax<(256*4))").HistoND({"BremPResC4", ";Electron energy [GeV] ;#theta [mrad]",4, bins,min,max }, {"ePredC","ThetaReconC","CalPhiC","PhiResC"});
  auto BremPhi4D5   = dfBrem_2.Filter("(ymin>(-256*5))&&(ymax<(256*5))").HistoND({"BremPResC5", ";Electron energy [GeV] ;#theta [mrad]",4, bins,min,max }, {"ePredC","ThetaReconC","CalPhiC","PhiResC"});
  auto BremPhi4D6   = dfBrem_2.Filter("(ymin>(-256*6))&&(ymax<(256*6))").HistoND({"BremPResC6", ";Electron energy [GeV] ;#theta [mrad]",4, bins,min,max }, {"ePredC","ThetaReconC","CalPhiC","PhiResC"});
  
  int totalRows = eBins*tBins*pBins;
  ROOT::RDataFrame d(totalRows);

  int slots = d.GetNSlots();

  auto d2 = d.Define("EBin",[&eBins](ULong64_t i){return i%eBins;},{"rdfentry_"})
    .Define("EMin",[&eMin,&eBins,&eSize](ULong64_t i){return eMin+i*eSize+(eSize/2);},{"EBin"})
    .Define("TBin",[&tBins,&eBins](ULong64_t i){return ((ULong64_t)(i/eBins)%tBins);},{"rdfentry_"})
    .Define("TMin",[&tMin,&tSize](ULong64_t i){return tMin+i*tSize+(tSize/2);},{"TBin"})
    .Define("PBin",[&pBins,&tBins,&eBins](ULong64_t i){return ((ULong64_t)(i/(eBins*tBins))%pBins);},{"rdfentry_"})
    .Define("PMin",[&pMin,&pSize](ULong64_t i){return pMin+i*pSize+(pSize/2);},{"PBin"})
    .Define("CountsRMS",GetCounts2(Phi4D.GetPtr()),{"EBin","TBin","PBin"}) //Real
    .Define("Counts","CountsRMS.counts")
    .Define("RMS","CountsRMS.rms")
    .Define("CountsRMS1",GetCounts2(Phi4D1.GetPtr()),{"EBin","TBin","PBin"})
    .Define("Counts1","CountsRMS1.counts")
    .Define("RMS1","CountsRMS1.rms")
    .Define("CountsRMS2",GetCounts2(Phi4D2.GetPtr()),{"EBin","TBin","PBin"})
    .Define("Counts2","CountsRMS2.counts")
    .Define("RMS2","CountsRMS2.rms")
    .Define("CountsRMS3",GetCounts2(Phi4D3.GetPtr()),{"EBin","TBin","PBin"})
    .Define("Counts3","CountsRMS3.counts")
    .Define("RMS3","CountsRMS3.rms")
    .Define("CountsRMS4",GetCounts2(Phi4D4.GetPtr()),{"EBin","TBin","PBin"})
    .Define("Counts4","CountsRMS4.counts")
    .Define("RMS4","CountsRMS4.rms")
    .Define("CountsRMS5",GetCounts2(Phi4D5.GetPtr()),{"EBin","TBin","PBin"})
    .Define("Counts5","CountsRMS5.counts")
    .Define("RMS5","CountsRMS5.rms")
    .Define("CountsRMS6",GetCounts2(Phi4D6.GetPtr()),{"EBin","TBin","PBin"})
    .Define("Counts6","CountsRMS6.counts")
    .Define("RMS6","CountsRMS6.rms")
    .Define("RCountsRMS",GetCounts2(RPhi4D.GetPtr()),{"EBin","TBin","PBin"}) //QR Recon
    .Define("RCounts","RCountsRMS.counts")
    .Define("RRMS","RCountsRMS.rms")
    .Define("RCountsRMS1",GetCounts2(RPhi4D1.GetPtr()),{"EBin","TBin","PBin"})
    .Define("RCounts1","RCountsRMS1.counts")
    .Define("RRMS1","RCountsRMS1.rms")
    .Define("RCountsRMS2",GetCounts2(RPhi4D2.GetPtr()),{"EBin","TBin","PBin"})
    .Define("RCounts2","RCountsRMS2.counts")
    .Define("RRMS2","RCountsRMS2.rms")
    .Define("RCountsRMS3",GetCounts2(RPhi4D3.GetPtr()),{"EBin","TBin","PBin"})
    .Define("RCounts3","RCountsRMS3.counts")
    .Define("RRMS3","RCountsRMS3.rms")
    .Define("RCountsRMS4",GetCounts2(RPhi4D4.GetPtr()),{"EBin","TBin","PBin"})
    .Define("RCounts4","RCountsRMS4.counts")
    .Define("RRMS4","RCountsRMS4.rms")
    .Define("RCountsRMS5",GetCounts2(RPhi4D5.GetPtr()),{"EBin","TBin","PBin"})
    .Define("RCounts5","RCountsRMS5.counts")
    .Define("RRMS5","RCountsRMS5.rms")
    .Define("RCountsRMS6",GetCounts2(RPhi4D6.GetPtr()),{"EBin","TBin","PBin"})
    .Define("RCounts6","RCountsRMS6.counts")
    .Define("RRMS6","RCountsRMS6.rms")
    .Define("BCountsRMS",GetCounts2(BremPhi4D.GetPtr()),{"EBin","TBin","PBin"}) //Brems Recon
    .Define("BCounts","BCountsRMS.counts")
    .Define("BRMS","BCountsRMS.rms")
    .Define("BCountsRMS1",GetCounts2(BremPhi4D1.GetPtr()),{"EBin","TBin","PBin"})
    .Define("BCounts1","BCountsRMS1.counts")
    .Define("BRMS1","BCountsRMS1.rms")
    .Define("BCountsRMS2",GetCounts2(BremPhi4D2.GetPtr()),{"EBin","TBin","PBin"})
    .Define("BCounts2","BCountsRMS2.counts")
    .Define("BRMS2","BCountsRMS2.rms")
    .Define("BCountsRMS3",GetCounts2(BremPhi4D3.GetPtr()),{"EBin","TBin","PBin"})
    .Define("BCounts3","BCountsRMS3.counts")
    .Define("BRMS3","BCountsRMS3.rms")
    .Define("BCountsRMS4",GetCounts2(BremPhi4D4.GetPtr()),{"EBin","TBin","PBin"})
    .Define("BCounts4","BCountsRMS4.counts")
    .Define("BRMS4","BCountsRMS4.rms")
    .Define("BCountsRMS5",GetCounts2(BremPhi4D5.GetPtr()),{"EBin","TBin","PBin"})
    .Define("BCounts5","BCountsRMS5.counts")
    .Define("BRMS5","BCountsRMS5.rms")
    .Define("BCountsRMS6",GetCounts2(BremPhi4D6.GetPtr()),{"EBin","TBin","PBin"})
    .Define("BCounts6","BCountsRMS6.counts")
    .Define("BRMS6","BCountsRMS6.rms")
    .Define("Ratio", "RCounts/(RCounts+BCounts*11.65)")
    .Define("Ratio1","RCounts1/(RCounts1+BCounts1*11.65)")
    .Define("Ratio2","RCounts2/(RCounts2+BCounts2*11.65)")
    .Define("Ratio3","RCounts3/(RCounts3+BCounts3*11.65)")
    .Define("Ratio4","RCounts4/(RCounts4+BCounts4*11.65)")
    .Define("Ratio5","RCounts5/(RCounts5+BCounts5*11.65)")
    .Define("Ratio6","RCounts6/(RCounts6+BCounts6*11.65)");

//   TString countFilter="Counts>10Counts>10";

//   auto TotalCounts  = d2.Histo1D({"TotalEnergyCountsR", ";#Energy [GeV]", eBins, eMin, eMax },"EMin","Counts");
//   auto ReconCounts  = d2.Filter(countFilter+"&&RMS<20").Histo1D({"PhiReconEnergyCounts", ";#Energy [GeV]", eBins, eMin, eMax },"EMin","Counts");
//   auto ReconCounts1 = d2.Filter(countFilter+"&&RMS1<20").Histo1D({"PhiReconEnergyCounts1", ";#Energy [GeV]", eBins, eMin, eMax },"EMin","Counts");
//   auto ReconCounts2 = d2.Filter(countFilter+"&&RMS2<20").Histo1D({"PhiReconEnergyCounts2", ";#Energy [GeV]", eBins, eMin, eMax },"EMin","Counts");
//   auto ReconCounts3 = d2.Filter(countFilter+"&&RMS3<20").Histo1D({"PhiReconEnergyCounts3", ";#Energy [GeV]", eBins, eMin, eMax },"EMin","Counts");
//   auto ReconCounts4 = d2.Filter(countFilter+"&&RMS4<20").Histo1D({"PhiReconEnergyCounts4", ";#Energy [GeV]", eBins, eMin, eMax },"EMin","Counts");
//   auto ReconCounts5 = d2.Filter(countFilter+"&&RMS5<20").Histo1D({"PhiReconEnergyCounts5", ";#Energy [GeV]", eBins, eMin, eMax },"EMin","Counts");
//   auto ReconCounts6 = d2.Filter(countFilter+"&&RMS6<20").Histo1D({"PhiReconEnergyCounts6", ";#Energy [GeV]", eBins, eMin, eMax },"EMin","Counts");

  auto TotalCounts  = d2.Histo1D({"TotalEnergyCountsR", ";Energy [GeV]", eBins, eMin, eMax },"EMin","Counts");
  auto ReconCounts  = d2.Filter("Counts>10&&RMS<30").Histo1D({"PhiReconEnergyCounts", ";Energy [GeV]", eBins, eMin, eMax },"EMin","Counts");
  auto ReconCounts1 = d2.Filter("Counts1>10&&RMS1<30").Histo1D({"PhiReconEnergyCounts1", ";Energy [GeV]", eBins, eMin, eMax },"EMin","Counts1");
  auto ReconCounts2 = d2.Filter("Counts2>10&&RMS2<30").Histo1D({"PhiReconEnergyCounts2", ";Energy [GeV]", eBins, eMin, eMax },"EMin","Counts2");
  auto ReconCounts3 = d2.Filter("Counts3>10&&RMS3<30").Histo1D({"PhiReconEnergyCounts3", ";Energy [GeV]", eBins, eMin, eMax },"EMin","Counts3");
  auto ReconCounts4 = d2.Filter("Counts4>10&&RMS4<30").Histo1D({"PhiReconEnergyCounts4", ";Energy [GeV]", eBins, eMin, eMax },"EMin","Counts4");
  auto ReconCounts5 = d2.Filter("Counts5>10&&RMS5<30").Histo1D({"PhiReconEnergyCounts5", ";Energy [GeV]", eBins, eMin, eMax },"EMin","Counts5");
  auto ReconCounts6 = d2.Filter("Counts6>10&&RMS6<30").Histo1D({"PhiReconEnergyCounts6", ";Energy [GeV]", eBins, eMin, eMax },"EMin","Counts6");

  auto RCounts      = d2.Histo1D({"TotalEnergyCountsR", ";Energy [GeV]", eBins, eMin, eMax },"EMin","RCounts");
  auto RReconCounts = d2.Filter("Counts>10&&RMS<30")  .Histo1D({"PhiReconEnergyCounts",  ";Energy [GeV]", eBins, eMin, eMax },"EMin","RCounts" );
  auto RCounts1     = d2.Filter("Counts1>10&&RMS1<30").Histo1D({"PhiReconEnergyCounts1", ";Energy [GeV]", eBins, eMin, eMax },"EMin","RCounts1");
  auto RCounts2     = d2.Filter("Counts2>10&&RMS2<30").Histo1D({"PhiReconEnergyCounts2", ";Energy [GeV]", eBins, eMin, eMax },"EMin","RCounts2");
  auto RCounts3     = d2.Filter("Counts3>10&&RMS3<30").Histo1D({"PhiReconEnergyCounts3", ";Energy [GeV]", eBins, eMin, eMax },"EMin","RCounts3");
  auto RCounts4     = d2.Filter("Counts4>10&&RMS4<30").Histo1D({"PhiReconEnergyCounts4", ";Energy [GeV]", eBins, eMin, eMax },"EMin","RCounts4");
  auto RCounts5     = d2.Filter("Counts5>10&&RMS5<30").Histo1D({"PhiReconEnergyCounts5", ";Energy [GeV]", eBins, eMin, eMax },"EMin","RCounts5");
  auto RCounts6     = d2.Filter("Counts6>10&&RMS6<30").Histo1D({"PhiReconEnergyCounts6", ";Energy [GeV]", eBins, eMin, eMax },"EMin","RCounts6");

  auto BCounts      = d2.Histo1D({"TotalEnergyCountsR", ";Energy [GeV]", eBins, eMin, eMax },"EMin","BCounts");
  auto BReconCounts = d2.Filter("Counts>10&&RMS<30")  .Histo1D({"PhiReconEnergyCounts",  ";Energy [GeV]", eBins, eMin, eMax },"EMin","BCounts" );
  auto BCounts1     = d2.Filter("Counts1>10&&RMS1<30").Histo1D({"PhiReconEnergyCounts1", ";Energy [GeV]", eBins, eMin, eMax },"EMin","BCounts1");
  auto BCounts2     = d2.Filter("Counts2>10&&RMS2<30").Histo1D({"PhiReconEnergyCounts2", ";Energy [GeV]", eBins, eMin, eMax },"EMin","BCounts2");
  auto BCounts3     = d2.Filter("Counts3>10&&RMS3<30").Histo1D({"PhiReconEnergyCounts3", ";Energy [GeV]", eBins, eMin, eMax },"EMin","BCounts3");
  auto BCounts4     = d2.Filter("Counts4>10&&RMS4<30").Histo1D({"PhiReconEnergyCounts4", ";Energy [GeV]", eBins, eMin, eMax },"EMin","BCounts4");
  auto BCounts5     = d2.Filter("Counts5>10&&RMS5<30").Histo1D({"PhiReconEnergyCounts5", ";Energy [GeV]", eBins, eMin, eMax },"EMin","BCounts5");
  auto BCounts6     = d2.Filter("Counts6>10&&RMS6<30").Histo1D({"PhiReconEnergyCounts6", ";Energy [GeV]", eBins, eMin, eMax },"EMin","BCounts6");

  auto R2ReconCounts = d2.Filter("Counts>10&&RMS<30&&Ratio>0.5")  .Histo1D({"PhiReconEnergyCounts",  ";Energy [GeV]", eBins, eMin, eMax },"EMin","RCounts" );
  auto R2Counts1     = d2.Filter("Counts1>10&&RMS1<30&&Ratio>0.5").Histo1D({"PhiReconEnergyCounts1", ";Energy [GeV]", eBins, eMin, eMax },"EMin","RCounts1");
  auto R2Counts2     = d2.Filter("Counts2>10&&RMS2<30&&Ratio>0.5").Histo1D({"PhiReconEnergyCounts2", ";Energy [GeV]", eBins, eMin, eMax },"EMin","RCounts2");
  auto R2Counts3     = d2.Filter("Counts3>10&&RMS3<30&&Ratio>0.5").Histo1D({"PhiReconEnergyCounts3", ";Energy [GeV]", eBins, eMin, eMax },"EMin","RCounts3");
  auto R2Counts4     = d2.Filter("Counts4>10&&RMS4<30&&Ratio>0.5").Histo1D({"PhiReconEnergyCounts4", ";Energy [GeV]", eBins, eMin, eMax },"EMin","RCounts4");
  auto R2Counts5     = d2.Filter("Counts5>10&&RMS5<30&&Ratio>0.5").Histo1D({"PhiReconEnergyCounts5", ";Energy [GeV]", eBins, eMin, eMax },"EMin","RCounts5");
  auto R2Counts6     = d2.Filter("Counts6>10&&RMS6<30&&Ratio>0.5").Histo1D({"PhiReconEnergyCounts6", ";Energy [GeV]", eBins, eMin, eMax },"EMin","RCounts6");

//   auto B2Counts      = d2.Histo1D({"TotalEnergyCountsR", ";Energy [GeV]", eBins, eMin, eMax },"EMin","BCounts");
//   auto B2ReconCounts = d2.Filter("Counts>10&&RMS<30&&Ratio>0.1")  .Histo1D({"PhiReconEnergyCounts",  ";Energy [GeV]", eBins, eMin, eMax },"EMin","BCounts" );
//   auto B2Counts1     = d2.Filter("Counts1>10&&RMS1<30&&Ratio>0.1").Histo1D({"PhiReconEnergyCounts1", ";Energy [GeV]", eBins, eMin, eMax },"EMin","BCounts1");
//   auto B2Counts2     = d2.Filter("Counts2>10&&RMS2<30&&Ratio>0.1").Histo1D({"PhiReconEnergyCounts2", ";Energy [GeV]", eBins, eMin, eMax },"EMin","BCounts2");
//   auto B2Counts3     = d2.Filter("Counts3>10&&RMS3<30&&Ratio>0.1").Histo1D({"PhiReconEnergyCounts3", ";Energy [GeV]", eBins, eMin, eMax },"EMin","BCounts3");
//   auto B2Counts4     = d2.Filter("Counts4>10&&RMS4<30&&Ratio>0.1").Histo1D({"PhiReconEnergyCounts4", ";Energy [GeV]", eBins, eMin, eMax },"EMin","BCounts4");
//   auto B2Counts5     = d2.Filter("Counts5>10&&RMS5<30&&Ratio>0.1").Histo1D({"PhiReconEnergyCounts5", ";Energy [GeV]", eBins, eMin, eMax },"EMin","BCounts5");
//   auto B2Counts6     = d2.Filter("Counts6>10&&RMS6<30&&Ratio>0.1").Histo1D({"PhiReconEnergyCounts6", ";Energy [GeV]", eBins, eMin, eMax },"EMin","BCounts6");


  //Proper Efficiency
  TotalCounts->Sumw2(0);
  TotalCounts->Sumw2();
  ReconCounts->Sumw2(0);
  ReconCounts->Sumw2();
  ReconCounts1->Sumw2(0);
  ReconCounts1->Sumw2();
  ReconCounts2->Sumw2(0);
  ReconCounts2->Sumw2();
  ReconCounts3->Sumw2(0);
  ReconCounts3->Sumw2();
  ReconCounts4->Sumw2(0);
  ReconCounts4->Sumw2();
  ReconCounts5->Sumw2(0);
  ReconCounts5->Sumw2();
  ReconCounts6->Sumw2(0);
  ReconCounts6->Sumw2();

  RCounts->Sumw2(0);
  RCounts->Sumw2();
  RReconCounts->Sumw2(0);
  RReconCounts->Sumw2();
  RCounts1->Sumw2(0);
  RCounts1->Sumw2();
  RCounts2->Sumw2(0);
  RCounts2->Sumw2();
  RCounts3->Sumw2(0);
  RCounts3->Sumw2();
  RCounts4->Sumw2(0);
  RCounts4->Sumw2();
  RCounts5->Sumw2(0);
  RCounts5->Sumw2();
  RCounts6->Sumw2(0);
  RCounts6->Sumw2();

  //  double meanQR     = 0.00374;
  double meanQR     = 0.5;
  double meanBrems  = 11.65;
  double PRratio    = meanBrems/meanQR;

  BCounts->Scale(PRratio);
  BReconCounts->Scale(PRratio);
  BCounts1->Scale(PRratio);
  BCounts2->Scale(PRratio);
  BCounts3->Scale(PRratio);
  BCounts4->Scale(PRratio);
  BCounts5->Scale(PRratio);
  BCounts6->Scale(PRratio);
  
  BCounts->Sumw2(0);
  BCounts->Sumw2();
  BReconCounts->Sumw2(0);
  BReconCounts->Sumw2();
  BCounts1->Sumw2(0);
  BCounts1->Sumw2();
  BCounts2->Sumw2(0);
  BCounts2->Sumw2();
  BCounts3->Sumw2(0);
  BCounts3->Sumw2();
  BCounts4->Sumw2(0);
  BCounts4->Sumw2();
  BCounts5->Sumw2(0);
  BCounts5->Sumw2();
  BCounts6->Sumw2(0);
  BCounts6->Sumw2();

  TH1D* efficiency  = (TH1D*)ReconCounts->Clone("Efficiency");
  efficiency->Divide(ReconCounts.GetPtr(),TotalCounts.GetPtr());
  efficiency->SetTitle("Phi Reconstruction Efficiency");
  TH1D* efficiency1 = (TH1D*)ReconCounts1->Clone("Efficiency1");
  efficiency1->Divide(ReconCounts1.GetPtr(),TotalCounts.GetPtr());
  efficiency1->SetTitle("Efficiencies1");
  TH1D* efficiency2 = (TH1D*)ReconCounts2->Clone("Efficiency2");
  efficiency2->Divide(ReconCounts2.GetPtr(),TotalCounts.GetPtr());
  efficiency2->SetTitle("Efficiencies2");
  TH1D* efficiency3 = (TH1D*)ReconCounts3->Clone("Efficiency3");
  efficiency3->Divide(ReconCounts3.GetPtr(),TotalCounts.GetPtr());
  efficiency3->SetTitle("Efficiencies3");
  TH1D* efficiency4 = (TH1D*)ReconCounts4->Clone("Efficiency4");
  efficiency4->Divide(ReconCounts4.GetPtr(),TotalCounts.GetPtr());
  efficiency4->SetTitle("Efficiencies4");
  TH1D* efficiency5 = (TH1D*)ReconCounts5->Clone("Efficiency5");
  efficiency5->Divide(ReconCounts5.GetPtr(),TotalCounts.GetPtr());
  efficiency5->SetTitle("Efficiencies5");
  TH1D* efficiency6 = (TH1D*)ReconCounts6->Clone("Efficiency6");
  efficiency6->Divide(ReconCounts6.GetPtr(),TotalCounts.GetPtr());
  efficiency6->SetTitle("Efficiencies6");

  
  TCanvas* can = new TCanvas("can","can",2200,1000);
  can->Divide(3,1);

  can->cd(1);
  efficiency->SetMinimum(0);
  efficiency->SetMaximum(1);
  efficiency->Draw("hist");
  efficiency1->SetLineColor(kRed);
  efficiency1->Draw("same hist");
  efficiency2->SetLineColor(kBlue);
  efficiency2->Draw("same hist");
  efficiency3->SetLineColor(kGreen);
  efficiency3->Draw("same hist");
  efficiency4->SetLineColor(kYellow);
  efficiency4->Draw("same hist");
  efficiency5->SetLineColor(kOrange);
  efficiency5->Draw("same hist");
  efficiency6->SetLineColor(kGray);
  efficiency6->Draw("same hist");
    
  auto legend = new TLegend(0.2,0.6,0.45,0.88);
  legend->AddEntry(efficiency,"Simulation Accepted","l");
  legend->AddEntry(efficiency6,"16.9 cm height","l");
  legend->AddEntry(efficiency5,"14.1 cm height","l");
  legend->AddEntry(efficiency4,"11.3 cm height","l");
  legend->AddEntry(efficiency3,"8.4 cm height","l");
  legend->AddEntry(efficiency2,"5.6 cm height","l");
  legend->AddEntry(efficiency1,"2.8 cm height","l");
  legend->Draw();

//   auto legend = new TLegend(0.1,0.7,0.48,0.9);
//   legend->AddEntry(efficiency,"simulated","l");
//   legend->AddEntry(efficiency1,"1Timepix4","l");
//   legend->AddEntry(efficiency2,"2Timepix4","l");
//   legend->AddEntry(efficiency3,"3Timepix4","l");
//   legend->AddEntry(efficiency4,"4Timepix4","l");
//   legend->AddEntry(efficiency5,"5Timepix4","l");
//   legend->AddEntry(efficiency6,"6Timepix4","l");
//   legend->Draw();

  TH1D* fracEff1 = (TH1D*)efficiency1->Clone("FracEff1");
  fracEff1->Divide(fracEff1,efficiency);
  TH1D* fracEff2 = (TH1D*)efficiency2->Clone("FracEff2");
  fracEff2->Divide(fracEff2,efficiency);
  TH1D* fracEff3 = (TH1D*)efficiency3->Clone("FracEff3");
  fracEff3->Divide(fracEff3,efficiency);
  TH1D* fracEff4 = (TH1D*)efficiency4->Clone("FracEff4");
  fracEff4->Divide(fracEff4,efficiency);
  TH1D* fracEff5 = (TH1D*)efficiency5->Clone("FracEff5");
  fracEff5->Divide(fracEff5,efficiency);
  TH1D* fracEff6 = (TH1D*)efficiency6->Clone("FracEff6");
  fracEff6->Divide(fracEff6,efficiency);

  can->cd(2);
  fracEff1->SetLineColor(kRed);
  fracEff1->Draw("hist");
  fracEff2->SetLineColor(kBlue);
  fracEff2->Draw("same hist");
  fracEff3->SetLineColor(kGreen);
  fracEff3->Draw("same hist");
  fracEff4->SetLineColor(kYellow);
  fracEff4->Draw("same hist");
  fracEff5->SetLineColor(kOrange);
  fracEff5->Draw("same hist");
  fracEff6->SetLineColor(kGray);
  fracEff6->Draw("same hist");

//   auto legend2 = new TLegend(0.1,0.7,0.48,0.9);
//   legend2->AddEntry(fracEff1,"1Timepix4","l");
//   legend2->AddEntry(fracEff2,"2Timepix4","l");
//   legend2->AddEntry(fracEff3,"3Timepix4","l");
//   legend2->AddEntry(fracEff4,"4Timepix4","l");
//   legend2->AddEntry(fracEff5,"5Timepix4","l");
//   legend2->AddEntry(fracEff6,"6Timepix4","l");
//   legend2->Draw();
  

  
  ROOT::VecOps::RVec<double> xpoints2 = {0,1,2,3,4,5,6};
  ROOT::VecOps::RVec<double> ypoints2 = {ReconCounts->Integral(),ReconCounts1->Integral(),ReconCounts2->Integral(),ReconCounts3->Integral(),ReconCounts4->Integral(),ReconCounts5->Integral(),ReconCounts6->Integral()};
  ypoints2/=TotalCounts->Integral();

  TGraph eff2(xpoints2.size(), &xpoints2[0], &ypoints2[0]);
  can->cd(3);
  eff2.SetMinimum(0);
  eff2.SetMarkerStyle(2);
  eff2.DrawClone("AP");
 
  can->SaveAs(outNamepng);


  TCanvas* can0 = new TCanvas("can0","can0",1000,1000);
  
  efficiency->SetMinimum(0);
  efficiency->SetMaximum(1);
  efficiency->Draw("hist");
//   efficiency1->SetLineColor(kRed);
//   efficiency1->Draw("same hist");
//   efficiency2->SetLineColor(kBlue);
//   efficiency2->Draw("same hist");
//   efficiency3->SetLineColor(kGreen);
//   efficiency3->Draw("same hist");
//   efficiency4->SetLineColor(kYellow);
//   efficiency4->Draw("same hist");
//   efficiency5->SetLineColor(kOrange);
//   efficiency5->Draw("same hist");
//   efficiency6->SetLineColor(kGray);
//   efficiency6->Draw("same hist");
//   legend->Draw();
  can0->SaveAs(outNamepng0);

  TCanvas* can2 = new TCanvas("can2","can2",2200,1000);
  can2->Divide(4,1);


  TH1D* bremRatio  = (TH1D*)ReconCounts->Clone("TotalRatio");
  TH1D* bremSum    = (TH1D*)ReconCounts->Clone("TotalSum");
  bremSum->Add(TotalCounts.GetPtr(),BCounts.GetPtr());
  bremRatio->Divide(TotalCounts.GetPtr(),bremSum);
  bremRatio->SetTitle("TotalRatio");
  TH1D* bremReconRatio  = (TH1D*)ReconCounts->Clone("ReconRatio");
  TH1D* bremReconSum    = (TH1D*)ReconCounts->Clone("ReconSum");
  bremReconSum->Add(RReconCounts.GetPtr(),BReconCounts.GetPtr());
  bremReconRatio->Divide(RReconCounts.GetPtr(),bremReconSum);
  bremReconRatio->SetTitle("ReconRatio");
  TH1D* bremRatio1  = (TH1D*)ReconCounts->Clone("Ratio1");
  TH1D* bremSum1    = (TH1D*)ReconCounts->Clone("Sum1");
  bremSum1->Add(RCounts1.GetPtr(),BCounts1.GetPtr());
  bremRatio1->Divide(RCounts1.GetPtr(),bremSum1);
  bremRatio1->SetTitle("Ratio1");
  TH1D* bremRatio2  = (TH1D*)ReconCounts->Clone("Ratio2");
  TH1D* bremSum2    = (TH1D*)ReconCounts->Clone("Sum2");
  bremSum2->Add(RCounts2.GetPtr(),BCounts2.GetPtr());
  bremRatio2->Divide(RCounts2.GetPtr(),bremSum2);
  bremRatio2->SetTitle("Ratio2");
  TH1D* bremRatio3  = (TH1D*)ReconCounts->Clone("Ratio3");
  TH1D* bremSum3    = (TH1D*)ReconCounts->Clone("Sum3");
  bremSum3->Add(RCounts3.GetPtr(),BCounts3.GetPtr());
  bremRatio3->Divide(RCounts3.GetPtr(),bremSum3);
  bremRatio3->SetTitle("Ratio3");
  TH1D* bremRatio4  = (TH1D*)ReconCounts->Clone("Ratio4");
  TH1D* bremSum4    = (TH1D*)ReconCounts->Clone("Sum4");
  bremSum4->Add(RCounts4.GetPtr(),BCounts4.GetPtr());
  bremRatio4->Divide(RCounts4.GetPtr(),bremSum4);
  bremRatio4->SetTitle("Ratio4");
  TH1D* bremRatio5  = (TH1D*)ReconCounts->Clone("Ratio5");
  TH1D* bremSum5    = (TH1D*)ReconCounts->Clone("Sum5");
  bremSum5->Add(RCounts5.GetPtr(),BCounts5.GetPtr());
  bremRatio5->Divide(RCounts5.GetPtr(),bremSum5);
  bremRatio5->SetTitle("Ratio5");
  TH1D* bremRatio6  = (TH1D*)ReconCounts->Clone("Ratio6");
  TH1D* bremSum6    = (TH1D*)ReconCounts->Clone("Sum6");
  bremSum6->Add(RCounts6.GetPtr(),BCounts6.GetPtr());
  bremRatio6->Divide(RCounts6.GetPtr(),bremSum6);
  bremRatio6->SetTitle("Ratio6");

  can2->cd(1);
//   bremReconRatio->SetMinimum(0);
//   bremReconRatio->SetMaximum(1);
  bremReconRatio->Draw("hist");
  bremRatio1->SetLineColor(kRed);
  bremRatio1->Draw("same hist");
  bremRatio2->SetLineColor(kBlue);
  bremRatio2->Draw("same hist");
  bremRatio3->SetLineColor(kGreen);
  bremRatio3->Draw("same hist");
  bremRatio4->SetLineColor(kYellow);
  bremRatio4->Draw("same hist");
  bremRatio5->SetLineColor(kOrange);
  bremRatio5->Draw("same hist");
  bremRatio6->SetLineColor(kGray);
  bremRatio6->Draw("same hist");
  
  auto legend2 = new TLegend(0.1,0.6,0.48,0.9);
  legend2->AddEntry(bremReconRatio,"simulated","l");
  legend2->AddEntry(bremRatio1,"1Timepix4","l");
  legend2->AddEntry(bremRatio2,"2Timepix4","l");
  legend2->AddEntry(bremRatio3,"3Timepix4","l");
  legend2->AddEntry(bremRatio4,"4Timepix4","l");
  legend2->AddEntry(bremRatio5,"5Timepix4","l");
  legend2->AddEntry(bremRatio6,"6Timepix4","l");
  legend2->Draw();


  //Other stuff
  TH1D* BfracEff1 = (TH1D*)bremRatio1->Clone("BfracEff1");
  BfracEff1->Divide(bremRatio1,bremReconRatio);
  TH1D* BfracEff2 = (TH1D*)bremRatio2->Clone("BfracEff2");
  BfracEff2->Divide(bremRatio2,bremReconRatio);
  TH1D* BfracEff3 = (TH1D*)bremRatio3->Clone("BfracEff3");
  BfracEff3->Divide(bremRatio3,bremReconRatio);
  TH1D* BfracEff4 = (TH1D*)bremRatio4->Clone("BfracEff4");
  BfracEff4->Divide(bremRatio4,bremReconRatio);
  TH1D* BfracEff5 = (TH1D*)bremRatio5->Clone("BfracEff5");
  BfracEff5->Divide(bremRatio5,bremReconRatio);
  TH1D* BfracEff6 = (TH1D*)bremRatio6->Clone("BfracEff6");
  BfracEff6->Divide(bremRatio6,bremReconRatio);

  can2->cd(2);

  BfracEff1->SetLineColor(kRed);
  BfracEff1->Draw("hist");
  BfracEff2->SetLineColor(kBlue);
  BfracEff2->Draw("same hist");
  BfracEff3->SetLineColor(kGreen);
  BfracEff3->Draw("same hist");
  BfracEff4->SetLineColor(kYellow);
  BfracEff4->Draw("same hist");
  BfracEff5->SetLineColor(kOrange);
  BfracEff5->Draw("same hist");
  BfracEff6->SetLineColor(kGray);
  BfracEff6->Draw("same hist");

  can2->cd(3);

  TH1D* brem2Ratio  = (TH1D*)ReconCounts->Clone("TotalRatio");
  brem2Ratio->Divide(R2ReconCounts.GetPtr(),RReconCounts.GetPtr());
  brem2Ratio->SetTitle("TotalRatio");
  TH1D* brem2Ratio1  = (TH1D*)ReconCounts->Clone("Ratio1");
  brem2Ratio1->Divide(R2Counts1.GetPtr(),RCounts1.GetPtr());
  brem2Ratio1->SetTitle("Ratio1");
  TH1D* brem2Ratio2  = (TH1D*)ReconCounts->Clone("Ratio2");
  brem2Ratio2->Divide(R2Counts2.GetPtr(),RCounts2.GetPtr());
  brem2Ratio2->SetTitle("Ratio2");
  TH1D* brem2Ratio3  = (TH1D*)ReconCounts->Clone("Ratio3");
  brem2Ratio3->Divide(R2Counts3.GetPtr(),RCounts3.GetPtr());
  brem2Ratio3->SetTitle("Ratio3");
  TH1D* brem2Ratio4  = (TH1D*)ReconCounts->Clone("Ratio4");
  brem2Ratio4->Divide(R2Counts4.GetPtr(),RCounts4.GetPtr());
  brem2Ratio4->SetTitle("Ratio4");
  TH1D* brem2Ratio5  = (TH1D*)ReconCounts->Clone("Ratio5");
  brem2Ratio5->Divide(R2Counts5.GetPtr(),RCounts5.GetPtr());
  brem2Ratio5->SetTitle("Ratio5");
  TH1D* brem2Ratio6  = (TH1D*)ReconCounts->Clone("Ratio6");
  brem2Ratio6->Divide(R2Counts6.GetPtr(),RCounts6.GetPtr());
  brem2Ratio6->SetTitle("Ratio6");

  brem2Ratio->SetMinimum(0);
  brem2Ratio->SetMaximum(1);
  brem2Ratio->Draw("hist");
  brem2Ratio1->SetLineColor(kRed);
  brem2Ratio1->Draw("same hist");
  brem2Ratio2->SetLineColor(kBlue);
  brem2Ratio2->Draw("same hist");
  brem2Ratio3->SetLineColor(kGreen);
  brem2Ratio3->Draw("same hist");
  brem2Ratio4->SetLineColor(kYellow);
  brem2Ratio4->Draw("same hist");
  brem2Ratio5->SetLineColor(kOrange);
  brem2Ratio5->Draw("same hist");
  brem2Ratio6->SetLineColor(kGray);
  brem2Ratio6->Draw("same hist");


  //Other stuff
  TH1D* B2fracEff1 = (TH1D*)brem2Ratio1->Clone("B2fracEff1");
  B2fracEff1->Divide(brem2Ratio1,brem2Ratio);
  TH1D* B2fracEff2 = (TH1D*)bremRatio2->Clone("B2fracEff2");
  B2fracEff2->Divide(brem2Ratio2,brem2Ratio);
  TH1D* B2fracEff3 = (TH1D*)bremRatio3->Clone("B2fracEff3");
  B2fracEff3->Divide(brem2Ratio3,brem2Ratio);
  TH1D* B2fracEff4 = (TH1D*)bremRatio4->Clone("B2fracEff4");
  B2fracEff4->Divide(brem2Ratio4,brem2Ratio);
  TH1D* B2fracEff5 = (TH1D*)bremRatio5->Clone("B2fracEff5");
  B2fracEff5->Divide(brem2Ratio5,brem2Ratio);
  TH1D* B2fracEff6 = (TH1D*)bremRatio6->Clone("B2fracEff6");
  B2fracEff6->Divide(brem2Ratio6,brem2Ratio);

  can2->cd(4);

  B2fracEff1->SetMinimum(0);
  B2fracEff1->SetMaximum(1);
  B2fracEff1->SetLineColor(kRed);
  B2fracEff1->Draw("hist");
  B2fracEff2->SetLineColor(kBlue);
  B2fracEff2->Draw("same hist");
  B2fracEff3->SetLineColor(kGreen);
  B2fracEff3->Draw("same hist");
  B2fracEff4->SetLineColor(kYellow);
  B2fracEff4->Draw("same hist");
  B2fracEff5->SetLineColor(kOrange);
  B2fracEff5->Draw("same hist");
  B2fracEff6->SetLineColor(kGray);
  B2fracEff6->Draw("same hist");

  

//   can2->cd(3);
//   auto ratio2D  = d2.Filter("Counts>10&&RMS<30") .Histo2D({"BremRatio",   ";#Ratio [GeV]", 100, 0, 1, eBins, eMin, eMax },"Ratio","EMin" );
//   auto ratio  = d2.Filter("Counts>10&&RMS<30") .Histo1D({"BremRatio",   ";#Ratio [GeV]", 40, 0, 1 },"Ratio" );
//   auto ratio1 = d2.Filter("Counts1>10&&RMS1<30") .Histo1D({"BremRatio1",  ";#Ratio [GeV]", 40, 0, 1 },"Ratio1" );
//   auto ratio2 = d2.Filter("Counts2>10&&RMS2<30") .Histo1D({"BremRatio2",  ";#Ratio [GeV]", 40, 0, 1 },"Ratio2" );
//   auto ratio3 = d2.Filter("Counts3>10&&RMS3<30") .Histo1D({"BremRatio3",  ";#Ratio [GeV]", 40, 0, 1 },"Ratio3" );
//   auto ratio4 = d2.Filter("Counts4>10&&RMS4<30") .Histo1D({"BremRatio4",  ";#Ratio [GeV]", 40, 0, 1 },"Ratio4" );
//   auto ratio5 = d2.Filter("Counts5>10&&RMS5<30") .Histo1D({"BremRatio5",  ";#Ratio [GeV]", 40, 0, 1 },"Ratio5" );
//   auto ratio6 = d2.Filter("Counts6>10&&RMS6<30") .Histo1D({"BremRatio6",  ";#Ratio [GeV]", 40, 0, 1 },"Ratio6" );


//   ratio2D->DrawClone("colz");

//   can2->cd(4);

//   ratio->DrawClone("hist");
//   ratio1->SetLineColor(kRed);
//   ratio1->DrawClone("same hist");
//   ratio2->SetLineColor(kBlue);
//   ratio2->DrawClone("same hist");
//   ratio3->SetLineColor(kGreen);
//   ratio3->DrawClone("same hist");
//   ratio4->SetLineColor(kYellow);
//   ratio4->DrawClone("same hist");
//   ratio5->SetLineColor(kOrange);
//   ratio5->DrawClone("same hist");
//   ratio6->SetLineColor(kGray);
//   ratio6->DrawClone("same hist");
  

  can2->SaveAs(outNamepng2);

}
