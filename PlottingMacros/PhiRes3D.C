#include "TString.h"
#include "TCanvas.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include <vector>
#include "TH1.h"
//#include "ROOT.h"


using namespace std;

std::string tag = "Decor-NormTE-4Layers-Bigger-4STEP4";

std::vector<TString> fileNames = {"/scratch/EIC/Results/ML-Out/"+tag+"_real_ETP.root","/home/simon/Analysis/EIC_Analysis/Reg_cell-"+tag+".root"};

struct countRMS{
  ROOT::RDF::RResultPtr<unsigned long long> counts;
  ROOT::RDF::RResultPtr<double>  rms;
};
// struct countRMS{
//   int counts;
//   double rms;
// };

struct GetCounts{
  GetCounts(ROOT::RDF::RNode dframe): df(dframe){}

  countRMS operator()(double emin, double emax,double tmin, double tmax,double pmin, double pmax){
    auto tempframe = df.Filter([&emin,&emax,&tmin,&tmax,&pmin,&pmax](double e,float t,double p){return ((e>emin)&&(e<emax)&&(t>tmin)&&(t<tmax)&&(p>pmin)&&(p<pmax));},{"eE","ThetaGen","Phi"});
    auto count = tempframe.Count();
    auto rms   = tempframe.StdDev("PhiRes");
    countRMS retVal = {count,rms};
    return retVal; 
  };

private:
  ROOT::RDF::RNode df;

};


void PhiRes3D(){

  TString outNamepng  = "plots/PResolutions-"+tag+".png";

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

  ROOT::RDataFrame dfcheet("dataset/TestTree",fileNames[0]);

  //Cheating
  auto dfcheet_2 = dfcheet.Define("CalPhi","atan2(DNN_CPU.sin_phiV_,DNN_CPU.cos_phiV_)*TMath::RadToDeg()")
    .Define("Phi","atan2(sin_phiV_,cos_phiV_)*TMath::RadToDeg()")
    .Define("PhiRes","(Phi-CalPhi)")
    .Define("ThetaGen","thetaE*1000")
    .Define("ThetaRecon","DNN_CPU.thetaE*1000")
    .Define("ThetaRes","(ThetaGen-ThetaRecon)")
    .Define("ePred","DNN_CPU.eE")
    .Define("ERes","100*(eE-DNN_CPU.eE)/eE");

//   ROOT::RDataFrame dfcheet("predict",fileNames[1]);

//   auto dfcheet_2 = dfcheet.Define("CalPhi","pPred")
//     .Define("Phi","phiV")
//     .Define("PhiRes","(Phi-CalPhi)")
//     .Define("ThetaGen","thetaE*1000")
//     .Define("ThetaRecon","tPred*1000")
//     .Define("ThetaRes","(ThetaGen-ThetaRecon)")
//     .Define("ERes","100*(eE-ePred)/eE");

  int    eBins     = 12;
  int    pBins     = 15;
  int    tBins     = 22;
  int    resBins   = 50;

  double eMin = 6;
  double eMax = 18;
  double pMin = -180;
  double pMax = 180;
  double tMin = 0.0;
  double tMax = 0.011*1000;


  double PResRange = 1.5*TMath::RadToDeg();
  double TResRange = 0.002*1000;
  double EResRange = 0.05*100;

  auto Phi1D    = dfcheet_2.Filter("ThetaGen>1").Histo1D({"PRes", ";#phi gen-recon [rad]", resBins, -PResRange, PResRange }, "PhiRes");
  auto Phi2D    = dfcheet_2.Filter("ThetaGen>1").Histo2D({"Phi",    "Phi Reconstruction (Theta > 1 mrad);#phi gen [rad];#phi recon [rad]",      pBins,   pMin,   pMax,   pBins,   pMin,   pMax,  }, "Phi",   "CalPhi");
  auto Phi2D2   = dfcheet_2.Filter("ThetaGen>1").Histo2D({"PRes", ";#phi gen-recon [rad];Electron energy [GeV] ",   resBins, -PResRange, PResRange, eBins, eMin, eMax },"PhiRes",  "eE");

  const std::vector< int > bins = {eBins,tBins,pBins,resBins};
  const std::vector< double > min  = {eMin,tMin,pMin,-PResRange};
  const std::vector< double > max  = {eMax,tMax,pMax, PResRange};
 
  auto Phi4D    = dfcheet_2.HistoND({"PRes", ";Electron energy [GeV] ;#theta [mrad]",4, bins,min,max }, {"eE","ThetaGen","Phi","PhiRes"});

  ROOT::VecOps::RVec<double> thetaBins;
  ROOT::VecOps::RVec<double> phiBins;
  ROOT::VecOps::RVec<double> energyBins;
  ROOT::VecOps::RVec<TString> energyFilters;
  ROOT::VecOps::RVec<TString> thetaFilters;
  ROOT::VecOps::RVec<TString> phiFilters;
  double tSize = (tMax-tMin)/tBins;
  double pSize = (pMax-pMin)/pBins;
  double eSize = (eMax-eMin)/eBins;

//   for(int i=0; i<eBins; i++){
//     energyBins.push_back(eMin+(energySize*(i+0.5)));
//     TString energyFilter;
//     energyFilter.Form("eE>%f&&eE<%f",eMin+energySize*(i),eMin+energySize*(i+1));
//     energyFilters.push_back(energyFilter);
//   }
//   for(int i=0; i<tBins; i++){
//     thetaBins.push_back(thetaMin+(thetaSize*(i+0.5)));
//     TString thetaFilter;
//     thetaFilter.Form("ThetaGen>%f&&ThetaGen<%f",thetaMin+thetaSize*(i),thetaMin+thetaSize*(i+1));
//     thetaFilters.push_back(thetaFilter);
//   }
//   for(int i=0; i<pBins; i++){
//     phiBins.push_back(phiMin+(energySize*(i+0.5)));
//     TString phiFilter;
//     phiFilter.Form("Phi>%f&&Phi<%f",phiMin+phiSize*(i),phiMin+phiSize*(i+1));
//     phiFilters.push_back(phiFilter);
//   }

  int totalRows = eBins*tBins*pBins;

  GetCounts getCounts(dfcheet_2);

  ROOT::RDataFrame d(totalRows);

  auto d2 = d.Define("EMin",[&eMin,&eSize,&eBins](ULong64_t i){return eMin+(i%eBins)*eSize;},{"rdfentry_"})
    .Define("EMax",[&eSize](double i){return i+eSize;},{"EMin"})
    .Define("TMin",[&tMin,&tSize,&tBins,&eBins](ULong64_t i){return tMin+((int)(i/eBins)%tBins)*tSize;},{"rdfentry_"})
    .Define("TMax",[&tSize](double i){return i+tSize;},{"TMin"})
    .Define("PMin",[&pMin,&pSize,&pBins,&tBins,&eBins](ULong64_t i){return pMin+((int)(i/(eBins*tBins))%pBins)*pSize;},{"rdfentry_"})
    .Define("PMax",[&pSize](double i){return i+pSize;},{"PMin"})
    .Define("CountsRMS",GetCounts(dfcheet_2),{"EMin","EMax","TMin","TMax","PMin","PMax"})
    .Define("rCounts","CountsRMS.counts")
    .Define("rRMS","CountsRMS.rms")
    .Define("Counts","rCounts.GetValue()")
    .Define("RMS","rRMS.GetValue()");

  d2.Snapshot("rms","RMS.root",{"EMin","EMax","TMin","TMax","PMin","PMax","Counts","RMS"});


//   Root::VecOps::RVec<ROOT::RDF::RResultPtr<unsigned long long>> counts;
//   Root::VecOps::RVec<ROOT::RDF::RResultPtr<double>> rmss;

//   for(auto eFilter:energyFilters){
//     auto df2 = dfcheet_2.Filter(eFilter.Data());
//     for(auto tFilter:thetaFilters){
//       auto df3 = df2.Filter(tFilter.Data());
//       for(auto pFilter:phiFilters){
// 	auto df4 = df3.Filter(pFilter.Data());
// 	auto count = df4.Count();
// 	counts.push_back(count);
// 	auto rms = df4.StdDev("PhiRes");
// 	rmss.push_back(rms);
//       }
//     }
//   }

//   for(int i=0; i<rmss.size(); i++){
//     std::cout << rmss[i].GetValue() << " " << counts[i].GetValue()<< std::endl;
//   }

//   for(int i=0; i<=tBins; i++){
//     thetaBins.push_back(thetaMin+(thetaSize*i));
//   }
//   for(int i=0; i<=pBins; i++){
//     phiBins.push_back(phiMin+(phiSize*i));
//   }

//   ROOT::VecOps::RVec<double> Theta;
//   ROOT::VecOps::RVec<double> Phi;
//   ROOT::VecOps::RVec<double> E;
//   ROOT::VecOps::RVec<double> N;
//   ROOT::VecOps::RVec<double> StdDev;

//   for(int i=0; i<eBins; i++){
//   }
  

 //  (Phi3D.GetPtr())->FitSlicesZ("L");

//   TH2D* phiSigma    = (TH2D*)(gDirectory->Get("PRes_2"));



//   can->cd(3);  
//   phiSigma->SetTitle("Phi Resolutions [rad]");
//   phiSigma->SetStats(0);
//   phiSigma->SetMaximum(0.3);
//   phiSigma->Draw("colz");
   
  can->cd(6);
  gPad->SetLogz();
  Phi2D2.GetPtr()->Draw("colz");
  Phi2D2->SetStats(0);
  
  can->SaveAs(outNamepng);


}
