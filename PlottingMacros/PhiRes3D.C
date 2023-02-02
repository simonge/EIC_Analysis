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

std::vector<TString> fileNames = {"/scratch/EIC/Results/ML-Out/"+tag+"_real_ETP.root","/home/simon/Analysis/EIC_Analysis/NeuralNetwork/Reg_point-"+tag+".root"};

// struct countRMS{
//   ROOT::RDF::RResultPtr<unsigned long long> counts;
//   ROOT::RDF::RResultPtr<double>  rms;
// };
struct countRMS{
  double counts;
  double rms;
};

// struct GetCounts{
//   GetCounts(ROOT::RDF::RNode dframe): df(dframe){}

//   countRMS operator()(double emin, double emax,double tmin, double tmax,double pmin, double pmax){
//     TString filter;
//     filter.Form("eE>%f&&eE<%f&&ThetaGen>%f&&ThetaGen<%f&&Phi>%f&&Phi<%f",emin,emax,tmin,tmax,pmin,pmax);
//     auto tempframe = df.Filter(filter.Data());
//     //    auto tempframe = df.Filter([&emin,&emax,&tmin,&tmax,&pmin,&pmax](double e,float t,double p){return ((e>emin)&&(e<emax)&&(t>tmin)&&(t<tmax)&&(p>pmin)&&(p<pmax));},{"eE","ThetaGen","Phi"});
//     auto count = tempframe.Count();
//     auto rms   = tempframe.StdDev("PhiRes");
//     countRMS retVal = {count,rms};
//     return retVal; 
//   };

// private:
//   ROOT::RDF::RNode df;

// };

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
      TFitResultPtr r = HistProj->Fit("gaus","LS");
      rms   = r->Parameter(2);
      //auto rms   = HistProj->GetStdDev();
    }
    
    countRMS retVal = {count,rms};
    return retVal; 
  };

private:
  THnD* Hist;

};


void PhiRes3D(){
  //ROOT::EnableImplicitMT();

  TString outNamepng  = "plots/PEfficiencies3-"+tag+".png";

  gStyle->SetStatW(0.3);
  gStyle->SetStatColor(0);
  gStyle->SetStatBorderSize(1);
  //gStyle->SetOptStat(0);
  gStyle->SetOptStat(1111);
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

  TCanvas* can = new TCanvas("can","can",2200,1000);
  can->Divide(4,2);

//   ROOT::RDataFrame dfcheet("dataset/TestTree",fileNames[0]);

//   //Cheating
//   auto dfcheet_2 = dfcheet.Define("CalPhi","atan2(DNN_CPU.sin_phiV_,DNN_CPU.cos_phiV_)*TMath::RadToDeg()")
//     .Define("Phi","atan2(sin_phiV_,cos_phiV_)*TMath::RadToDeg()")
//     .Define("PhiRes","(Phi-CalPhi)")
//     .Define("ThetaGen","thetaE*1000")
//     .Define("ThetaRecon","DNN_CPU.thetaE*1000")
//     .Define("ThetaRes","(ThetaGen-ThetaRecon)")
//     .Define("ePred","DNN_CPU.eE")
//     .Define("ERes","100*(eE-DNN_CPU.eE)/eE");

  ROOT::RDataFrame dfcheet("predict",fileNames[1]);

  auto dfcheet_2 = dfcheet//.Filter("fit_minchi2<0.002")
    .Define("Phi","phiV*TMath::RadToDeg()")
    .Define("ThetaGen","thetaE*1000")
    .Define("CalPhiR","pPredR*TMath::RadToDeg()")
    .Define("PhiResR","(Phi-CalPhiR)")
    .Define("ThetaReconR","tPredR*1000")
    .Define("ThetaResR","(ThetaGen-ThetaReconR)")
    .Define("EResR","100*(eE-ePredR)/eE")
    .Define("CalPhiC","pPredC*TMath::RadToDeg()")
    .Define("PhiResC","(Phi-CalPhiC)")
    .Define("ThetaReconC","tPredC*1000")
    .Define("ThetaResC","(ThetaGen-ThetaReconC)")
    .Define("EResC","100*(eE-ePredC)/eE")
    .Define("CalPhiF","pPredF*TMath::RadToDeg()")
    .Define("PhiResF","(Phi-CalPhiF)")
    .Define("ThetaReconF","tPredF*1000")
    .Define("ThetaResF","(ThetaGen-ThetaReconF)")
    .Define("EResF","100*(eE-ePredF)/eE");

  int    eBins     = 40;
  int    pBins     = 20;
  int    tBins     = 40;
//   int    eBins     = 24;
//   int    pBins     = 30;
//   int    tBins     = 44;
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
 
  auto Phi4DR    = dfcheet_2.HistoND({"PResR", ";Electron energy [GeV] ;#theta [mrad]",4, bins,min,max }, {"eE","ThetaGen","Phi","PhiResR"});
  auto Phi4DC    = dfcheet_2.HistoND({"PResC", ";Electron energy [GeV] ;#theta [mrad]",4, bins,min,max }, {"eE","ThetaGen","Phi","PhiResC"});
  auto Phi4DF    = dfcheet_2.HistoND({"PResF", ";Electron energy [GeV] ;#theta [mrad]",4, bins,min,max }, {"eE","ThetaGen","Phi","PhiResF"});
  
  int totalRows = eBins*tBins*pBins;
  ROOT::RDataFrame d(totalRows);

  int slots = d.GetNSlots();

  auto d2 = d.Define("EBin",[&eBins](ULong64_t i){return i%eBins;},{"rdfentry_"})
    .Define("EMin",[&eMin,&eBins,&eSize](ULong64_t i){return eMin+i*eSize+(eSize/2);},{"EBin"})
    .Define("TBin",[&tBins,&eBins](ULong64_t i){return ((ULong64_t)(i/eBins)%tBins);},{"rdfentry_"})
    .Define("TMin",[&tMin,&tSize](ULong64_t i){return tMin+i*tSize+(tSize/2);},{"TBin"})
    .Define("PBin",[&pBins,&tBins,&eBins](ULong64_t i){return ((ULong64_t)(i/(eBins*tBins))%pBins);},{"rdfentry_"})
    .Define("PMin",[&pMin,&pSize](ULong64_t i){return pMin+i*pSize+(pSize/2);},{"PBin"})
    .Define("CountsRMSR",GetCounts2(Phi4DR.GetPtr()),{"EBin","TBin","PBin"})
    .Define("Counts","CountsRMSR.counts")
    .Define("RMSR","CountsRMSR.rms")
    .Define("CountsRMSC",GetCounts2(Phi4DC.GetPtr()),{"EBin","TBin","PBin"})
    .Define("RMSC","CountsRMSC.rms")
    .Define("CountsRMSF",GetCounts2(Phi4DF.GetPtr()),{"EBin","TBin","PBin"})
    .Define("RMSF","CountsRMSF.rms");

  auto TotalCounts  = d2.Histo1D({"TotalEnergyCountsR", ";#Energy [GeV]", eBins, eMin, eMax },"EMin","Counts");
  auto ReconCountsR = d2.Filter("Counts>5&&RMSR<20").Histo1D({"PhiReconEnergyCountsR", ";#Energy [GeV]", eBins, eMin, eMax },"EMin","Counts");
  auto ReconCountsC = d2.Filter("Counts>5&&RMSC<20").Histo1D({"PhiReconEnergyCountsC", ";#Energy [GeV]", eBins, eMin, eMax },"EMin","Counts");
  auto ReconCountsF = d2.Filter("Counts>5&&RMSF<20").Histo1D({"PhiReconEnergyCountsF", ";#Energy [GeV]", eBins, eMin, eMax },"EMin","Counts");

  //Proper Efficiency
  TotalCounts->Sumw2(0);
  TotalCounts->Sumw2();
  ReconCountsR->Sumw2(0);
  ReconCountsR->Sumw2();
  ReconCountsC->Sumw2(0);
  ReconCountsC->Sumw2();
  ReconCountsF->Sumw2(0);
  ReconCountsF->Sumw2();

  TH1D* efficiencyR = (TH1D*)ReconCountsR->Clone("EfficiencyR");
  efficiencyR->Divide(ReconCountsR.GetPtr(),TotalCounts.GetPtr());
  efficiencyR->SetTitle("EfficienciesR");
  TH1D* efficiencyC = (TH1D*)ReconCountsC->Clone("EfficiencyC");
  efficiencyC->Divide(ReconCountsC.GetPtr(),TotalCounts.GetPtr());
  efficiencyC->SetTitle("EfficienciesC");
  TH1D* efficiencyF = (TH1D*)ReconCountsF->Clone("EfficiencyF");
  efficiencyF->Divide(ReconCountsF.GetPtr(),TotalCounts.GetPtr());
  efficiencyF->SetTitle("EfficienciesF");

  
  auto TotalCounts2D = d2.Histo2D({"TotalEnergyCounts", ";#Energy [GeV]", eBins, eMin, eMax, tBins, tMin, tMax },"EMin","TMin","Counts");
  auto ReconCounts2DR = d2.Filter("Counts>5&&RMSR<20").Histo2D({"PhiReconEnergyCountsR", ";#Energy [GeV]", eBins, eMin, eMax, tBins, tMin, tMax },"EMin","TMin","Counts");
  auto ReconCounts2DC = d2.Filter("Counts>5&&RMSC<20").Histo2D({"PhiReconEnergyCountsC", ";#Energy [GeV]", eBins, eMin, eMax, tBins, tMin, tMax },"EMin","TMin","Counts");
  auto ReconCounts2DF = d2.Filter("Counts>5&&RMSF<20").Histo2D({"PhiReconEnergyCountsF", ";#Energy [GeV]", eBins, eMin, eMax, tBins, tMin, tMax },"EMin","TMin","Counts");
   
  TH1D* efficiency3R = (TH1D*)ReconCounts2DR->Clone("EfficiencyR");
  efficiency3R->Divide(ReconCounts2DR.GetPtr(),TotalCounts2D.GetPtr());
  efficiency3R->SetTitle("EfficienciesR");
  TH1D* efficiency3C = (TH1D*)ReconCounts2DC->Clone("EfficiencyC");
  efficiency3C->Divide(ReconCounts2DC.GetPtr(),TotalCounts2D.GetPtr());
  efficiency3C->SetTitle("EfficienciesC");
  TH1D* efficiency3F = (TH1D*)ReconCounts2DF->Clone("EfficiencyF");
  efficiency3F->Divide(ReconCounts2DF.GetPtr(),TotalCounts2D.GetPtr());
  efficiency3F->SetTitle("EfficienciesF");

  d2.Snapshot("rms","RMS2.root",{"EBin","TBin","PBin","EMin","TMin","PMin","Counts","RMSR","RMSC","RMSF"});
   

  can->cd(1);
  efficiencyR->SetMinimum(0);
  efficiencyR->SetMaximum(1);
  efficiencyR->Draw("hist");
  efficiencyC->SetLineColor(kRed);
  efficiencyC->Draw("same hist");
  efficiencyF->SetLineColor(kBlue);
  efficiencyF->Draw("same hist");

  can->cd(2);
  gPad->SetLogz();
  ReconCounts2DR->DrawClone("colz");
  can->cd(3);
  gPad->SetLogz();
  ReconCounts2DC->DrawClone("colz");
  can->cd(4);
  gPad->SetLogz();
  ReconCounts2DF->DrawClone("colz");

  can->cd(5);
  gPad->SetLogz();
  TotalCounts2D->DrawClone("colz");

  can->cd(6);
  efficiency3R->Draw("colz");
  can->cd(7);
  efficiency3C->Draw("colz");
  can->cd(8);
  efficiency3F->Draw("colz");

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




//   ROOT::RDataFrame d(totalRows);

//   int slots = d.GetNSlots();
//   auto d2 = d.Define("EMin",[&eMin,&eSize,&eBins](ULong64_t i){return eMin+(i%eBins)*eSize;},{"rdfentry_"})
//     .Define("EMax",[&eSize](double i){return i+eSize;},{"EMin"})
//     .Define("TMin",[&tMin,&tSize,&tBins,&eBins](ULong64_t i){return tMin+((int)(i/eBins)%tBins)*tSize;},{"rdfentry_"})
//     .Define("TMax",[&tSize](double i){return i+tSize;},{"TMin"})
//     .Define("PMin",[&pMin,&pSize,&pBins,&tBins,&eBins](ULong64_t i){return pMin+((int)(i/(eBins*tBins))%pBins)*pSize;},{"rdfentry_"})
//     .Define("PMax",[&pSize](double i){return i+pSize;},{"PMin"})
//     .Define("CountsRMS",GetCounts(dfcheet_2),{"EMin","EMax","TMin","TMax","PMin","PMax"})
//     .Define("rCounts","CountsRMS.counts")
//     .Define("rRMS","CountsRMS.rms")
//     .Define("Counts","rCounts.GetValue()")
//     .Define("RMS","rRMS.GetValue()");

//   d2.Snapshot("rms","RMS.root",{"EMin","EMax","TMin","TMax","PMin","PMax","Counts","RMS"});


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
   
  
  can->SaveAs(outNamepng);


}
