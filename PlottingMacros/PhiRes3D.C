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
    auto rms   = HistProj->GetStdDev();
    countRMS retVal = {count,rms};
    return retVal; 
  };

private:
  THnD* Hist;

};


void PhiRes3D(){
  //ROOT::EnableImplicitMT();

  TString outNamepng  = "plots/PResolutions-"+tag+".png";

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

  auto Energy       = dfcheet_2.Histo1D({"TotalEnergy", ";#phi gen-recon [rad]", eBins, eMin, eMax }, "eE");
  auto FilterEnergy = dfcheet_2.Filter("ThetaGen>1").Histo1D({"FilterEnergy", ";#phi gen-recon [rad]", eBins, eMin, eMax }, "eE");

  auto Phi1D    = dfcheet_2.Filter("ThetaGen>1").Histo1D({"PRes", ";#phi gen-recon [rad]", resBins, -PResRange, PResRange }, "PhiRes");
  auto Phi2D    = dfcheet_2.Filter("ThetaGen>1").Histo2D({"Phi",    "Phi Reconstruction (Theta > 1 mrad);#phi gen [rad];#phi recon [rad]",      pBins,   pMin,   pMax,   pBins,   pMin,   pMax,  }, "Phi",   "CalPhi");
  auto Phi2D2   = dfcheet_2.Filter("ThetaGen>1").Histo2D({"PRes", ";#phi gen-recon [rad];Electron energy [GeV] ",   resBins, -PResRange, PResRange, eBins, eMin, eMax },"PhiRes",  "eE");

  const std::vector< int >    bins = {eBins,tBins,pBins,resBins};
  const std::vector< double > min  = {eMin,tMin,pMin,-PResRange};
  const std::vector< double > max  = {eMax,tMax,pMax, PResRange};
  double tSize = (tMax-tMin)/tBins;
  double pSize = (pMax-pMin)/pBins;
  double eSize = (eMax-eMin)/eBins;
 
  auto Phi4D    = dfcheet_2.HistoND({"PRes", ";Electron energy [GeV] ;#theta [mrad]",4, bins,min,max }, {"eE","ThetaGen","Phi","PhiRes"});
  
  int totalRows = eBins*tBins*pBins;
  ROOT::RDataFrame d(totalRows);

  int slots = d.GetNSlots();

  auto d2 = d.Define("EBin",[&eBins](ULong64_t i){return i%eBins;},{"rdfentry_"})
    .Define("EMin",[&eMin,&eBins,&eSize](ULong64_t i){return eMin+i*eSize+(eSize/2);},{"EBin"})
    .Define("TBin",[&tBins,&eBins](ULong64_t i){return ((ULong64_t)(i/eBins)%tBins);},{"rdfentry_"})
    .Define("TMin",[&tMin,&tSize](ULong64_t i){return tMin+i*tSize+(tSize/2);},{"TBin"})
    .Define("PBin",[&pBins,&tBins,&eBins](ULong64_t i){return ((ULong64_t)(i/(eBins*tBins))%pBins);},{"rdfentry_"})
    .Define("PMin",[&pMin,&pSize](ULong64_t i){return pMin+i*pSize+(pSize/2);},{"PBin"})
    .Define("CountsRMS",GetCounts2(Phi4D.GetPtr()),{"EBin","TBin","PBin"})
    .Define("Counts","CountsRMS.counts")
    .Define("RMS","CountsRMS.rms");

  auto TotalCounts = d2.Histo1D({"TotalEnergyCounts", ";#Energy [GeV]", eBins, eMin, eMax },"EMin","Counts");
  auto ReconCounts = d2.Filter("Counts>5&&RMS<30").Histo1D({"PhiReconEnergyCounts", ";#Energy [GeV]", eBins, eMin, eMax },"EMin","Counts");

  //Proper Efficiency
  TotalCounts->Sumw2(0);
  ReconCounts->Sumw2(0);
  TotalCounts->Sumw2();
  ReconCounts->Sumw2();

  TH1D* efficiency = (TH1D*)ReconCounts->Clone("Efficiency");
  efficiency->Divide(ReconCounts.GetPtr(),TotalCounts.GetPtr());
  efficiency->SetTitle("Efficiencies");

  can->cd(1);
  efficiency->SetMinimum(0);
  efficiency->SetMaximum(1);
  efficiency->Draw();

  //ThetaCut Efficiency
  FilterEnergy->Sumw2(0);
  Energy->Sumw2(0);
  FilterEnergy->Sumw2();
  Energy->Sumw2();

  TH1D* efficiency2 = (TH1D*)FilterEnergy->Clone("Efficiency2");
  efficiency2->Divide(FilterEnergy.GetPtr(),Energy.GetPtr());
  efficiency2->SetTitle("Efficiencies2");

  can->cd(2);
  efficiency2->SetMinimum(0);
  efficiency2->SetMaximum(1);
  efficiency2->Draw();

  can->cd(3);
  gPad->SetLogz();
  Phi2D2.GetPtr()->DrawClone("colz");
  //  Phi2D2->SetStats(0);
  
  auto TotalCounts2D = d2.Histo2D({"TotalEnergyCounts", ";#Energy [GeV]", eBins, eMin, eMax, tBins, tMin, tMax },"EMin","TMin","Counts");
  auto ReconCounts2D = d2.Filter("Counts>5&&RMS<30").Histo2D({"PhiReconEnergyCounts", ";#Energy [GeV]", eBins, eMin, eMax, tBins, tMin, tMax },"EMin","TMin","Counts");
   
  TH1D* efficiency3 = (TH1D*)ReconCounts2D->Clone("Efficiency");
  efficiency3->Divide(ReconCounts2D.GetPtr(),TotalCounts2D.GetPtr());
  efficiency3->SetTitle("Efficiencies");

  can->cd(4);
  gPad->SetLogz();
  ReconCounts2D->DrawClone("colz");

  can->cd(5);
  gPad->SetLogz();
  TotalCounts2D->DrawClone("colz");

  can->cd(6);
  efficiency3->Draw("colz");

  d2.Snapshot("rms","RMS.root",{"EBin","TBin","PBin","EMin","TMin","PMin","Counts","RMS"});
  
  

  ROOT::VecOps::RVec<double> thetaBins;
  ROOT::VecOps::RVec<double> phiBins;
  ROOT::VecOps::RVec<double> energyBins;
  ROOT::VecOps::RVec<TString> energyFilters;
  ROOT::VecOps::RVec<TString> thetaFilters;
  ROOT::VecOps::RVec<TString> phiFilters;

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
