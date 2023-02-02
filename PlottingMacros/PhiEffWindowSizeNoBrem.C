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
std::string tag = "FitTrain_uniform";

std::vector<TString> fileNames = {"/scratch/EIC/Results/ML-Out/"+tag+"_ETP.root","/home/simon/Analysis/EIC_Analysis/NeuralNetwork/Reg_point-"+tag+".root"};

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

  TString outNamepng  = "plots/PSizeEfficiencies-"+tag+".png";

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

  int    eBins     = 40;
  int    pBins     = 20;
  int    tBins     = 40;
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
    .Define("RMS6","CountsRMS6.rms");

//   TString countFilter="Counts>10Counts>10";

//   auto TotalCounts  = d2.Histo1D({"TotalEnergyCountsR", ";#Energy [GeV]", eBins, eMin, eMax },"EMin","Counts");
//   auto ReconCounts  = d2.Filter(countFilter+"&&RMS<20").Histo1D({"PhiReconEnergyCounts", ";#Energy [GeV]", eBins, eMin, eMax },"EMin","Counts");
//   auto ReconCounts1 = d2.Filter(countFilter+"&&RMS1<20").Histo1D({"PhiReconEnergyCounts1", ";#Energy [GeV]", eBins, eMin, eMax },"EMin","Counts");
//   auto ReconCounts2 = d2.Filter(countFilter+"&&RMS2<20").Histo1D({"PhiReconEnergyCounts2", ";#Energy [GeV]", eBins, eMin, eMax },"EMin","Counts");
//   auto ReconCounts3 = d2.Filter(countFilter+"&&RMS3<20").Histo1D({"PhiReconEnergyCounts3", ";#Energy [GeV]", eBins, eMin, eMax },"EMin","Counts");
//   auto ReconCounts4 = d2.Filter(countFilter+"&&RMS4<20").Histo1D({"PhiReconEnergyCounts4", ";#Energy [GeV]", eBins, eMin, eMax },"EMin","Counts");
//   auto ReconCounts5 = d2.Filter(countFilter+"&&RMS5<20").Histo1D({"PhiReconEnergyCounts5", ";#Energy [GeV]", eBins, eMin, eMax },"EMin","Counts");
//   auto ReconCounts6 = d2.Filter(countFilter+"&&RMS6<20").Histo1D({"PhiReconEnergyCounts6", ";#Energy [GeV]", eBins, eMin, eMax },"EMin","Counts");

  auto TotalCounts  = d2.Histo1D({"TotalEnergyCountsR", ";#Energy [GeV]", eBins, eMin, eMax },"EMin","Counts");
  auto ReconCounts  = d2.Filter("Counts>10&&RMS<30").Histo1D({"PhiReconEnergyCounts", ";#Energy [GeV]", eBins, eMin, eMax },"EMin","Counts");
  auto ReconCounts1 = d2.Filter("Counts1>10&&RMS1<30").Histo1D({"PhiReconEnergyCounts1", ";#Energy [GeV]", eBins, eMin, eMax },"EMin","Counts1");
  auto ReconCounts2 = d2.Filter("Counts2>10&&RMS2<30").Histo1D({"PhiReconEnergyCounts2", ";#Energy [GeV]", eBins, eMin, eMax },"EMin","Counts2");
  auto ReconCounts3 = d2.Filter("Counts3>10&&RMS3<30").Histo1D({"PhiReconEnergyCounts3", ";#Energy [GeV]", eBins, eMin, eMax },"EMin","Counts3");
  auto ReconCounts4 = d2.Filter("Counts4>10&&RMS4<30").Histo1D({"PhiReconEnergyCounts4", ";#Energy [GeV]", eBins, eMin, eMax },"EMin","Counts4");
  auto ReconCounts5 = d2.Filter("Counts5>10&&RMS5<30").Histo1D({"PhiReconEnergyCounts5", ";#Energy [GeV]", eBins, eMin, eMax },"EMin","Counts5");
  auto ReconCounts6 = d2.Filter("Counts6>10&&RMS6<30").Histo1D({"PhiReconEnergyCounts6", ";#Energy [GeV]", eBins, eMin, eMax },"EMin","Counts6");

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

  TH1D* efficiency  = (TH1D*)ReconCounts->Clone("Efficiency");
  efficiency->Divide(ReconCounts.GetPtr(),TotalCounts.GetPtr());
  efficiency->SetTitle("Efficiencies");
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
  
  auto legend = new TLegend(0.1,0.7,0.48,0.9);
  legend->AddEntry(efficiency,"simulated","l");
  legend->AddEntry(efficiency1,"1Timepix4","l");
  legend->AddEntry(efficiency2,"2Timepix4","l");
  legend->AddEntry(efficiency3,"3Timepix4","l");
  legend->AddEntry(efficiency4,"4Timepix4","l");
  legend->AddEntry(efficiency5,"5Timepix4","l");
  legend->AddEntry(efficiency6,"6Timepix4","l");
  legend->Draw();

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

//   ROOT::RDataFrame dfBrem("dataset/TestTree",fileNames[0]);

//   //Cheating
//   auto dfBrem_2 = dfBrem
//     .Define("Phi","phiV*TMath::RadToDeg()")
//     .Define("ThetaGen","thetaE*1000")
//     .Define("CalPhiC","pPredC*TMath::RadToDeg()")
//     .Define("PhiResC","(Phi-CalPhiC)")
//     .Define("ThetaReconC","tPredC*1000")
//     .Define("ThetaResC","(ThetaGen-ThetaReconC)")
//     .Define("EResC","100*(eE-ePredC)/eE");

//   auto BremPhi4D    = dfBrem_2.HistoND({"BremPResC", ";Electron energy [GeV] ;#theta [mrad]",4, bins,min,max }, {"eE","ThetaGen","Phi","PhiResC"});
//   auto BremPhi4D1   = dfBrem_2.Filter("(ymin>(-256*1))&&(ymax<(256*1))").HistoND({"BremPResC1", ";Electron energy [GeV] ;#theta [mrad]",4, bins,min,max }, {"eE","ThetaGen","Phi","PhiResC"});
//   auto BremPhi4D2   = dfBrem_2.Filter("(ymin>(-256*2))&&(ymax<(256*2))").HistoND({"BremPResC2", ";Electron energy [GeV] ;#theta [mrad]",4, bins,min,max }, {"eE","ThetaGen","Phi","PhiResC"});
//   auto BremPhi4D3   = dfBrem_2.Filter("(ymin>(-256*3))&&(ymax<(256*3))").HistoND({"BremPResC3", ";Electron energy [GeV] ;#theta [mrad]",4, bins,min,max }, {"eE","ThetaGen","Phi","PhiResC"});
//   auto BremPhi4D4   = dfBrem_2.Filter("(ymin>(-256*4))&&(ymax<(256*4))").HistoND({"BremPResC4", ";Electron energy [GeV] ;#theta [mrad]",4, bins,min,max }, {"eE","ThetaGen","Phi","PhiResC"});
//   auto BremPhi4D5   = dfBrem_2.Filter("(ymin>(-256*5))&&(ymax<(256*5))").HistoND({"BremPResC5", ";Electron energy [GeV] ;#theta [mrad]",4, bins,min,max }, {"eE","ThetaGen","Phi","PhiResC"});
//   auto BremPhi4D6   = dfBrem_2.Filter("(ymin>(-256*6))&&(ymax<(256*6))").HistoND({"BremPResC6", ";Electron energy [GeV] ;#theta [mrad]",4, bins,min,max }, {"eE","ThetaGen","Phi","PhiResC"});

//   ROOT::RDataFrame dBrem(totalRows);

//   int slots = dBrem.GetNSlots();

//   auto dBrem2 = dBrem.Define("EBin",[&eBins](ULong64_t i){return i%eBins;},{"rdfentry_"})
//     .Define("EMin",[&eMin,&eBins,&eSize](ULong64_t i){return eMin+i*eSize+(eSize/2);},{"EBin"})
//     .Define("TBin",[&tBins,&eBins](ULong64_t i){return ((ULong64_t)(i/eBins)%tBins);},{"rdfentry_"})
//     .Define("TMin",[&tMin,&tSize](ULong64_t i){return tMin+i*tSize+(tSize/2);},{"TBin"})
//     .Define("PBin",[&pBins,&tBins,&eBins](ULong64_t i){return ((ULong64_t)(i/(eBins*tBins))%pBins);},{"rdfentry_"})
//     .Define("PMin",[&pMin,&pSize](ULong64_t i){return pMin+i*pSize+(pSize/2);},{"PBin"})
//     .Define("CountsRMS",GetCounts2(BremPhi4D.GetPtr()),{"EBin","TBin","PBin"})
//     .Define("Counts","CountsRMS.counts")
//     .Define("RMS","CountsRMS.rms")
//     .Define("CountsRMS1",GetCounts2(BremPhi4D1.GetPtr()),{"EBin","TBin","PBin"})
//     .Define("Counts1","CountsRMS1.counts")
//     .Define("RMS1","CountsRMS1.rms")
//     .Define("CountsRMS2",GetCounts2(BremPhi4D2.GetPtr()),{"EBin","TBin","PBin"})
//     .Define("Counts2","CountsRMS2.counts")
//     .Define("RMS2","CountsRMS2.rms")
//     .Define("CountsRMS3",GetCounts2(BremPhi4D3.GetPtr()),{"EBin","TBin","PBin"})
//     .Define("Counts3","CountsRMS3.counts")
//     .Define("RMS3","CountsRMS3.rms")
//     .Define("CountsRMS4",GetCounts2(BremPhi4D4.GetPtr()),{"EBin","TBin","PBin"})
//     .Define("Counts4","CountsRMS4.counts")
//     .Define("RMS4","CountsRMS4.rms")
//     .Define("CountsRMS5",GetCounts2(BremPhi4D5.GetPtr()),{"EBin","TBin","PBin"})
//     .Define("Counts5","CountsRMS5.counts")
//     .Define("RMS5","CountsRMS5.rms")
//     .Define("CountsRMS6",GetCounts2(BremPhi4D6.GetPtr()),{"EBin","TBin","PBin"})
//     .Define("Counts6","CountsRMS6.counts")
//     .Define("RMS6","CountsRMS6.rms");


}
