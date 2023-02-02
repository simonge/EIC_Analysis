#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>
 
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TH2F.h"
#include "TMath.h"
#include "Math/Vector3D.h"
using namespace ROOT::Math;
 
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodBase.h"
#include "ROOT/RDataFrame.hxx"
 
using namespace TMVA;
 
//-----------------------------------------------------------------------------------------
// Get Prediction Functor
//-----------------------------------------------------------------------------------------
struct getPrediction{
  getPrediction(TMVA::MethodBase* kl,float* yp,float* zp,float* xv,float* yv) : method(kl), yP(yp), zP(zp), xV(xv), yV(yv){}
  
  ROOT::VecOps::RVec<float> operator()(const XYZVector& pos, const XYZVector& vec) {
    ROOT::VecOps::RVec<float> values;
      
    *yP = pos.Y();
    *zP = pos.Z();
    *xV = vec.X();
    *yV = vec.Y();
    values = method->GetRegressionValues();
    

    return values;
  };
  
private: 
  TMVA::MethodBase* method;
  float* yP;
  float* zP;
  float* xV;
  float* yV;
};

void RunNetworkNew( TString myMethodList = "" )
{
   //---------------------------------------------------------------
   // This loads the library
   TMVA::Tools::Instance();
 
   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;
 
   Use["DNN_CPU"] = 1;
 
   // --------------------------------------------------------------------------------------------------
 
   // --- Create the Reader object
 
   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );
 
   // Create a set of variables and declare them to the reader
   // - the variable names MUST corresponds in name and type to those given in the weight file(s) used
   Float_t yP, zP, xV, yV;
//    reader->AddVariable( "cell_cut[0].fCoordinates.fY",    &yP     );
//    reader->AddVariable( "cell_cut[0].fCoordinates.fZ",    &zP     );
//    reader->AddVariable( "cell_vector[0].fCoordinates.fX", &xV     );
//    reader->AddVariable( "cell_vector[0].fCoordinates.fY", &yV     );

   reader->AddVariable( "real_cut[0].fCoordinates.fY",    &yP     );
   reader->AddVariable( "real_cut[0].fCoordinates.fZ",    &zP     );
   reader->AddVariable( "real_vector[0].fCoordinates.fX", &xV     );
   reader->AddVariable( "real_vector[0].fCoordinates.fY", &yV     );
 
   // Spectator variables declared in the training have to be added to the reader, too
   Float_t  eE, logQ2;
   reader->AddSpectator( "eE",     &eE );
   reader->AddSpectator( "logQ2",  &logQ2 );
 
   // --- Book the MVA methods
 
   TString dir    = "dataset/weights/";
   TString prefix = "RealHits4layer";
   TString tag    = "tempClusterBrems";
   //   TString tag    = "tempClusterQR";
   //TString tag    = "clusterFrontWindow";
   //   TString tag    = "FrontWindow";
   //   TString prefix = "TMVARegression";
   TString ofile = "Reg_point2-"+prefix+tag+".root";
 
   TString methodName = "DNN_CPU";//it->first + " method";
   TString weightfile = dir + prefix + "_" + "DNN_CPU" + ".weights.xml";
   reader->BookMVA( methodName, weightfile );
   TMVA::MethodBase* kl = dynamic_cast<TMVA::MethodBase*>(reader->FindMVA(methodName));

   
   // Prepare input tree (this must be replaced by your data source)
   // in this example, there is a toy tree with signal and one with background events
   // we'll later on use only the "signal" events for the test in this example.
   //
   TString fname = "/scratch/EIC/Analysis/"+tag+".root";
   ROOT::RDataFrame df("temp",fname);
 
   TStopwatch sw;
   sw.Start();
   
   auto df2 = df//.Filter("(Tag1_4||Tag2_4)&&iFilter")
     //     .Define("values",getPrediction(kl,&yP,&zP,&xV,&yV),{"cell_cut","cell_vector"})
     .Filter("fit_cut.size()>0")
     .Define("indexMin","ArgMin(fit_chi2)")  //Fitted position
     .Define("fit_minchi2","fit_chi2[indexMin]")  //Fitted position
     .Filter("fit_minchi2<0.01")  //Fitted position
     .Define("Rposition","real_cut[0]") //Exact position
     .Define("Rvector",  "real_vector[0]")
     .Define("Cposition","cell_cut[0]") //Cell position
     .Define("Cvector",  "cell_vector[0]")
     .Define("Fposition","fit_cut[indexMin]")  //Fitted position
     .Define("Fvector","fit_vector[indexMin]")
     .Define("xmin","fit_minX[indexMin]")
     .Define("xmax","fit_maxX[indexMin]")
     .Define("ymin","fit_minY[indexMin]")
     .Define("ymax","fit_maxY[indexMin]")
     .Define("Rvalues",getPrediction(kl,&yP,&zP,&xV,&yV),{"Rposition","Rvector"})
     .Define("Cvalues",getPrediction(kl,&yP,&zP,&xV,&yV),{"Cposition","Cvector"})
     .Define("Fvalues",getPrediction(kl,&yP,&zP,&xV,&yV),{"Fposition","Fvector"})
     .Define("ePredR","Rvalues[0]")
     .Define("tPredR","Rvalues[1]")
     .Define("pPredR","atan2(Rvalues[3],Rvalues[2])")
     .Define("ePredC","Cvalues[0]")
     .Define("tPredC","Cvalues[1]")
     .Define("pPredC","atan2(Cvalues[3],Cvalues[2])")
     .Define("ePredF","Fvalues[0]")
     .Define("tPredF","Fvalues[1]")
     .Define("pPredF","atan2(Fvalues[3],Fvalues[2])");
//   std::vector<RResultPtr*> hists;
   auto Renergy = df2.Histo2D({"Renergy","Renergy",100,0,18,100,0,18},     "ePredR","eE");
   auto Rtheta  = df2.Histo2D({"Rtheta", "Rtheta", 100,0,0.01,100,0,0.01}, "tPredR","thetaE");
   auto Rphi    = df2.Histo2D({"Rphi",   "Rphi",   100,-TMath::Pi(),TMath::Pi(),100,-TMath::Pi(),TMath::Pi()}, "pPredR","phiV");
   auto Cenergy = df2.Histo2D({"Cenergy","Cenergy",100,0,18,100,0,18},     "ePredC","eE");
   auto Ctheta  = df2.Histo2D({"Ctheta", "Ctheta", 100,0,0.01,100,0,0.01}, "tPredC","thetaE");
   auto Cphi    = df2.Histo2D({"Cphi",   "Cphi",   100,-TMath::Pi(),TMath::Pi(),100,-TMath::Pi(),TMath::Pi()}, "pPredC","phiV");
   auto Fenergy = df2.Histo2D({"Fenergy","Fenergy",100,0,18,100,0,18},     "ePredF","eE");
   auto Ftheta  = df2.Histo2D({"Ftheta", "Ftheta", 100,0,0.01,100,0,0.01}, "tPredF","thetaE");
   auto Fphi    = df2.Histo2D({"Fphi",   "Fphi",   100,-TMath::Pi(),TMath::Pi(),100,-TMath::Pi(),TMath::Pi()}, "pPredF","phiV");
      
   // --- Write histograms

   TFile *target  = new TFile( ofile,"RECREATE" );
   
   Renergy->Write();
   Rtheta->Write();
   Rphi->Write();
   Cenergy->Write();
   Ctheta->Write();
   Cphi->Write();
   Fenergy->Write();
   Ftheta->Write();
   Fphi->Write();
   //   for (auto hist:hists) hist->Write();
   target->Close();

   ROOT::RDF::RSnapshotOptions opts;
   opts.fMode = "UPDATE";

   df2.Snapshot("predict",ofile,{"eE","thetaE","phiV","ePredR","tPredR","pPredR","ePredC","tPredC","pPredC","ePredF","tPredF","pPredF","xmin","xmax","ymin","ymax","fit_minchi2"},opts);


   sw.Stop();
   std::cout << "--- End of event loop: "; sw.Print();
   
   std::cout << "--- Created root file: \"" << target->GetName()
             << "\" containing the MVA output histograms" << std::endl;
 
   delete reader;
 
   std::cout << "==> TMVARegressionApplication is done!" << std::endl << std::endl;
}
 
int main( int argc, char** argv )
{
   // Select methods (don't look at this code - not of interest)
   TString methodList;
   for (int i=1; i<argc; i++) {
      TString regMethod(argv[i]);
      if(regMethod=="-b" || regMethod=="--batch") continue;
      if (!methodList.IsNull()) methodList += TString(",");
      methodList += regMethod;
   }
   RunNetworkNew(methodList);
   return 0;
}
