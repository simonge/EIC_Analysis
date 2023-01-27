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
  
  ROOT::VecOps::RVec<float> operator()(const ROOT::VecOps::RVec<XYZVector>& pos, const ROOT::VecOps::RVec<XYZVector>& vec) {
    ROOT::VecOps::RVec<float> values;

    if(vec.size()>0 && pos.size()>0){
      
      *yP = pos[0].Y();
      *zP = pos[0].Z();
      *xV = vec[0].X();
      *yV = vec[0].Y();
      values = method->GetRegressionValues();
    }

    return values;
  };
  
private: 
  TMVA::MethodBase* method;
  float* yP;
  float* zP;
  float* xV;
  float* yV;
};

void RunNetwork( TString myMethodList = "" )
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
   TString tag    = "";
   //   TString tag    = "FrontWindow";
   //   TString prefix = "TMVARegression";
   TString ofile = "Reg_point-"+prefix+tag+".root";
 
   TString methodName = "DNN_CPU";//it->first + " method";
   TString weightfile = dir + prefix + "_" + "DNN_CPU" + ".weights.xml";
   reader->BookMVA( methodName, weightfile );
   TMVA::MethodBase* kl = dynamic_cast<TMVA::MethodBase*>(reader->FindMVA(methodName));

   
   // Prepare input tree (this must be replaced by your data source)
   // in this example, there is a toy tree with signal and one with background events
   // we'll later on use only the "signal" events for the test in this example.
   //
   TString fname = "/scratch/EIC/Analysis/temp"+tag+".root";
   ROOT::RDataFrame df("temp",fname);
 
   TStopwatch sw;
   sw.Start();
   
   auto df2 = df.Filter("(Tag1_4||Tag2_4)&&iFilter")
     //     .Define("values",getPrediction(kl,&yP,&zP,&xV,&yV),{"cell_cut","cell_vector"})
     .Define("values",getPrediction(kl,&yP,&zP,&xV,&yV),{"real_cut","real_vector"})
     .Define("ePred","values[0]")
     .Define("tPred","values[1]")
     .Define("pPred","atan2(values[3],values[2])");

//   std::vector<RResultPtr*> hists;
   auto energy = df2.Histo2D({"energy","energy",100,0,18,100,0,18},     "ePred","eE");
   auto theta  = df2.Histo2D({"theta", "theta", 100,0,0.01,100,0,0.01}, "tPred","thetaE");
   auto phi    = df2.Histo2D({"phi",   "phi",   100,-TMath::Pi(),TMath::Pi(),100,-TMath::Pi(),TMath::Pi()}, "pPred","phiV");
      
   // --- Write histograms

   TFile *target  = new TFile( ofile,"RECREATE" );
   
   energy->Write();
   theta->Write();
   phi->Write();
   //   for (auto hist:hists) hist->Write();
   target->Close();

   ROOT::RDF::RSnapshotOptions opts;
   opts.fMode = "UPDATE";

   df2.Snapshot("predict",ofile,{"eE","thetaE","phiV","ePred","tPred","pPred"},opts);


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
   RunNetwork(methodList);
   return 0;
}
