
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
 
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
 
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/TMVARegGui.h"
 
 
using namespace TMVA;
 
void FPRegression( TString myMethodList = "" )
{

   TString fname = "/scratch/EIC/Analysis/FP_Tagger_Test_events3.root";
   TString outfileName( "/scratch/EIC/Results/ML-Out/FP_Reg_qr_out18x275_Cheat_E_NoVertex.root" );

    //---------------------------------------------------------------
   // This loads the library
   TMVA::Tools::Instance();
 
   ROOT::EnableImplicitMT(8);
 
   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;
 
   // Neural Network
   Use["DNN_CPU"]         = 1;

   // ---------------------------------------------------------------
 
   std::cout << std::endl;
   std::cout << "==> Start TMVARegression" << std::endl;
 
   // Select methods (don't look at this code - not of interest)
   if (myMethodList != "") {
      for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;
 
      std::vector<TString> mlist = gTools().SplitString( myMethodList, ',' );
      for (UInt_t i=0; i<mlist.size(); i++) {
         std::string regMethod(mlist[i].Data());
 
         if (Use.find(regMethod) == Use.end()) {
            std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
            for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
            std::cout << std::endl;
            return;
         }
         Use[regMethod] = 1;
      }
   }
 
   // --------------------------------------------------------------------------------------------------
 
   // Here the preparation phase begins
 
   // Create a new root output file
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
 
   // Create the factory object. Later you can choose the methods
   TMVA::Factory *factory = new TMVA::Factory( "TMVARegression", outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;G;N;N+G;P+N+G:AnalysisType=Regression" );
 
 
   TMVA::DataLoader *dataloader=new TMVA::DataLoader("dataset");
   // If you wish to modify default settings
   // (please check "src/Config.h" to see all available global options)
   //
   //     (TMVA::gConfig().GetVariablePlotting()).fTimesRMS = 8.0;
   //     (TMVA::gConfig().GetIONames()).fWeightFileDir = "myWeightDirectory";
 
   // Define the input variables that shall be used for the MVA training
   // note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
   // [all types of expressions that can also be parsed by TTree::Draw( "expression" )]
//     dataloader->AddVariable( "TMath::Floor(Xpix11[0]/2)", "Xpix11", "units", 'mm' );
//     dataloader->AddVariable( "TMath::Floor(Ypix11[0]/2)", "Ypix11", "units", 'mm' );
//     dataloader->AddVariable( "TMath::Floor(Xpix12[0]/2)", "Xpix12", "units", 'mm' );
//     dataloader->AddVariable( "TMath::Floor(Ypix12[0]/2)", "Ypix12", "units", 'mm' );
  dataloader->AddVariable( "real_position_x[0]", "real_position_x", "units", 'mm' );
  dataloader->AddVariable( "real_position_y[0]", "real_position_y", "units", 'mm' );
  dataloader->AddVariable( "real_position_z[0]", "real_position_z", "units", 'mm' );
  dataloader->AddVariable( "real_vector[0].fX", "real_vector.x", "units", 'mm' );
  dataloader->AddVariable( "real_vector[0].fY", "real_vector.y", "units", 'mm' );
  dataloader->AddVariable( "real_vector[0].fZ", "real_vector.z", "units", 'mm' );
//   dataloader->AddVariable( "vertex.fX", "VertexX", "units", 'mm' );
//   dataloader->AddVariable( "vertex.fY", "VertexY", "units", 'mm' );
//   dataloader->AddVariable( "vertex.fZ", "VertexZ", "units", 'mm' );
 
   // You can add so-called "Spectator variables", which are not used in the MVA training,
   // but will appear in the final "TestTree" produced by TMVA. This TestTree will contain the
   // input variables, the response values of all trained MVAs, and the spectator variables
  dataloader->AddSpectator( "eE" );
  dataloader->AddSpectator( "logQ2" );
    dataloader->AddSpectator( "qPhi" );
    dataloader->AddSpectator( "qTheta" );
    //    dataloader->AddSpectator( "logQ2",  "Spectator 1", "units", 'F' );
//    dataloader->AddSpectator( "spec2:=var1*3",  "Spectator 2", "units", 'F' );
 
   // Add the variable carrying the regression target
   // dataloader->AddTarget( "Xpix11" );
//    dataloader->AddTarget( "Ypix11" );
    dataloader->AddTarget( "eE" );
    //         dataloader->AddTarget( "cos(qPhi)" );
    //   dataloader->AddTarget( "sin(qPhi)" );
//     dataloader->AddTarget( "qTheta" );
    //dataloader->AddTarget( "logQ2" );
//    dataloader->AddTarget( "qx" );
//    dataloader->AddTarget( "qy" );
//    dataloader->AddTarget( "qz" );
 

   // It is also possible to declare additional targets for multi-dimensional regression, ie:
   //     factory->AddTarget( "fvalue2" );
   // BUT: this is currently ONLY implemented for MLP
 
   // Read training and test data (see TMVAClassification for reading ASCII files)
   // load the signal and background event samples from ROOT trees
   TFile *input(0);
   //TString fname = "qr_out18x275-2.root";
   //TString fname = "/scratch/EIC/Analysis/psi2s_10_100.root";
   if (!gSystem->AccessPathName( fname )) {
      input = TFile::Open( fname ); // check if file in local directory exists
   }
   if (!input) {
      std::cout << "ERROR: could not open data file" << std::endl;
      exit(1);
   }
   std::cout << "--- TMVARegression           : Using input file: " << input->GetName() << std::endl;
 
   // Register the regression tree
 
   TTree *regTree1 = (TTree*)input->Get("detectors");
//    TTree *regTree2 = (TTree*)input->Get("detector2");
 
   // global event weights per tree (see below for setting event-wise weights)
   Double_t regWeight  = 1.0;
 
   // You can add an arbitrary number of regression trees
   dataloader->AddRegressionTree( regTree1, regWeight );
 
   // This would set individual event weights (the variables defined in the
   // expression need to exist in the original TTree)
   //dataloader->SetWeightExpression( "var1", "Regression" );
 
   // Apply additional cuts on the signal and background samples (can be different)
   TCut mycut = "eE<18 && eE>2"; 
 
   dataloader->PrepareTrainingAndTestTree(mycut,"nTrain_Regression=0:nTest_Regression=0:SplitMode=Random:NormMode=NumEvents:!V");

 
 
   if (Use["DNN_CPU"] or Use["DNN_GPU"]) {
 
     //     TString layoutString("Layout=TANH|200,TANH|100,TANH|100,LINEAR");
     //     TString layoutString("Layout=TANH|50,TANH|30,TANH|20,TANH|20,LINEAR");
     TString layoutString("Layout=TANH|64,TANH|64,TANH|32,TANH|16,TANH|8,TANH|4,LINEAR");
     //     TString layoutString("Layout=TANH|400,TANH|8,TANH|16,TANH|8,TANH|4,LINEAR");
     // TString layoutString("Layout=LINEAR|6,TANH|6,LINEAR");
 
 
      TString trainingStrategyString("TrainingStrategy=");
 
      trainingStrategyString +="LearningRate=1e-4,Momentum=0.2,ConvergenceSteps=2000,BatchSize=400,TestRepetitions=1,WeightDecay=0,Regularization=L2,Optimizer=Adam";
      //      trainingStrategyString +="LearningRate=1e-3,Momentum=0.3,ConvergenceSteps=10,BatchSize=100,TestRepetitions=1,WeightDecay=0,Regularization=L2,Optimizer=Adam";
 
      TString nnOptions("!H:V:ErrorStrategy=SUMOFSQUARES:VarTransform=G:WeightInitialization=XAVIERUNIFORM:Architecture=GPU");
      //      TString nnOptions("!H:V:ErrorStrategy=SUMOFSQUARES:VarTransform=G+N:WeightInitialization=XAVIERUNIFORM:Architecture=GPU");
      nnOptions.Append(":");
      nnOptions.Append(layoutString);
      nnOptions.Append(":");
      nnOptions.Append(trainingStrategyString);
 
      factory->BookMethod(dataloader, TMVA::Types::kDL, "DNN_CPU", nnOptions); // NN
   }

   // --------------------------------------------------------------------------------------------------
 
   // Now you can tell the factory to train, test, and evaluate the MVAs
 
   // Train MVAs using the set of training events
   factory->TrainAllMethods();
 
   // Evaluate all MVAs using the set of test events
   factory->TestAllMethods();
 
   // Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();
 
   // --------------------------------------------------------------
 
   // Save the output
   outputFile->Close();
 
   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVARegression is done!" << std::endl;
 
   delete factory;
   delete dataloader;
 
   // Launch the GUI for the root macros
   if (!gROOT->IsBatch()) TMVA::TMVARegGui( outfileName );
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
   FPRegression(methodList);
   return 0;
}
