
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
 
void TaggerRegression( TString myMethodList = "" )
{

   TString fname = "/scratch/EIC/Analysis/temp.root";
   //   TString fname = "/scratch/EIC/Analysis/FP_Tagger_Test_events3.root";
   //   TString outfileName( "/scratch/EIC/Results/ML-Out/test_cell_ETP.root" );
   TString tag = "4layerP";
   TString outfileName( "/scratch/EIC/Results/ML-Out/"+tag+"_real_ETP.root" );

    //---------------------------------------------------------------
   // This loads the library
   TMVA::Tools::Instance();
 
   ROOT::EnableImplicitMT(8);
  
   // --------------------------------------------------------------------------------------------------
    // Here the preparation phase begins
    // Create a new root output file
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
 
   // Create the factory object. Later you can choose the methods
   TMVA::Factory *factory = new TMVA::Factory( "RealHits"+tag, outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:AnalysisType=Regression" );
   //   TMVA::Factory *factory = new TMVA::Factory( "TMVARegression", outputFile,
   //                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;G;N;N+G;P+N+G:AnalysisType=Regression" );
 
 
   TMVA::DataLoader *dataloader=new TMVA::DataLoader("dataset");
   // If you wish to modify default settings
   // (please check "src/Config.h" to see all available global options)
   //
   //     (TMVA::gConfig().GetVariablePlotting()).fTimesRMS = 8.0;
   //     (TMVA::gConfig().GetIONames()).fWeightFileDir = "myWeightDirectory";
 
   // Define the input variables that shall be used for the MVA training
   // note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
   // [all types of expressions that can also be parsed by TTree::Draw( "expression" )]


   // Variables with cell binning
   //    dataloader->AddVariable( "cell_cut[0].fY", "real_position_y", "units", 'F' );
   //    dataloader->AddVariable( "cell_cut[0].fZ", "real_position_z", "units", 'F' );
   //    dataloader->AddVariable( "cell_vector[0].fX", "real_vector.x", "units", 'F' );
   //    dataloader->AddVariable( "cell_vector[0].fY", "real_vector.y", "units", 'F' );

   // Variables with no cell binning
   dataloader->AddVariable( "real_cut[0].fCoordinates.fY",    "real_position_y", "units", 'F' );
   dataloader->AddVariable( "real_cut[0].fCoordinates.fZ",    "real_position_z", "units", 'F' );
   dataloader->AddVariable( "real_vector[0].fCoordinates.fX", "real_vector.x",   "units", 'F' );
   dataloader->AddVariable( "real_vector[0].fCoordinates.fY", "real_vector.y",   "units", 'F' );


   //dataloader->AddVariable( "real_vector[0].fZ", "real_vector.z", "units", 'F' );

//   dataloader->AddVariable( "vertex.fX", "VertexX", "units", 'mm' );
//   dataloader->AddVariable( "vertex.fY", "VertexY", "units", 'mm' );
//   dataloader->AddVariable( "vertex.fZ", "VertexZ", "units", 'mm' );
 
   // You can add so-called "Spectator variables", which are not used in the MVA training,
   // but will appear in the final "TestTree" produced by TMVA. This TestTree will contain the
   // input variables, the response values of all trained MVAs, and the spectator variables
  dataloader->AddSpectator( "eE" );
  dataloader->AddSpectator( "logQ2" );
//   dataloader->AddSpectator( "qPhi" );
//   dataloader->AddSpectator( "qTheta" );
  //    dataloader->AddSpectator( "logQ2",  "Spectator 1", "units", 'F' );
  //    dataloader->AddSpectator( "spec2:=var1*3",  "Spectator 2", "units", 'F' );
  
  // Add the variable carrying the regression target
  dataloader->AddTarget( "eE" );
  dataloader->AddTarget( "thetaE" );
  dataloader->AddTarget( "cos(phiV)" );
  dataloader->AddTarget( "sin(phiV)" );
 

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
 
   TTree *regTree1 = (TTree*)input->Get("temp");
   regTree1->SetEntries(2000000);

   // global event weights per tree (see below for setting event-wise weights)
   Double_t regWeight  = 1.0;
 
   // You can add an arbitrary number of regression trees
   dataloader->AddRegressionTree( regTree1, regWeight );
 
   // This would set individual event weights (the variables defined in the
   // expression need to exist in the original TTree)
   //dataloader->SetWeightExpression( "var1", "Regression" );
   // Apply additional cuts on the signal and background samples (can be different)
   TCut mycut = "(Tag1_4||Tag2_4)&&iFilter"; 
   //TCut mycut = ""; 
//    dataloader->AddTree( regTree1, "Regression", regWeight,mycut );
 
   dataloader->PrepareTrainingAndTestTree(mycut,"nTrain_Regression=0:nTest_Regression=0:SplitMode=Random:NormMode=NumEvents:!V");

   //  TString layoutString("Layout=TANH|64,TANH|32,TANH|16,LINEAR");//USED 3layer
   TString layoutString("Layout=TANH|200,TANH|200,TANH|100,TANH|50,LINEAR");//USED 1layer
   //   TString layoutString("Layout=TANH|100,TANH|50,TANH|25,TANH|16,LINEAR");//USED 2layer
   //   TString layoutString("Layout=TANH|64,TANH|32,TANH|32,LINEAR");//USED 2layer
   //TString layoutString("Layout=TANH|1024,TANH|64,TANH|32,TANH|16,LINEAR");//USED 0.145751
     
   TString trainingStrategyString("TrainingStrategy=");
   
   //   trainingStrategyString +="LearningRate=1e-4,Momentum=0.01,MaxEpochs=4000,ConvergenceSteps=1000,BatchSize=50,TestRepetitions=1,Regularization=L2,Optimizer=Adam";
   trainingStrategyString +="LearningRate=2e-4,Momentum=0.0,MaxEpochs=4000,ConvergenceSteps=200,BatchSize=200,TestRepetitions=5,Regularization=L2,Optimizer=Adam";
   //   trainingStrategyString +="LearningRate=1e-4,Momentum=0.01,DropConfig=0.0+0.0+0.0+0.0+0.0,ConvergenceSteps=1000,BatchSize=200,TestRepetitions=1,WeightDecay=0.0,Regularization=L2,Optimizer=Adam";
   //      trainingStrategyString +="LearningRate=1e-3,Momentum=0.3,ConvergenceSteps=10,BatchSize=100,TestRepetitions=1,WeightDecay=0,Regularization=L2,Optimizer=Adam";
   
   TString nnOptions("!H:V:ErrorStrategy=SUMOFSQUARES:VarTransform=G,N(thetaE),P:WeightInitialization=XAVIERUNIFORM:Architecture=GPU");
   //   TString nnOptions("!H:V:ErrorStrategy=SUMOFSQUARES:VarTransform=G+N:WeightInitialization=XAVIERUNIFORM:Architecture=GPU");
   //      TString nnOptions("!H:V:ErrorStrategy=SUMOFSQUARES:VarTransform=G+N:WeightInitialization=XAVIERUNIFORM:Architecture=GPU");
   nnOptions.Append(":");
   nnOptions.Append(layoutString);
   nnOptions.Append(":");
   nnOptions.Append(trainingStrategyString);
   
   factory->BookMethod(dataloader, TMVA::Types::kDL, "DNN_CPU", nnOptions); // NN

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
   TaggerRegression(methodList);
   return 0;
}
