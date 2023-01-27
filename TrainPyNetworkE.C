
 
// need to add the current directory (from where we are running this macro)
// to the include path for Cling
 
#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"

#include "TMVA/MethodDNN.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/TMVARegGui.h"
#include "TMVA/PyMethodBase.h"
using namespace TMVA::Experimental;
 
void TrainPyNetworkE(){
 

   ROOT::EnableImplicitMT(12);

   gSystem->Setenv("KERAS_BACKEND", "tensorflow");
   // for using Keras
   TMVA::PyMethodBase::PyInitialize();

   TString outfileName = "pyOut.root";
   TString inputName   = "/scratch/EIC/Analysis/temp.root";

   auto outputFile = TFile::Open( outfileName, "RECREATE" );
 
   TString typeName = "TEST";
   
   TMVA::Factory factory( typeName,outputFile,"!V:ROC:!Silent:Color:DrawProgressBar:AnalysisType=Regression" );
   //   TMVA::Factory factory( typeName, outputFile,"!V:ROC:!Silent:Color:DrawProgressBar:AnalysisType=Regression" );


   TString dataFolderName = "dataset"; 
   TMVA::DataLoader *dataloader=new TMVA::DataLoader(dataFolderName);

   // Variables with cell binning
   //    dataloader->AddVariable( "cell_cut[0].fY", "real_position_y", "units", 'F' );
   //    dataloader->AddVariable( "cell_cut[0].fZ", "real_position_z", "units", 'F' );
   //    dataloader->AddVariable( "cell_vector[0].fX", "real_vector.x", "units", 'F' );
   //    dataloader->AddVariable( "cell_vector[0].fY", "real_vector.y", "units", 'F' );

   //   Variables with no cell binning
      dataloader->AddVariable( "(float)real_cut[0].fCoordinates.fY",    "real_position_y", "units", 'F' );
    dataloader->AddVariable( "(float)real_cut[0].fCoordinates.fZ",    "real_position_z", "units", 'F' );
   dataloader->AddVariable( "real_vector[0].fCoordinates.fX", "real_vector_x",   "units", 'F' );
    dataloader->AddVariable( "(float)real_vector[0].fCoordinates.fY", "real_vector_y",   "units", 'F' );

  dataloader->AddTarget( "eE" );
  //  dataloader->AddTarget( "thetaE" );
  //   dataloader->AddTarget( "cos(phiV)" );
  //   dataloader->AddTarget( "sin(phiV)" );


  TFile *input(0);
  input = TFile::Open( inputName ); // check if file in local directory exists
  
  std::cout << "--- TMVARegression           : Using input file: " << input->GetName() << std::endl;
  
  TTree *regTree1 = (TTree*)input->Get("temp");
  regTree1->SetEntries(200000);
  dataloader->AddRegressionTree( regTree1 );

   TCut mycut = "(Tag1_4||Tag2_4)&&iFilter"; 
   //TCut mycut = ""; 
//    dataloader->AddTree( regTree1, "Regression", regWeight,mycut );
 
   dataloader->PrepareTrainingAndTestTree(mycut,"nTrain_Regression=0:nTest_Regression=0:SplitMode=Random:NormMode=NumEvents:!V");

  factory.BookMethod( dataloader, TMVA::Types::kPyKeras, "PyKeras",
		     "H:!V:VarTransform=D:FilenameModel=modelRegression.h5:tf.keras:"
		     "FilenameTrainedModel=trainedModelRegression.h5:NumEpochs=2000:BatchSize=1024:"
		     "GpuOptions=allow_growth=True");

  
  factory.TrainAllMethods();
   
  factory.TestAllMethods();
  
  factory.EvaluateAllMethods();

  outputFile->Close();
   
}
