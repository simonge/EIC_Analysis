
 
// need to add the current directory (from where we are running this macro)
// to the include path for Cling
R__ADD_INCLUDE_PATH($PWD)
#include <cmath>
#include "trainedModelRegression.hxx"
#include "TMVA/SOFIEHelpers.hxx"
#include "TMVA/RSofieReader.hxx"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodBase.h"
#include "ROOT/RDataFrame.hxx"
#include "TMVA/RInferenceUtils.hxx"

#include "TCanvas.h"
using namespace ROOT::Math;

 
using namespace TMVA::Experimental;
 
void RunPyNetworkE(){
 
   std::string inputName = "/scratch/EIC/Analysis/temp.root";
 
   ROOT::EnableImplicitMT(10);
 
   ROOT::RDataFrame df("temp", inputName);
   auto df1 = df.Filter("(Tag1_4||Tag2_4)&&iFilter")
     .Define("Pos0","real_cut[0]")
     .Define("PosY","(float)Pos0.y()")
     .Define("PosZ","(float)Pos0.z()")
     .Define("Vec0","real_vector[0]")
     .Define("VecX","(float)Vec0.x()")
     .Define("VecY","(float)Vec0.y()");

   int nslots = df1.GetNSlots();
   std::cout << "Running using " << nslots << " threads" << std::endl;
   RSofieReader model("trainedModelRegression.h5");
   //RSofieReader model("modelRegression.h5");

   std::vector<float> testinput = {0.5};
   auto output = model.Compute(testinput);

   std::cout << output[0] << std::endl;

   df1 = df1.Define("DNN_Value", Compute<1,float>(model),
			{"VecX"});
//    df1 = df1.Define("DNN_Value", Compute<4,float>(model),
// 			{"PosY","PosZ","VecX","VecY"});
//    df1 = df1.DefineSlot("DNN_Value", fun,
//  			{"VecY","VecX","PosZ","PosY"});
//    df1 = df1.DefineSlot("DNN_Value", SofieFunctor<4, TMVA_SOFIE_trainedModelRegression::Session>(nslots),
// 			{"PosY","PosZ","VecX","VecY"});
   auto h1 = df1.Histo2D({"h_sig", "", 400, 0, 18, 400, 0, 18}, "DNN_Value","eE");
 
   

   df1.Snapshot("temp","DNNout.root",{"eE","DNN_Value","PosY","PosZ","VecX","VecY"});
   
   auto c1 = new TCanvas();
//    gStyle->SetOptStat(0);
 
   h1->DrawClone("colz");
   c1->BuildLegend();
 
}
