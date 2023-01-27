#include <cstdlib>
#include <vector>
#include <iostream>
 
using namespace TMVA::Experimental;
                

void ConvertToSOFIE( std::string modelFile = "trainedModelRegression.h5" )
{

  bool verbose=0;

     //Parsing the saved Keras .h5 file into RModel object
  SOFIE::RModel model = SOFIE::PyKeras::Parse(modelFile);
 
 
  //Generating inference code
  model.Generate();
  // generate output header. By default it will be modelName.hxx
  model.OutputGenerated();
 
  if(!verbose) return;

  //Printing required input tensors
  std::cout<<"\n\n";
  model.PrintRequiredInputTensors();
  
  //Printing initialized tensors (weights)
  std::cout<<"\n\n";
  model.PrintInitializedTensors();
 
  //Printing intermediate tensors
  std::cout<<"\n\n";
  model.PrintIntermediateTensors();
  
  //Printing generated inference code
  std::cout<<"\n\n";
  model.PrintGenerated();
  
  return;
}
