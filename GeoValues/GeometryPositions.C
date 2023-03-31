
#include "DD4hep/Detector.h"
#include "DDRec/CellIDPositionConverter.h"
 
using namespace ROOT::Math;
  
//-----------------------------------------------------------------------------------------
// Main Analysis Call
//-----------------------------------------------------------------------------------------

void GeometryPositions(std::string compactName = "/home/simon/EIC/epic/epic_ip6.xml"){

  // -------------------------
  // Get the DD4hep instance
  // Load the compact XML file
  // Initialize the position converter tool
  dd4hep::Detector& detector = dd4hep::Detector::getInstance();
  detector.fromXML(compactName);
  dd4hep::rec::CellIDPositionConverter cellid_converter(detector);
  auto taggers = detector.detector("BackwardsTaggerStation");
//   auto place = tagger.placement().matrix();
//   std::cout << place << std::endl;
  auto taggersPos  = taggers.placement().position();
  TGeoRotation    mat = taggers.placement().matrix();
  TGeoTranslation pos = taggers.placement().matrix();

  auto pos2 = pos*mat;

  auto tagger1 = taggers.child("Tagger1");
  auto tagger2 = taggers.child("Tagger2");
  auto tag1Pos = taggersPos+tagger1.placement().position();
  
  for(auto [tname,tagger]:taggers.children()){
    for(auto [cname,child]:tagger.children()){
      auto align =  child.nominal().worldTransformation();
      align.Print();
    }

  }
  
}

//--------------------------------------------------------------------------------------------
