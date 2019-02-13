#include <iostream>
#include <string>

// #include <TROOT.h>
// #include <TFile.h>
// #include <TGeoManager.h>
// #include <TGLViewer.h>
// #include <TPad.h>
// #include <TClass.h>

// #include <G4Box.hh>
// #include <G4PVPlacement.hh>
// #include <G4SystemOfUnits.hh>
#include <G4GDMLParser.hh>

void geometry() {
  std::string filename = "/home/dan/studies/detector.gdml";
  
  std::string fieldShellName = "volFieldShell";
  std::string cathodeName = "volCathode";
  std::string leftPixelPlaneName = "volLeftPixelPlane";
  std::string rightPixelPlaneName = "volRightPixelPlane";  

  G4GDMLParser parser;
  parser.Read(filename);  
}
