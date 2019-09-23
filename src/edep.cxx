#include <TChain.h>
#include <string>

const std::string dataFileName = "../data/aho.root";

int main()
{
  TChain * ch0 = new TChain("sparse3d_mcst_tree");
  TChain * ch1 = new TChain("sparse3d_mcst_dx_tree");
  TChain * ch2 = new TChain("sparse3d_mcst_dedx_tree");  

  for ( TChain * ch: {ch0, ch1, ch2} ) {
    ch -> AddFile(dataFileName.c_str());
    ch -> GetEntry(0);
  }
  
  return 0;
}
