#include <iostream>
#include <fstream>

#include "pad.h"
#include "physics.h"
#include "boundary.h"

pad::pad(double total_time, std::vector <double> center_pos)
{
  int nSteps = total_time/dt;
  for ( int i = 0; i < nSteps; i++ ) {
    response.push_back(0);
  }
  center = center_pos;
}

void pad::add_response(double q,
		       double t0,
		       std::vector <double> charge_position,
		       field <double> * potential,
		       field <double> * weighting,
		       boundary bound)
{
  int startingInd = t0/dt;

  std::cout << "induction event starting at t0 = " << t0 << '\n'
	    << "on pad center x: " << center[0]
	    << "\ty: " << center[1] << std::endl;
  
  // the sensitive pad in the weighting field is at 0,0,0,
  // so place the drifted charge relative to that
  std::vector <double> relative_pos = charge_position - center;
  std::vector <double> proj_pos = {charge_position[0], charge_position[1], 0};
  
  if ( mag( proj_pos - center) < 1.4 ) {
    path * driftPath;
    drift_path(relative_pos, potential, bound, driftPath);

    // driftPath -> print_to_file("path.dat");
  
    std::vector <double> pSeries = ramo_induction(q, driftPath, weighting);

    delete driftPath;
    
    std::ofstream outFile("pSeries.dat");
    for ( uint i = 0; i < pSeries.size(); i++ ) {
      outFile << i*dt << ','
	      << pSeries[i] << '\n';
    }
    outFile.close();
    
    for ( int i = 0; i < pSeries.size(); i++ ) {
      response[i + startingInd] += pSeries[i];
    }
  }
  else {
    std::cout << "Pad out of range!  No response!" << std::endl;
  }
}

void pad::add_response(double q,
		       double t0,
		       std::vector <double> init_position,
		       path * driftPath,
		       field <double> * potential,
		       field <double> * weighting,
		       boundary bound)
{
  int startingInd = t0/dt;

  std::cout << "induction event starting at t0 = " << t0 << '\n'
	    << "on pad center x: " << center[0]
	    << "\ty: " << center[1] << std::endl;
  
  // the sensitive pad in the weighting field is at 0,0,0,
  // so place the drifted charge relative to that
  // std::vector <double> relative_pos = init_position - center;
  std::vector <double> proj_pos = {init_position[0], init_position[1], 0};
  
  if ( mag(proj_pos - center) < 1. ) {
    driftPath -> shift(-1*center);
    
    std::vector <double> pSeries = ramo_induction(q, driftPath, weighting);

    driftPath -> shift(center);
    
    for ( int i = 0; i < pSeries.size(); i++ ) {
      response[i + startingInd] += pSeries[i];
    }
  }
  else {
    std::cout << "Pad out of range!  No response!" << std::endl;
  }
}

void pad::print_to_file(std::string filename)
{
  std::ofstream outFile (filename.c_str());
  for ( uint i = 0; i < response.size(); i++ ) {
    double t = i*dt;
    outFile << t << ','
	    << response[i] << '\n';
  }
  outFile.close();
}
