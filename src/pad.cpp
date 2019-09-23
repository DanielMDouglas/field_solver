#include <iostream>
#include <fstream>

#include "pad.h"
#include "physics.h"
#include "boundary.h"

pad::pad(double total_time, std::vector <double> center_pos)
{
  int nSteps = total_time/dt;
  for ( int i = 0; i < nSteps; i++ ) {
    current.push_back(0);
  }
  center = center_pos;
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

  // std::cout << "induction event starting at t0 = " << t0 << '\n'
  // 	    << "on pad center x: " << center[0]
  // 	    << "\ty: " << center[1] << std::endl;
  
  // the sensitive pad in the weighting field is at 0,0,0,
  // so place the drifted charge relative to that
  // std::vector <double> relative_pos = init_position - center;
  std::vector <double> proj_pos = {init_position[0], init_position[1], 0};
  
  if ( mag(proj_pos - center) < 1. ) {
    driftPath -> shift(-1*center);
    
    std::vector <double> pSeries = ramo_induction(q, driftPath, weighting);

    driftPath -> shift(center);
    
    for ( int i = 0; i < pSeries.size(); i++ ) {
      current[i + startingInd] += pSeries[i];
    }
  }
  // else {
  //   std::cout << "Pad out of range!  No response!" << std::endl;
  // }
}

void pad::calculate_output()
{
  // bullshit numbers!
  double discriminator_threshold = 1e-18; // C
  double hold_time = 2; // us
  
  std::vector <double> collected_charge;
  double t = 0;
  double sum = 0;
  double hold_counter = 0;
  for ( int i = 0; i < current.size(); i++ ) {
    sum += current[i]*dt;
    collected_charge.push_back(sum);
    if ( sum >= discriminator_threshold ) {
      hold_counter += dt;
      if ( hold_counter >= hold_time ) {
	std::vector <double> hit = {t, sum};
	output.push_back(hit);
	std::cout << "ADC outputs " << sum << "\t at " << t << std::endl;
	sum = 0;
	hold_counter = 0;
      }
    }
    t += dt;
  }
}

void pad::print_current_to_file(std::string filename)
{
  std::ofstream outFile (filename.c_str());
  for ( uint i = 0; i < current.size(); i++ ) {
    double t = i*dt;
    outFile << t << ','
	    << current[i] << '\n';
  }
  outFile.close();
}

void pad::print_output_to_file(std::string filename)
{
  std::ofstream outFile;
  outFile.open (filename.c_str(), std::ofstream::app);
  for ( uint i = 0; i < output.size(); i++ ) {
    std::cout << output[i][1] << std::endl;
    outFile << output[i][0] << ','
	    << output[i][1] << ','
	    << center[0] << ','
	    << center[1] << '\n';
  }
  outFile.close();
}
