#include <cmath>
#include <iostream>

#include <TF2.h>

#include "physics.h"
#include "volume.h"

std::vector <double> E(std::vector <double> pos, field <double> * V)
{
  return -1*(V -> interpolate_grad(pos));
}

std::vector <double> driftV(std::vector <double> eField, double temperature)
{
  // Drift Velocity as a function of Electric Field and LAr Temperature
  // from : W. Walkowiak, NIM A 449 (2000) 288-294
  //
  // Efield should have units of kV/cm
  // Temperature should have units of Kelvin

  double magE = mag(eField);
  std::vector <double> dir = norm(eField);
  if ( magE > 4.0 ) {
    magE = 4.0;
  }

  if ( temperature < 87.0 || temperature > 94.0 ) {
    std::cout << "DriftVelocity Warning!: Temperature value of "
	      << temperature
	      << " K is outside of range covered by drift velocity"
	      << " parameterization.  Returned value may not be"
	      << " correct"
	      << std::endl;
  }

  double tShift = -87.203 + temperature;
  double xFit = 0.0938163 - 0.0052563*tShift - 0.0001470*tShift*tShift;
  double uFit = 5.18406 + 0.01448*tShift - 0.003497*tShift*tShift - 0.000516*tShift*tShift*tShift;
  double vd;

  // Icarus Parameter set
  // Use as default
  double P1 = -0.04640; // K^-1
  double P2 = 0.01712; // K^-1
  double P3 = 1.88125; // (kV/cm)^-1
  double P4 = 0.99408; // kV/cm
  double P5 = 0.01172; // (kV/cm)^-P6
  double P6 = 4.20214;
  double T0 = 105.749; // K

  // Walkowiak Parameter set
  double P1W = -0.01481; // K^-1
  double P2W = 0.0075; // K^-1
  double P3W = 0.141; // (kV/cm)^-1
  double P4W = 12.4; // kV/cm
  double P5W = 1.627; // (kV/cm)^-P6
  double P6W = 0.317;
  double T0W = 90.371; // K

  // From Craig Thorne . . . currently not documented
  // smooth transition from linear at small fields to
  // icarus fit at most fields to Walkowiak at very high fields
  if ( magE < xFit ) vd = magE*uFit;
  else if ( magE < 0.619 ) {
    vd = ((P1*(temperature-T0)+1)
	  *(P3*magE*log(1+P4/magE) + P5*pow(magE, P6))
	  +P2*(temperature-T0));
  }
  else if ( magE < 0.699 ) {
    vd = 12.5*(magE - 0.619)*((P1W*(temperature-T0W)+1)
			      *(P3W*magE*log(1+P4W/magE) + P5W*pow(magE, P6W))
			      +P2W*(temperature-T0W)) +
      12.5*(0.699 - magE)*((P1*(temperature-T0)+1)
			   *(P3*magE*log(1+P4/magE) + P5*pow(magE, P6))
			   +P2*(temperature-T0));
  }
  else {
    vd = ((P1W*(temperature-T0W)+1)
	  *(P3W*magE*log(1+P4W/magE) + P5W*pow(magE, P6W))
	  +P2W*(temperature-T0W));
  }

  vd /= 10;

  return -vd*dir; // cm/us
}

// path drift_path(std::vector <double> init_pos, field * V, boundary b)
// {
//   const double dt = 1.e-4; // us
//   const int max_iter = 1e7;

//   double t = 0;
//   path trajectory = path(dt);
//   std::vector <double> pos = init_pos; // cm
  
//   std::cout << "Initial pos: " << "(" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")" << std::endl;
  
//   while ( true ) {
//     trajectory.steps.push_back(pos);
//     if ( (pos[0] < b.Xmin) or (pos[0] > b.Xmax) or
// 	 (pos[1] < b.Ymin) or (pos[1] > b.Ymax) or
// 	 (pos[2] < b.Zmin) or (pos[2] > b.Zmax) ) {
//       std::cout << "particle left boundary!" << std::endl;
//       trajectory.fate = "OOB";
//       break;
//     }
//     else if ( t >= dt*max_iter ) {
//       std::cout << "particle timed out!" << std::endl;
//       trajectory.fate = "OOT";
//       break;
//     }
//     else if ( b.is_in_boundary(pos[0], pos[1], pos[2]) ) {
//       std::cout << "particle terminated in a volume!" << std::endl;
//       trajectory.fate = "volume";
//       break;
//     }
//     else {
//       // assume T = boiling point of Ar for now
//       pos = pos + driftV(E(pos, V), 87.302)*dt;
//       t += dt;
//     }
//   }
  
//   trajectory.arrivalTime = t;

//   std::cout << "Final pos: " << "(" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")" << std::endl;
//   std::cout << "Arrival time: " << t << std::endl;
//   std::cout << "Fate: " << trajectory.fate << std::endl;
//   std::cout << std::endl;

//   return trajectory;
// }

void drift_path(std::vector <double> init_pos, field <double> * V, boundary b, path *& trajectory)
{
  const int max_iter = 1e7;

  double t = 0;
  trajectory = new path(dt);
  std::vector <double> pos = init_pos; // cm
  
  std::cout << "Initial pos: " << "(" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")" << std::endl;
  
  while ( true ) {
    // check OOB
    if ( (pos[0] < b.Xmin) or (pos[0] > b.Xmax) or
	 (pos[1] < b.Ymin) or (pos[1] > b.Ymax) or
	 (pos[2] < b.Zmin) or (pos[2] > b.Zmax) ) {
      std::cout << "particle left boundary!" << std::endl;
      trajectory -> fate = "OOB";
      break;
    }
    // check OOT
    if ( t >= dt*max_iter ) {
      std::cout << "particle timed out!" << std::endl;
      trajectory -> fate = "OOT";
      break;
    }
    // check boundary
    if ( b.is_in_conductor(pos[0], pos[1], pos[2]) ) {
      std::cout << "particle terminated in an electrode!" << std::endl;
      trajectory -> fate = "conductor";
      break;
    }
    // // check volume
    // if ( b.is_in_volume(pos[0], pos[1], pos[2]) ) {
    //   std::cout << "particle terminated in a volume!" << std::endl;
    //   trajectory -> fate = "volume";
    //   break;
    // }
    // move forward one step
    else {
      // assume T = boiling point of Ar for now
      std::vector <double> Efield = E(pos, V);
      std::vector <double> vel = driftV(Efield, Tb);
      
      trajectory -> pos.push_back(pos);
      trajectory -> vel.push_back(vel);
       
      pos = pos + vel*dt;
      t += dt;
    }
  }
  
  trajectory -> arrivalTime = t;

  std::cout << "Final pos: " << "(" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")" << std::endl;
  std::cout << "Arrival time: " << t << std::endl;
  std::cout << "Fate: " << trajectory -> fate << std::endl;
  std::cout << std::endl;
}

std::vector <double> ramo_induction(double q, path * driftPath, field <double> * weight)
{

  std::vector <double> pSeries;
  
  for ( int i = 1; i < driftPath -> pos.size(); i++ ) {
    std::vector <double> pos = driftPath -> pos[i];
    std::vector <double> vel = driftPath -> vel[i];
    std::vector <double> eField = E(pos, weight);

    double induced_current = q*dot(vel, eField); // C / us

    pSeries.push_back(induced_current);
    
    // std::cout << i*driftPath -> dt << '\t' << induced_current*1.e15 << std::endl;
  }

  return pSeries;
}

void sample(std::vector <double> depPos, double * sampleT, std::vector <double> * samplePos)
{
  // [0] - n0     - number of electrons in the deposition
  // [1] - DT     - transverse diffusion coefficient (cm^2/us)
  // [2] - DL     - longitudinal diffusion coefficient (cm^2/us)
  // [3] - z      - distance from deposit to drift plane (cm)
  // [4] - v      - drift velocity (cm/us)
  // [5] - lambda - inverse mfp (cm^-1)
  
  //  x  - t      - time from deposit (us)
  //  y  - rho    - radial param (cm)
  
  TF2 * pdf = new TF2("pdf", "([0]/(4*pi*[1]*x*sqrt(4*pi*[2]*x)))*exp(-pow([3] - [4]*x, 2)/(4*[2]*x) - [5]*[4]*x)*exp(-pow(y, 2)/(4*[1]*x))", 0, 1000, 0, 1000);

  // TF2 * pdf = new TF2("pdf", "[0]/(4*pi*[1]*sqrt(4*pi*[2]))*exp(-pow(x, 2)/(4*[2]))", 0, 1000, 0, 1000);

  // TF2 * pdf = new TF2("pdf", "[0]*[1]*[2]*[3]*[4]*[5]*x*y");
  
  pdf -> SetParameter(0, 1); // n0
  // pdf -> SetParameter(1, 16.3e3); // DT
  // pdf -> SetParameter(2, 6.2e3); // DL
  pdf -> SetParameter(1, 12.0e-6); // DT
  pdf -> SetParameter(2, 7.2e-6); // DL
  pdf -> SetParameter(3, depPos[2] - 2.5); // z
  pdf -> SetParameter(4, mag(driftV(0.5*zhat, Tb))); // v
  pdf -> SetParameter(5, 0); // lambda

  std::cout << mag(driftV(0.5*zhat, Tb)) << std::endl;
  
  double rho;
  double t;
  pdf -> GetRandom2(t, rho);
  
  double theta = 0;
  sampleT = &t;
  samplePos = new std::vector <double> {depPos[0] + rho*cos(theta),
					depPos[1] + rho*sin(theta),
					2.5};
  
  std::cout << t << '\t' << rho << std::endl;

  // return std::vector <double> {depPos[0], depPos[1], 2.5};
}
