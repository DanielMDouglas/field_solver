#include <iostream>
#include <fstream>
#include <vector>
// #include <pthread.h>
#include <thread>

#include "drift.h"

const double q = -1.602e-19; // C
const double m = 9.109e-31; // kg
const double mu = 1.e4; // cm^2 V^-1 s^-1
const double dt = 1.e-11; // s

const int max_iter = 1e4;
const double ds = 0.01; // differentiation distance

std::vector <double> Fe(std::vector <double> pos, scalarField V)
{
  double Ex = (V.interpolate(pos - ds*xhat) - V.interpolate(pos + ds*xhat))/(2*ds);
  double Ey = (V.interpolate(pos - ds*yhat) - V.interpolate(pos + ds*yhat))/(2*ds);
  double Ez = (V.interpolate(pos - ds*zhat) - V.interpolate(pos + ds*zhat))/(2*ds);

  return std::vector <double> {1.e7*q*Ex, 1.e7*q*Ey, 1.e7*q*Ez};
}

std::vector <double> Fd(std::vector <double> mom) {
  return (1.e4*q/(m*mu))*mom;
}

path drift_path(std::vector <double> init_pos, scalarField V, boundary b)
{
  double t = 0;
  path trajectory = path(dt);
  std::vector <double> pos = init_pos;
  std::vector <double> mom = {0, 0, 0};
  
  // std::cout << "Initial pos: " << "(" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")" << std::endl;
  
  while ( true ) {
    trajectory.steps.push_back(pos);
    if ( (pos[0] < -0.95) or (pos[0] > 0.95) or (pos[2] < -0.15) ) {
      // std::cout << "particle left boundary!" << std::endl;
      trajectory.fate = "OOB";
      break;
    }
    else if ( t >= dt*max_iter ) {
      // std::cout << "particle timed out!" << std::endl;
      trajectory.fate = "OOT";
      break;
    }
    else if ( b.is_in_boundary(pos[0], pos[1], pos[2]) ) {
      // std::cout << "particle terminated in a volume!" << std::endl;
      trajectory.fate = "volume";
      break;
    }
    else {
      std::vector <double> Fnet = Fe(pos, V) + Fd(mom);
      mom = mom + Fnet*dt;
      pos = pos + mom*(dt/m);
      t += dt;
    }
  }
  
  trajectory.arrivalTime = t;

  // std::cout << "Final pos: " << "(" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")" << std::endl;
  // std::cout << "Arrival time: " << t << std::endl;
  // std::cout << std::endl;

  return trajectory;
}

void drift_and_save(double xi, scalarField potential, boundary detector)
{
  path trajectory = drift_path(std::vector <double> {xi, 0.15, 1.5}, potential, detector);
  std::ostringstream stringStream;
  stringStream << "paths/drift_from_" << xi << ".dat";
  trajectory.print_to_file(stringStream.str());
}

int main()
{
  scalarField potential = scalarField("wire_cathode_300V.dat");
  boundary detector;

  int nPaths = 50;

  // std::vector <std::thread> threads {nPaths};
  std::thread * threads [nPaths];
  
  std::vector <double> xi_space = linspace(-0.9, 0.9, nPaths);
  // for ( double xi: linspace(-0.9, 0.9, nPaths) ) {
  for ( int i = 0; i < nPaths; i++ ) {
    double xi = xi_space[i];
    threads[i] = new std::thread (drift_and_save, xi, potential, detector);
  }

  for ( int i = 0; i < nPaths; i++ ) {
    threads[i] -> join();
  }
    
  return 0;
}
