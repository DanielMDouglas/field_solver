#include <G4GDMLParser.hh>

#include <boundary.h>

boundary::boundary()
{
  volumes[0] = new volume(0., 100., 0., 100., 0., 2., 0.);
  volumes[1] = new volume(0., 100., 0., 100., 98., 100., 0.);
  volumes[2] = new volume(0., 2., 0., 100., 0., 100., 0.);
  volumes[3] = new volume(98., 100., 0., 100., 0., 100., 0.);
  volumes[4] = new volume(0., 100., 0., 2., 0., 100., 0.);
  volumes[5] = new volume(0., 100., 98., 100., 0., 100., 0.);

  for ( uint i = 0; i < 6; i++ ) {
    if ( volumes[i] -> Xmin < Xmin ) {
      Xmin = volumes[i] -> Xmin;
    }
    if ( volumes[i] -> Xmax > Xmax ) {
      Xmax = volumes[i] -> Xmax;
    }
    if ( volumes[i] -> Ymin < Ymin ) {
      Ymin = volumes[i] -> Ymin;
    }
    if ( volumes[i] -> Ymax > Ymax ) {
      Ymax = volumes[i] -> Ymax;
    }
    if ( volumes[i] -> Zmin < Zmin ) {
      Zmin = volumes[i] -> Zmin;
    }
    if ( volumes[i] -> Zmax > Zmax ) {
      Zmax = volumes[i] -> Zmax;
    }
  }

  std::cout << Xmin << '\t' << Xmax << '\n'
	    << Ymin << '\t' << Ymax << '\n'
	    << Zmin << '\t' << Zmax << '\n'
	    << std::endl;
}

bool boundary::is_in_boundary(double x, double y, double z)
{
  bool is_in_any = false;
  for ( uint i = 0; i < 6; i++ ) {
    if ( volumes[i] -> is_in_boundary(x, y, z) ) {
      is_in_any = true;
    }
  }
  return is_in_any;
}
