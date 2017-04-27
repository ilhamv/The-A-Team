#ifndef _SHANNON_ENTROPY_MESH_HEADER_
#define _SHANNON_ENTROPY_MESH_HEADER_

#include <vector>
#include <cmath>

#include "Particle.h"

class Shannon_Entropy_Mesh
{
	private:
		double xmin, xmax, ymin, ymax, zmin, zmax;
		int x_nmesh, y_nmesh, z_nmesh;

		std::vector < std::vector < std::vector <double> > > mesh;

	public:
		Shannon_Entropy_Mesh( double x, double X, int xn, double y, double Y, int yn, double z, double Z, int zn ) :
			xmin(x), xmax(X), x_nmesh(xn), ymin(y), ymax(Y), y_nmesh(yn), zmin(z), zmax(Z), z_nmesh(zn)
		{
			std::vector <double> v;
			for ( int k = 0; k < z_nmesh; k++ ) { v.push_back(0.0); }
			std::vector < std::vector <double> > vv;
			for (int j = 0; j < y_nmesh; j++ ) { vv.push_back(v); }
			for (int i = 0; i < x_nmesh; i++ ) { mesh.push_back(vv); }
		}
		~Shannon_Entropy_Mesh() {};

		void clear();
		void update( Particle_t& P, double N );
		double entropy( double total_n );
};

#endif