#include "Shannon_Entropy_Mesh.h"

void Shannon_Entropy_Mesh::clear()
{
	for ( int i = 0; i < x_nmesh; i++ )
	{
		for ( int j = 0; j < y_nmesh; j++ )
		{
			for ( int k = 0; k < z_nmesh; k++ )
			{
				mesh.at(i).at(j).at(k) = 0.0;
			}
		}
	}
	return;
}

void Shannon_Entropy_Mesh::update( Particle_t& P, double N )
{
	int x_index = floor( ( ( P.pos().x - xmin ) * x_nmesh ) / (xmax - xmin) );
	int y_index = floor( ( ( P.pos().y - ymin ) * y_nmesh ) / (ymax - ymin) );
	int z_index = floor( ( ( P.pos().z - zmin ) * z_nmesh ) / (zmax - zmin) );

	mesh.at(x_index).at(y_index).at(z_index) += N;
	return;
}

double Shannon_Entropy_Mesh::entropy( double N )
{
	double e;
	for ( int i = 0; i < x_nmesh; i++ )
	{
		for ( int j = 0; j < y_nmesh; j++ )
		{
			for ( int k = 0; k < z_nmesh; k++ )
			{
				double p = mesh.at(i).at(j).at(k) / N;
				if ( p > 0.0 ) { e -= p * std::log(p); }
			}
		}
	}
	return e;
}