#include <cmath> // sqrt

#include "Distribution.h"
#include "Random.h"
#include "Particle.h"
#include "Source.h"
#include "Point.h"


Particle_t Normal_HPlanex_Source::getSource()
{
	Point_t p_pos( posx, 0.0, 0.0 );
	Point_t p_dir( sense, 0.0, 0.0 );
	Particle_t P( p_pos, p_dir );
	return P;
}


Particle_t Isotropic_HPlanex_Source::getSource()
{
	const double cos_t = sense * std::sqrt( Urand() );
	const double sin_t = std::sqrt( 1.0 - cos_t * cos_t );
	
	Point_t p_pos( posx, 0.0, 0.0 );
	Point_t p_dir( cos_t, 0.0, sin_t );
	Particle_t P( p_pos, p_dir );
	return P;
}


Particle_t Isotropic_Slabx_Source::getSource()
{
	Point_t p_pos( uniform->sample(), 0.0, 0.0 );
	Particle_t P( p_pos, isotropic.sample() );
	return P;
}


Particle_t Isotropic_Point_Source::getSource()
{
	Particle_t P( pos, isotropic.sample() );
	return P;
}


Particle_t DiskZ_Source::getSource()
{
	double x,y;
	// Rejection sampling: square --> circle
	do
	{
		x = 2.0 * Urand() - 1.0;
		y = 2.0 * Urand() - 1.0;
	}
	while ( x*x + y*y > 1.0 );

	x = x * r;
	y = y * r;
	Point_t pos( x0 + x, y0 + y, z0 );
	Point_t dir( 0.0, 0.0, sense );
	Particle_t P( pos, dir );
	return P;
}


Particle_t Sphere_Shell_Source::getSource()
{
	double x,y,z;
	// Rejection sampling: Cube --> Spherical shell
	do
	{
		x = 2.0 * Urand() - 1.0;
		y = 2.0 * Urand() - 1.0;
		z = 2.0 * Urand() - 1.0;
	}
	while ( ( x*x + y*y + z*z > 1.0 ) || ( x*x + y*y + z*z < risq ) );

	x = x * ro;
	y = y * ro;
	z = z * ro;

	Point_t pos( x0 + x, y0 + y, z + z0 );
	Point_t dir( isotropic.sample() );
	Particle_t P( pos, dir );
	return P;
}


Particle_t User_Source::getSource() 
{
  	Point_t p = dist_dir->sample();
	p.normalize();
	Particle_t P ( dist_pos->sample(), p );
  	return P;
}
