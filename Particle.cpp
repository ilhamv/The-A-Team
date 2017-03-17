#include <cmath>      // cos, sin, sqrt

#include "Region.h"
#include "Particle.h"
#include "Random.h"
#include "Point.h"
#include "Const.h"    // PI2


// Getters
Point_t                   Particle_t::pos()    const { return p_pos; }    // Position
Point_t                   Particle_t::dir()    const { return p_dir; }    // Direction
bool                      Particle_t::alive()  const { return p_alive; }  // Alive/dead flag
double                    Particle_t::weight() const { return p_weight; } // Weight
std::shared_ptr<Region_t> Particle_t::region() const { return p_region; } // Particle region


// Setters
void Particle_t::setDirection( const Point_t& p )                { p_dir = p; }    // Direction
void Particle_t::setWeight( const double w )                     { p_weight = w; } // Weight
void Particle_t::setRegion( const std::shared_ptr<Region_t>& R ) { p_region = R; } // Particle region


// Move particle a distance dmove along its current trajectory
void Particle_t::move( const double dmove ) 
{
	p_pos.x += p_dir.x * dmove;
	p_pos.y += p_dir.y * dmove;
	p_pos.z += p_dir.z * dmove;
}


// Scatter particle with scattering angle mu0
// Sample azimuth angle azi to get new direction x
// Scattering angle mu0 is sampled in and passed by the Reaction (see Reaction.h)
void Particle_t::scatter( const double mu0 )
{
	const double     azi = PI2 * Urand();
	const double cos_azi = std::cos(azi);
	const double sin_azi = std::sin(azi);
	const double       A = std::sqrt( 1.0 - mu0 * mu0 );
	
	if( p_dir.z != 1.0 )
	{
		const double       B = std::sqrt( 1.0 - p_dir.z * p_dir.z );
		const double       C = A / B;
		
		Point_t            q; // to store new direction point
        
		q.x = p_dir.x * mu0 + ( p_dir.x * p_dir.z * cos_azi - p_dir.y * sin_azi ) * C;
		q.y = p_dir.y * mu0 + ( p_dir.y * p_dir.z * cos_azi + p_dir.x * sin_azi ) * C;
		q.z = p_dir.z * mu0 - cos_azi * A * B;
		
		p_dir.x = q.x;
		p_dir.y = q.y;
		p_dir.z = q.z;
	}
	
	// If Omega = k, interchange z and y in the scattering formula
	else
	{
		const double       B = std::sqrt( 1.0 - p_dir.y * p_dir.y );
		const double       C = A / B;
		
		Point_t            q; // to store new direction point
        
		q.x = p_dir.x * mu0 + ( p_dir.x * p_dir.y * cos_azi - p_dir.z * sin_azi ) * C;
		q.z = p_dir.z * mu0 + ( p_dir.z * p_dir.y * cos_azi + p_dir.x * sin_azi ) * C;
		q.y = p_dir.y * mu0 - cos_azi * A * B;
		
		p_dir.x = q.x;
		p_dir.y = q.y;
		p_dir.z = q.z;
	}
}


// Kill particle
void Particle_t::kill()
{ 
	p_alive  = false; 
	p_weight = 0.0;
}


// Search and set particle region
void Particle_t::searchRegion( const std::vector<std::shared_ptr<Region_t>>& Region )
{
	for( const auto& R : Region )
	{
		// check if particle is in the current region R
		if ( R->testPoint( p_pos ) )
		{
			p_region = R;
			return;
		}
	}
}
