#include <cmath>      // cos, sin, sqrt

#include "Geometry.h"
#include "Particle.h"
#include "Random.h"
#include "Point.h"
#include "Const.h"    // PI2


// Getters
Point_t                   Particle_t::pos()    const { return p_pos; }    // Position
Point_t                   Particle_t::dir()    const { return p_dir; }    // Direction
bool                      Particle_t::alive()  const { return p_alive; }  // Alive/dead flag
double                    Particle_t::weight() const { return p_weight; } // Weight
double                    Particle_t::time()   const { return p_time; }   // Time
double                    Particle_t::energy() const { return p_energy; } // Energy
double                    Particle_t::speed()  const { return p_speed; }  // Speed
std::shared_ptr<Region_t> Particle_t::region() const { return p_region; } // Particle region


// Setters
void Particle_t::setDirection( const Point_t& p )                { p_dir    = p; } // Direction
void Particle_t::setWeight( const double w )                     { p_weight = w; } // Weight
void Particle_t::setRegion( const std::shared_ptr<Region_t>& R ) { p_region = R; } // Particle region
void Particle_t::setTime( const double t )                       { p_time   = t; } // Elapsed time
void Particle_t::setEnergy( const double E )                                       // Energy, and speed
{ 
	p_energy = E;
	p_speed  = std::sqrt( p_energy * 191313184.022 ) * 100.0; // cm/s
	// note, the constant above is 2.0 * ( 1.60218e-19 J/eV ) / ( 1.674929e-27 kg )
}
void Particle_t::setSpeed( const double v )                                        // Speed, and energy
{ 
	p_speed  = v; // cm/s
	p_energy = 5.22703129e-13 * v * v;
	// note, the constant above is 0.5 / ( 1.60218e-19 J/eV ) * ( 1.674929e-27 kg ) / ( 10000 m/cm )
}


// Move particle a distance dmove along its current trajectory
void Particle_t::move( const double dmove ) 
{
	p_pos.x += p_dir.x * dmove;
	p_pos.y += p_dir.y * dmove;
	p_pos.z += p_dir.z * dmove;
	
	// Advance particle time
	setTime( p_time + dmove / p_speed );
}


// Scatter particle with scattering angle mu0, with nucleus having mass A
// Sample azimuth angle azi to get new direction x
// Scattering angle mu0 is sampled in and passed by the Reaction (see Reaction.h)
// Scattering is trated in Center of mass (COM) frame
void Particle_t::scatter( const double mu0, const double A )
{
	// Particle velocity - LAB
	Point_t v_lab( p_speed * p_dir.x, p_speed * p_dir.y, p_speed * p_dir.z );
	
	// COM velocity
	const Point_t u( v_lab.x / ( 1.0 + A ), v_lab.y / ( 1.0 + A ), v_lab.z / ( 1.0 + A ) );
	
	// Particle velocity - COM
	Point_t v_c( v_lab.x - u.x, v_lab.y - u.y, v_lab.z - u.z );
	
	// Particle speed - COM
	const double speed_c = std::sqrt( v_c.x*v_c.x+ v_c.y*v_c.y+ v_c.z*v_c.z );

	// Particle initial direction - COM
	const Point_t dir_c( v_c.x / speed_c, v_c.y / speed_c, v_c.z / speed_c );
	
	// Sample azimuthal direction and prepare for scattering the direction in COM
	const double     azi = PI2 * Urand();
	const double cos_azi = std::cos(azi);
	const double sin_azi = std::sin(azi);
	const double      Ac = std::sqrt( 1.0 - mu0 * mu0 );
	Point_t dir_cNew; // To store final direction - COM

	if( dir_c.z != 1.0 )
	{
		const double       B = std::sqrt( 1.0 - dir_c.z * dir_c.z );
		const double       C = Ac / B;
		
		dir_cNew.x = dir_c.x * mu0 + ( dir_c.x * dir_c.z * cos_azi - dir_c.y * sin_azi ) * C;
		dir_cNew.y = dir_c.y * mu0 + ( dir_c.y * dir_c.z * cos_azi + dir_c.x * sin_azi ) * C;
		dir_cNew.z = dir_c.z * mu0 - cos_azi * Ac * B;
	}
	
	// If dir_c = k, interchange z and y in the scattering formula
	else
	{
		const double       B = std::sqrt( 1.0 - dir_c.y * dir_c.y );
		const double       C = Ac / B;
		
		Point_t            q; // to store new direction point
        
		q.x = dir_c.x * mu0 + ( dir_c.x * dir_c.y * cos_azi - dir_c.z * sin_azi ) * C;
		q.z = dir_c.z * mu0 + ( dir_c.z * dir_c.y * cos_azi + dir_c.x * sin_azi ) * C;
		q.y = dir_c.y * mu0 - cos_azi * Ac * B;
		
		dir_cNew.x = q.x;
		dir_cNew.y = q.y;
		dir_cNew.z = q.z;
	}

	// Final velocity - COM
	v_c.x = speed_c * dir_cNew.x;
	v_c.y = speed_c * dir_cNew.y;
	v_c.z = speed_c * dir_cNew.z;
	
	// Final velocity - LAB
	v_lab.x = v_c.x + u.x;
	v_lab.y = v_c.y + u.y;
	v_lab.z = v_c.z + u.z;

	// Final speed - LAB
	setSpeed( std::sqrt( v_lab.x*v_lab.x+ v_lab.y*v_lab.y+ v_lab.z*v_lab.z ) ); // Final energy is computed as well

	// Final direction - LAB
	p_dir.x = v_lab.x / p_speed;
	p_dir.y = v_lab.y / p_speed;
	p_dir.z = v_lab.z / p_speed;
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
