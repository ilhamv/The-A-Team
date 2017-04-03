#include <vector>
#include <memory>
#include <cmath>
#include <cstring>

#include "Region.h"
#include "Const.h" // MAX, Max_less
#include "Point.h"


// Getters
std::string Region_t::name()       { return r_name; }       // name
double      Region_t::volume()     { return r_volume; }     // volume
double      Region_t::importance() { return r_importance; } // importance
// Get macroXsec of the contained material
double      Region_t::SigmaT()  { return material->SigmaT(); }
double      Region_t::SigmaS()  { return material->SigmaS(); }
double      Region_t::SigmaC()  { return material->SigmaC(); }
double      Region_t::SigmaF()  { return material->SigmaF(); }
double      Region_t::SigmaA()  { return material->SigmaA(); }


// Take in a pair of surface pointer and integer describing sense (must not be zero!)
// and append to vector of surfaces
void Region_t::addSurface( const std::shared_ptr< Surface_t >& S, const int sense )
{ surfaces.push_back( std::make_pair( S, sense ) ); }


// Add the material
void Region_t::setMaterial( const std::shared_ptr< Material_t >& M ) 
{ material = M; }


// Add estimator
void Region_t::addEstimator( const std::shared_ptr< Estimator_t >& E ) 
{ estimators.push_back( E ); }


// Test if point is inside the region
bool Region_t::testPoint( const Point_t& p )
{
  	// Loop over surfaces in region, if not on correct side return false
  	// if on correct side of all surfaces, particle is in the region and return true
  	for ( const auto& S : surfaces ) 
	{
    		// first = surface pointer, second = +/- 1 indicating sense
    		if ( S.first->eval( p ) * S.second < 0 ) { return false; }  
  	}
  	return true;
}


// Move particle and score any estimators
void Region_t::moveParticle( Particle_t& P, const double dmove )
{
	P.move( dmove );
	// Advance particle time
	P.setTime( P.time() + dmove / P.speed() );
	for ( const auto& e : estimators ) { e->score( P.weight(), dmove ); }
}


// Find the closest surface and travel distance for particle p to reach
std::pair< std::shared_ptr< Surface_t >, double > Region_t::surface_intersect( const Particle_t& P ) 
{
  	double dist = MAX;
	std::shared_ptr< Surface_t > S = nullptr;
  	for ( const auto& s : surfaces ) 
	{
    		double d = s.first->distance( P );
    		if ( d < dist )
		{ 
			dist = d;
			S    = s.first;
		}
  	}
	return std::make_pair( S, dist ); 
}


// Return particle collision distance
double Region_t::collision_distance()
{ return material->collision_distance_sample(); }
//
// Vacuum region: Return sligthly less than very large number for collision distance
// to ensure collision if no surface intersection
double Region_Vacuum::collision_distance()
{ return MAX_less; } // MAX_less = 0.9 MAX


// Collision
// Let the Material take care of the collision sample and reaction process
void Region_t::collision( Particle_t& P, std::stack< Particle_t >& Pbank )
{ material->collision_sample( P, Pbank ); }		
//
// Vacuum region: Kill particle at collision
void Region_Vacuum::collision( Particle_t& P, std::stack< Particle_t >& Pbank )
{ P.kill(); }		
