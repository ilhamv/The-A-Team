#include <vector> // vector
#include <memory> // shared_ptr
#include <cassert>

#include "Random.h"
#include "Nuclide.h"


// Getters
std::string Nuclide_t::name()   { return n_name; } // Name
// microXs
double Nuclide_t::sigmaS( const double E ) 
{ 
	for ( auto& r : reactions )
	{
		if ( r->type("scatter") ) { return r->xs( E ); }
	}
	// Nuclide doesn't have the reaction
	return 0.0;
}
double Nuclide_t::sigmaC( const double E ) 
{ 
	for ( auto& r : reactions )
	{
		if ( r->type("capture") ) { return r->xs( E ); }
	}
	// Nuclide doesn't have the reaction
	return 0.0;
}
double Nuclide_t::sigmaF( const double E ) 
{ 
	for ( auto& r : reactions )
	{
		if ( r->type("fission") ) { return r->xs( E ); }
	}
	// Nuclide doesn't have the reaction
	return 0.0;
}
double Nuclide_t::sigmaT( const double E )
{ 
	double sum = 0.0;

	for ( auto& r : reactions )
	{ sum += r->xs( E ); }

	return sum; 
}


// Add reaction
void Nuclide_t::addReaction( const std::shared_ptr< Reaction_t >& R ) 
{ 
	if ( R->type("scatter") ) { scatter = R; } // Attach pointer on scattering reaction
	reactions.push_back( R ); 
}


// Randomly sample a reaction type from the nuclide
std::shared_ptr< Reaction_t > Nuclide_t::reaction_sample( const double E ) 
{
	double u = sigmaT( E ) * Urand();
	double s = 0.0;
	for ( auto& r : reactions ) 
	{
		s += r->xs( E );
		if ( s > u ) { return r; }
	}
    assert( false ); // should never reach here
    return nullptr;
}


// Simulate scattering for scattering matrix MGXS
void Nuclide_t::simulate_scatter( Particle_t& P )
{
	std::stack<Particle_t> null;
	scatter->sample(P,null);
}
