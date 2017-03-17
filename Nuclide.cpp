#include <vector> // vector
#include <memory> // shared_ptr

#include "Random.h"
#include "Nuclide.h"


// Getters
std::string Nuclide_t::name()   { return n_name; } // Name
// microXs
double      Nuclide_t::sigmaT() { return n_sigmaT; }
double      Nuclide_t::sigmaS() { return n_sigmaS; }
double      Nuclide_t::sigmaC() { return n_sigmaC; }
double      Nuclide_t::sigmaF() { return n_sigmaF; }


// Add reaction
void Nuclide_t::addReaction( const std::shared_ptr< Reaction_t >& R ) 
{ 
	// Accumulate total microXsec
	n_sigmaT += R->xs();
	if      ( R->type("scatter") ) { n_sigmaS = R->xs(); }
	else if ( R->type("capture") ) { n_sigmaC = R->xs(); }
	else if ( R->type("fission") ) { n_sigmaF = R->xs(); }
	reactions.push_back( R );
}


// Randomly sample a reaction type from the nuclide
std::shared_ptr< Reaction_t > Nuclide_t::reaction_sample() 
{
	double u = n_sigmaT * Urand();
	double s = 0.0;
	for ( auto r : reactions ) 
	{
		s += r->xs();
		if ( s > u ) { return r; }
	}
}
