#include <vector>  // vector
#include <memory>  // shared_ptr

#include "Random.h"
#include "Particle.h"
#include "Reaction.h"
#include "Nuclide.h"
#include "Material.h"


// Getters
std::string Material_t::name()   { return m_name; }   // Name
// macroXsec
double      Material_t::SigmaT() { return m_SigmaT; }
double      Material_t::SigmaS() { return m_SigmaS; }
double      Material_t::SigmaC() { return m_SigmaC; }
double      Material_t::SigmaF() { return m_SigmaF; }
double      Material_t::SigmaA() { return m_SigmaA; }


// Add a pair of nuclide and its total macroXs
// the supplied variable are the nuclide and its nuclide density
void Material_t::addNuclide( const std::shared_ptr< Nuclide_t >& Nuclide, double N ) 
{ 
	double newSigmaT =  Nuclide->sigmaT() * N; // Nuclide total macroXs
	m_SigmaT += newSigmaT;                // Accumulate material total macroXs
	m_SigmaS += Nuclide->sigmaS() * N; // Accumulate material scatter macroXs
	m_SigmaC += Nuclide->sigmaC() * N; // Accumulate material capture macroXs
	m_SigmaF += Nuclide->sigmaF() * N; // Accumulate material fission macroXs
	m_SigmaA  = m_SigmaF + m_SigmaC;          // Update absorption macroXs
	lambda  = -1.0 / m_SigmaT;            // Update lambda
	nuclides.push_back( std::make_pair( Nuclide, newSigmaT ) );
}


// Sample collision distance
double Material_t::collision_distance_sample()
{ return lambda * log( Urand() ); }


// Sample collided nuclide
std::shared_ptr< Nuclide_t > Material_t::nuclide_sample()
{
	double u = m_SigmaT * Urand();
	double s = 0.0;
	for ( auto& n : nuclides ) 
	{
		// first is pointer to nuclide, second is nuclide total macros XS
		s += n.second;
		if ( s > u ) { return n.first; }
	}
}


// Sample entire collision (nuclide, then nuclide reaction)
// Then, process the reaction on the Particle 
void Material_t::collision_sample( Particle_t& P, std::stack<Particle_t>& Pbank ) 
{
	// First sample nuclide
	std::shared_ptr< Nuclide_t >  N = nuclide_sample();

	// Now get the reaction
	std::shared_ptr< Reaction_t > R = N->reaction_sample();
	
	// Finally process the reaction on the Particle
	R->sample( P, Pbank );
}
