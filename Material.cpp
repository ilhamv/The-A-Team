#include <vector>  // vector
#include <memory>  // shared_ptr
#include <cassert>
#include <cmath>   //floor(double)
#include <iostream>

#include "Random.h"
#include "Particle.h"
#include "Reaction.h"
#include "Nuclide.h"
#include "Material.h"
#include "Const.h" // PI2


// Getters
std::string Material_t::name()   { return m_name; }   // Name
// macroXsec
double Material_t::SigmaS( const double E ) 
{ 
	double sum = 0.0;
	for ( auto& n : nuclides )
	{
		sum += n.first->sigmaS( E ) * n.second;
	}	
	return sum;
}
double Material_t::SigmaC( const double E ) 
{ 
	double sum = 0.0;
	for ( auto& n : nuclides )
	{
		sum += n.first->sigmaC( E ) * n.second;
	}	
	return sum;
}
double Material_t::SigmaF( const double E ) 
{ 
	double sum = 0.0;
	for ( auto& n : nuclides )
	{
		sum += n.first->sigmaF( E ) * n.second;
	}	
	return sum;
}
double Material_t::SigmaT( const double E ) 
{ 
	double sum = 0.0;
	for ( auto& n : nuclides )
	{
		sum += n.first->sigmaT( E ) * n.second;
	}	
	return sum;
}
double Material_t::nu( const double E )
{
	double sum = 0.0;
	double num = 0.0;
	for ( auto& n : nuclides )
	{
		sum += n.first->nu( E );
		if ( n.first->nu( E ) != 0 )
		{ num++; }
	}
	if ( num == 0 ) { return 0; }
	return sum / num;
}
double Material_t::nuSigmaF( const double E ) 
{ 
	double sum = 0.0;
	for ( auto& n : nuclides )
	{
		sum += n.first->nusigmaF( E ) * n.second;
	}	
	return sum;
}


// Add a nuclide
// the supplied variable are the nuclide and its nuclide density
void Material_t::addNuclide( const std::shared_ptr< Nuclide_t >& Nuclide, double N ) 
{ nuclides.push_back( std::make_pair( Nuclide, N ) ); }


// Sample collision distance
double Material_t::collision_distance_sample( const double E )
{ return - std::log( Urand() ) / SigmaT( E ); }


// Sample collided nuclide
std::shared_ptr< Nuclide_t > Material_t::nuclide_sample( const double E )
{
	double u = SigmaT( E ) * Urand();
	double s = 0.0;
	for ( auto& n : nuclides ) 
	{
		// first is pointer to nuclide, second is nuclide density
		s += n.first->sigmaT( E ) * n.second;
		if ( s > u ) { return n.first; }
	}
    assert( false ); // should never reach here
    return nullptr;
}

// Bank fission neutrons for k-eigenvalue simulations
void Material_t::bank_fission_neutrons ( Particle_t& P, double K, std::stack< Particle_t>& Fbank, std::shared_ptr <Shannon_Entropy_Mesh> shannon_mesh )
{
	//std::cout << "k: " << K << std::endl;
	//std::cout << "SigmaT: " << SigmaT( P.energy() ) << std::endl;
	double bank_nu = floor( ( P.weight() / K ) * nu ( P.energy() ) * ( SigmaF( P.energy() ) / SigmaT( P.energy() ) ) + Urand() );

	//std::cout << "Particle weight: " << P.weight() << std::endl;
	//std::cout << "Number of neutrons to add to fission bank: " << bank_nu << std::endl;

	if( bank_nu != 0 )
		shannon_mesh->update( P, bank_nu );

	//std::cout << "mesh updated" << std::endl;

	for ( int n = 0; n < bank_nu; n++ )
	{
		double r = Urand();
		double s = 0.0;
		for( auto& n : nuclides )
		{
			s += n.first->nusigmaF( P.energy() ) / nuSigmaF( P.energy() );
			if ( r < s )
			{
				// Sample polar cosine and azimuthal angle uniformly
				double mu  = 2.0 * Urand() - 1.0;
				double azi = PI2 * Urand();

				// Convert to Cartesian coordinates
				double c = std::sqrt( 1.0 - mu * mu );
				Point_t d;
				d.y = std::cos( azi ) * c;
				d.z = std::sin( azi ) * c;
				d.x = mu;

				Fbank.emplace( P.pos(), d, n.first->Chi( P.energy() ), P.time(), P.weight() );
				break;
			}
		}
	}
	return;
}

// Sample entire collision (nuclide, then nuclide reaction)
// Then, process the reaction on the Particle 
void Material_t::collision_sample( Particle_t& P, bool eigenvalue, double K, std::stack< Particle_t >& Pbank, std::stack< Particle_t >& Fbank, std::shared_ptr <Shannon_Entropy_Mesh> shannon_mesh ) 
{

	//std::cout << "about to bank" << std::endl;
	//First bank fission neutrons if running a k-eigenvalue simulation
	if ( eigenvalue ) { bank_fission_neutrons( P, K, Fbank, shannon_mesh ); }

	//std::cout << "fission neutrons banked" << std::endl;

	// Then sample nuclide
	std::shared_ptr< Nuclide_t >  N = nuclide_sample( P.energy() );

	//std::cout << "nuclide sampled" << std::endl;

	// Now get the reaction
	std::shared_ptr< Reaction_t > R = N->reaction_sample( P.energy() );
	
	// Finally process the reaction on the Particle
	R->sample( P, Pbank );

	return;
}
		

// Simulate scattering for scattering matrix MGXS
void Material_t::simulate_scatter( Particle_t& P )
{
	// Sample the scattering nuclide
	double u = SigmaS( P.energy() ) * Urand();
	double s = 0.0;
	for ( auto& n : nuclides ) 
	{
		// first is pointer to nuclide, second is nuclide density
		s += n.first->sigmaS( P.energy() ) * n.second;
		if ( s > u ) 
		{ 
			// Simulate the scatter
			return n.first->simulate_scatter( P );
	       	}
	}
	// There is no scattering nuclide
}
