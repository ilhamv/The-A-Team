#ifndef _NUCLIDE_HEADER_
#define _NUCLIDE_HEADER_

#include <memory>  // shared_ptr
#include <vector>  // vector
#include <cstring> // strcmp

#include "Reaction.h"


// Nuclide base class
class Nuclide_t
{
	private:
    		std::string                                  n_name;         // Name
		const double                                 n_A;            // Nuclide mass
		std::vector< std::shared_ptr< Reaction_t > > reactions;      // Reactions
		std::shared_ptr< Reaction_t >                    scatter;        // Scattering reaction

	public:
		 Nuclide_t( const std::string str, const double a ) : n_name(str), n_A(a) {};
           	~Nuclide_t() {};

		// Getters
		std::string name(); // Name
		// microXs
		double      sigmaT( const double E );
		double      sigmaS( const double E );
		double      sigmaC( const double E );
		double      sigmaF( const double E );
		double      sigmaScont( const double E );
		double      sigmaCcont( const double E );
		double      sigmaFcont( const double E );

		// Add reaction
		void   addReaction( const std::shared_ptr< Reaction_t >& R );
		
		// Randomly sample a reaction type from the nuclide
		std::shared_ptr<Reaction_t> reaction_sample( const double E );
		
		// Simulate scattering for scattering matrix MGXS
		void simulate_scatter( Particle_t& P );
};


#endif
