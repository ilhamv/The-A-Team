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
		double                                       n_sigmaT = 0.0; // Total microXsec
    		double                                       n_sigmaS;       // Scattering microXsec
    		double                                       n_sigmaC;       // Capture microXsec
    		double                                       n_sigmaF;       // Fission microXsec
		std::vector< std::shared_ptr< Reaction_t > > reactions;      // Reactions

	public:
		 Nuclide_t( const std::string str ) : n_name(str) {};
           	~Nuclide_t() {};

		// Getters
		std::string name(); // Name
		// microXs
		double      sigmaT();
		double      sigmaS();
		double      sigmaC();
		double      sigmaF();

		// Add reaction
		void   addReaction( const std::shared_ptr< Reaction_t >& R );
		
		// Randomly sample a reaction type from the nuclide
		std::shared_ptr<Reaction_t> reaction_sample();
};


#endif
