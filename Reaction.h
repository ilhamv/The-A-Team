#ifndef _REACTION_HEADER_
#define _REACTION_HEADER_

#include <vector>  // vector
#include <memory>  // shared_ptr
#include <stack>   // stack
#include <cstring> // string

#include "Particle.h"
#include "Distribution.h"
#include "Point.h"
#include "XSec.h"


// Reaction base class
class Reaction_t
{
	private:
		std::shared_ptr<XSec_t> r_xs; // Reaction microXs in E[eV]
	
	public:
		// Constructor: Pass the microXs
	       	 Reaction_t( std::shared_ptr<XSec_t> x ) : r_xs(x) {}; // Pass the microXs
		~Reaction_t() {};

		// Get the microXs
		virtual double xs( const double E ) final { return r_xs->xs(E); };

		// Sample the reaction process on the working particle and the particle bank
		virtual void sample( Particle_t& P, std::stack< Particle_t >& Pbank ) = 0;
		
		// Check type
		virtual bool type( const std::string s ) = 0;
};


// Capture reaction
class Capture_Reaction : public Reaction_t 
{
	public:
		// Constructor: Pass the microXs
		 Capture_Reaction( std::shared_ptr<XSec_t> x ) : Reaction_t(x) {}; // Pass the microXs
		~Capture_Reaction() {};

		// Kill the working particle upon reaction sample
		void  sample( Particle_t& P, std::stack< Particle_t >& Pbank );
		
		// Check type
		bool type( const std::string s );
};


// Scatter reaction
class Scatter_Reaction : public Reaction_t 
{
	private:
		std::shared_ptr< Distribution_t<double> > scatter_dist; // Scattering angle distribution
		const double                              A;            // Nuclide mass
	public:
		// Constructor: Pass the microXs
		 Scatter_Reaction( std::shared_ptr<XSec_t> x, const std::shared_ptr< Distribution_t<double> >& D, const double a ) :
			Reaction_t(x), scatter_dist(D), A(a) {}; // Pass the microXs and scattering angle distribution
		~Scatter_Reaction() {};

		// Scatter the working particle
		void  sample( Particle_t& P, std::stack< Particle_t >& Pbank );
		
		// Check type
		bool type( const std::string s ); 
};


// Fission reaction
class Fission_Reaction : public Reaction_t 
{
	private:
		std::shared_ptr< Distribution_t<int> > nu_dist; // Fission multiplicity distribution
		IsotropicDirection_Distribution isotropic;      // Isotropic distribution for emerging neutron
    std::shared_ptr< Distribution_t<double> > fissionEnergy_dist;
		
	public:
		// Constructor: Pass the microXs
		 Fission_Reaction( std::shared_ptr<XSec_t> x, const std::shared_ptr< Distribution_t<int> >& D, const std::shared_ptr< Distribution_t<double> >& W ) :
		 	Reaction_t(x), nu_dist(D), fissionEnergy_dist(W) {}; // Pass the microXs and fission multiplicity distribution
		~Fission_Reaction() {};

		// Sample fission multiplicity, then appropriately pushing new fission particles to the bank
		// --> Reaction.cpp
		void sample( Particle_t& P, std::stack< Particle_t >& Pbank );
		
		// Check type
		bool type( const std::string s );
};


#endif
