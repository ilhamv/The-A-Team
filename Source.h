#ifndef _SOURCE_HEADER_
#define _SOURCE_HEADER_

#include <vector>     // vector
#include <memory>     // shared_ptr

#include "Particle.h"
#include "Const.h"    // EPSILON
#include "Surface.h"
#include "Region.h"


// Particle source base class
class Source_t
{
  	public:
 		 Source_t() {};
		~Source_t() {};

		// Get the particle source
		virtual Particle_t getSource() = 0;
};


// Normal Hyperplane Source (1D x-direction use only!)
class Normal_HPlanex_Source : public Source_t
{
	private:
		// Source location and direction sense wrt. the surface
		const double posx; 
		const int    sense;
	
	public:
		// Source is slightly moved away from the surface (+EPSILON*direction)
 		 Normal_HPlanex_Source( const double p1, const int p2 ) : posx( p1 + EPSILON*p2 ), sense(p2) {};
		~Normal_HPlanex_Source() {};
		
		// Get the particle source
		Particle_t getSource();
};


// Isotropic Hyperplane Flux Source (1D x-direction use only!)
class Isotropic_HPlanex_Source : public Source_t
{
	private:
		// Source location and direction sense wrt. the surface
		const double posx;
	       	const int    sense;
	
	public:
		// Source is slightly moved away from the surface (+EPSILON*direction)
 		 Isotropic_HPlanex_Source( const double p1, const int p2 ) : posx( p1 + EPSILON*p2 ), sense(p2) {};
		~Isotropic_HPlanex_Source() {};
		
		// Get the particle source
		Particle_t getSource();
};


// Isotropic Uniform Slab Source (1D x-direction use only!)
class Isotropic_Slabx_Source : public Source_t
{
	private:
		// Direction and position distribution 
		IsotropicDirection_Distribution       isotropic; // Isotropic direction distribution
		std::shared_ptr<Uniform_Distribution> uniform;   // Uniform position distribution
	
	public:
		// Source is slightly moved away from the surface (+EPSILON*direction)
 		 Isotropic_Slabx_Source( const double p1, const double p2 )
		 { uniform = std::make_shared<Uniform_Distribution> ( p1, p2 ); }

		~Isotropic_Slabx_Source() {};

		// Get the particle source
		Particle_t getSource();
};


// Isotropic Point Source
class Isotropic_Point_Source : public Source_t
{
	private:
		Point_t                         pos;       // Position
		IsotropicDirection_Distribution isotropic; // Isotropic direction distribution

	public:
		Isotropic_Point_Source( const double p1, const double p2, const double p3 )
		{ pos.x = p1; pos.y = p2; pos.z = p3; }
		~Isotropic_Point_Source() {};

		// Get the particle source
		Particle_t getSource();
};


// DiskZ Source
class DiskZ_Source : public Source_t
{
	private:
		// Center position radius and direction sense of the source
		const double x0, y0, z0, r;
		const int    sense;

	public:
		 DiskZ_Source( const double p1, const double p2, const double p3, const double p4, const int p5 ) :
			x0(p1), y0(p2), z0(p3), r(p4), sense(p5) {};
		~DiskZ_Source() {};

		// Get the particle source with rejection sampling
		// Direct methods are way too costly (cos, sin, sqrt)
		// Acceptance probability is pretty good ~3.14/4
		Particle_t getSource();
};


// Spherical Shell Source
class Sphere_Shell_Source : public Source_t
{
	private:
		// Center position and outer radius
		const double x0, y0, z0, ro;
		// Inner radius normalized by the outer radius, then squared
		const double risq;
		IsotropicDirection_Distribution isotropic; // Isotropic direction distribution

	public:
		 Sphere_Shell_Source( const double p1, const double p2, const double p3, const double p4, const double p5 ) :
			x0(p1), y0(p2), z0(p3), risq( p4*p4/p5/p5 ), ro(p5) {};
		~Sphere_Shell_Source() {};

		// Get the particle source with rejection sampling
		// Direct methods are way too costly (acos, sin, sqrt)
		// Acceptance probability is almost half, ~3.14/6, yet it is still faster
		Particle_t getSource();
};


// User defined source
class User_Source : public Source_t
{
  	private:
    		// Position and direction distribution
		const std::shared_ptr< Distribution_t<Point_t> > dist_pos;
    		const std::shared_ptr< Distribution_t<Point_t> > dist_dir;
  	public:
     		User_Source( const std::shared_ptr< Distribution_t<Point_t> > pos, const std::shared_ptr< Distribution_t<Point_t> > dir )
			: dist_pos(pos), dist_dir(dir) {};
    		~User_Source() {};
    		Particle_t getSource();
};


// Source Bank
// collection of sources and its probability (or ratio)
// interface for every individual sources to the simulation
class Source_Bank
{
	private:
		std::vector< std::pair< std::shared_ptr<Source_t>, double > >  sources;
		double                                                         total = 0.0; // total probability (or ratio)
	
	public:
		 Source_Bank() {};
		~Source_Bank() {};
		
		void addSource( const std::shared_ptr<Source_t>& S, const double prob = 1.0 )
		{ 
			sources.push_back( std::make_pair( S, prob ) );
			total += prob;
		}
		
		// Get source
		// sources are sampled wrt to their probability
		// then, particle region is searched and set
		Particle_t getSource( const std::vector<std::shared_ptr<Region_t>>& Region )
		{
			const double xi = total * Urand();
			double s  = 0.0;
			for ( auto& So : sources ) 
			{
				// first is source, second is ratio
				s += So.second;
				if ( s > xi ) 
				{ 
					Particle_t P = So.first->getSource();
					P.searchRegion( Region );
					return P;
				}
			}
		}
};


#endif
