#ifndef _REGION_HEADER_
#define _REGION_HEADER_

#include <vector>  // vector
#include <memory>  // shared_ptr
#include <stack>   // stack
#include <cstring> // string

#include "Particle.h"
#include "Surface.h"
#include "Material.h"
#include "Estimator.h"
#include "Point.h"
#include "Const.h"  // MAX_less


// Forward declaration
class Surface_t;
class Estimator_t;

// Region base class
class Region_t
{
	private:
		const std::string                                             r_name;       // Region name
		const double                                                  r_volume;     // Region volume (required if flux estimator is employed)
		const double                                                  r_importance; // Region importance (for variance reduction)
		std::vector< std::pair< std::shared_ptr< Surface_t >, int > > surfaces;     // Pairs of bounding surface and its sense
		std::shared_ptr< Material_t >                                 material;     // Contained material
		std::vector< std::shared_ptr< Estimator_t > >                 estimators;   // Estimators attached

	public:
		// Constructor: pass name, importance, and volume
     		 Region_t( const std::string n, const double imp, const double vol ) :
			 r_name(n), r_volume(vol), r_importance(imp) {};
    		~Region_t() {};

    		// Getters
		std::string name();       // name
		double      volume();     // volume
		double      importance(); // importance
		// MacroXsec of the contained material
		double      SigmaT();
		double      SigmaS();
		double      SigmaC();
		double      SigmaF();
		double      SigmaA();
		
		// Set the material
    		void setMaterial( const std::shared_ptr< Material_t >& M );
		
		// Add a bounding surface
		void addSurface ( const std::shared_ptr< Surface_t  >& S, const int sense );
		
		// Add estimator
		void addEstimator( const std::shared_ptr< Estimator_t >& E ); 

    		// Test if particle is in the region
		bool testPoint( const Point_t& p );
		
		// Move particle and score any estimators
		void moveParticle( Particle_t& P, const double dmove );
		
		// Return the closest bounding surface and the corresponding particle hit distance 
		std::pair< std::shared_ptr< Surface_t >, double > surface_intersect( const Particle_t& P );
		
		// Return particle collision distance
		virtual double collision_distance();

		// Let the Material take care of the collision sample and reaction process
		virtual void collision( Particle_t& P, std::stack< Particle_t >& Pbank );
};


// Vacuum region
// assigned for regions having no material
// overriding the collision distance and collision function
class Region_Vacuum : public Region_t
{
	public:
		// Constructor: pass name, importance, and volume
     		 Region_Vacuum( const std::string n, const double imp, const double vol ) : Region_t(n, imp, vol) {};
    		~Region_Vacuum() {};
		
		// Return sligthly less than very large number for collision distance
		// to ensure collision if no surface intersection
		double collision_distance();

		// Kill particle at collision
		void collision( Particle_t& P, std::stack< Particle_t >& Pbank );
};


#endif
