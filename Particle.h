#ifndef _PARTICLE_HEADER_
#define _PARTICLE_HEADER_

#include <memory> // shared_ptr
#include <vector> // vector

#include "Point.h"


// Forward declarations
class Region_t;


// Particle base class
class Particle_t
{
	private:
		Point_t                   p_pos;          // Position
		Point_t                   p_dir;          // Direction
		bool                      p_alive = true; // Alive/dead flag
		double                    p_weight = 1.0; // Weight (for variance reduction)
		std::shared_ptr<Region_t> p_region;       // Particle region

  	public:
		// Constructor: Pass position and direction and initialize particle
		Particle_t( const Point_t& p1, const Point_t& p2 ) : p_pos(p1), p_dir(p2) {};
		~Particle_t() {};

    		// Getters
		Point_t                   pos()    const; // Position
    		Point_t                   dir()    const; // Direction 
		bool                      alive()  const; // Alive/dead flag
		double                    weight() const; // Weight
		std::shared_ptr<Region_t> region() const; // Particle region

		// Setters
		void setDirection( const Point_t& p );                   // Direction
		void setWeight   ( const double w );                     // Weight
		void setRegion   ( const std::shared_ptr<Region_t>& R ); // Particle region

		// Kill: live/die flag --> False, weight = 0
		void kill();		                                                   
		
		// Move particle a distance dmove along its current trajectory
    		void move( const double dmove );                                           
		
		// Scatter with scattering angle mu0
		void scatter( const double mu0 );                                          
		
		// Search and set particle region
		void searchRegion( const std::vector<std::shared_ptr<Region_t>>& Region ); 
};


#endif
