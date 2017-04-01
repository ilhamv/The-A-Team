#ifndef _SURFACE_HEADER_
#define _SURFACE_HEADER_

#include <memory>      // shared_ptr
#include <vector>      // vector
#include <cstring>     // string

#include "Const.h"     // EPSILON
#include "Particle.h"
#include "Estimator.h"
#include "Point.h"


// Forward declarations
class Region_t;
class Estimator_t;


// Surface class base
class Surface_t
{
	private:
		const std::string  s_name; // Surface name

	protected:
		std::vector< std::shared_ptr< Estimator_t > > estimators; // Estimators attached
		
		// Crossing the surface --> an epsilon kick to the working particle
		virtual void cross ( Particle_t& P ) final { P.move( EPSILON ); }
		
	public:
     		 Surface_t( const std::string n ) : s_name(n) {}; // Pass the ID#
    		~Surface_t() {};

    		// Get name
		virtual std::string name()  final { return s_name; };
		
		// Add an estimator
		virtual void addEstimator( const std::shared_ptr< Estimator_t >& E ) final 
		{ estimators.push_back( E ); }

		// Evaluate point location via the "S" equation
		virtual double eval( const Point_t& p ) = 0;

		// Get particle distance to surface
    		virtual double distance( const Particle_t& P ) = 0;
		
		// Hit surface implementation
		// (non-reflective)
		virtual void hit( Particle_t& P, const std::vector<std::shared_ptr<Region_t>>& Region );
};


////////////////////
// Plane Surfaces //
////////////////////



// Plane surface
class Plane_Surface : public Surface_t
{
  	protected:
    		const double a, b, c, d; // Parameters for S equation
  	public:
     		 Plane_Surface( const std::string n, const double pa, const double pb, const double pc, const double pd ) : 
			 Surface_t(n), a(pa), b(pb), c(pc), d(pd) {};
    		~Plane_Surface() {};

		// Employ the plane surface equation
		double eval    ( const Point_t& p );
     		double distance( const Particle_t& P );
};


// Plane surface - Reflective
// *Seems redundantly hard-coded, yet it avoids reflective flag check at every surface intersection
class Plane_Reflective : public Plane_Surface
{
	private:
		double modx, mody, modz;
	
	public:
     		 Plane_Reflective( const std::string n, const double pa, const double pb, const double pc, const double pd ) :
			Plane_Surface(n, pa, pb, pc, pd)
		{
			const double L = 2.0 / ( a*a + b*b + c*c );
			modx = L * a;
			mody = L * b;
			modz = L * c;
		};
    		~Plane_Reflective() {};

		// Reflect angle, then epsilon kick at hit
		void hit( Particle_t& P, const std::vector<std::shared_ptr<Region_t>>& Region );
};


// Sphere surface
class Sphere_Surface : public Surface_t 
{
	private:
    		const double x0, y0, z0, rad;
  	
	public:
     		 Sphere_Surface( const std::string n, const double p1, const double p2, const double p3, const double p4 ) : 
       			Surface_t(n), x0(p1), y0(p2), z0(p3), rad(p4) {};
    		~Sphere_Surface() {};

     		double eval    ( const Point_t& p );
     		double distance( const Particle_t& P );
};


// Infinite Cylinder-X surface
class CylinderX_Surface : public Surface_t
{
	private:
		const double y0, z0, rad;

	public:
		 CylinderX_Surface( const std::string n, const double p1, const double p2, const double p3 ) :
			Surface_t(n), y0(p1), z0(p2), rad(p3) {};
		~CylinderX_Surface() {};

		double eval    ( const Point_t& p );
		double distance( const Particle_t& P );
};


// Infinite Cylinder-Y surface
class CylinderY_Surface : public Surface_t
{
	private:
		const double x0, z0, rad;

	public:
		 CylinderY_Surface( const std::string n, const double p1, const double p2, const double p3 ) :
			Surface_t(n), x0(p1), z0(p2), rad(p3) {};
		~CylinderY_Surface() {};

		double eval    ( const Point_t& p );
		double distance( const Particle_t& P );
};


// Infinite Cylinder-Z surface
class CylinderZ_Surface : public Surface_t
{
	private:
		const double x0, y0, rad;

	public:
		 CylinderZ_Surface( const std::string n, const double p1, const double p2, const double p3 ) :
			Surface_t(n), x0(p1), y0(p2), rad(p3) {};
		~CylinderZ_Surface() {};

		double eval    ( const Point_t& p );
		double distance( const Particle_t& P );
};


#endif
