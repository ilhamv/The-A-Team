#include "Surface.h"
#include "Particle.h"
#include "Point.h"
#include "Const.h"  // MAX, EPSILON
#include "Solver.h" // sphere_quad, cylinder_quad



// Non-reflective surface hit implementation
void Surface_t::hit( Particle_t& P, const std::vector<std::shared_ptr<Region_t>>& Region )
{
	// Cross the surface (Epsilon kick)
	cross( P );
	// Search and set new region
	P.searchRegion( Region );
	// Score estimators 
	for ( auto e : estimators ) { e->score( P ); }
}



// Surface type methods: eval, distance and reflective
// eval          : Evaluate point location via the "S" equation
// distance      : Get particle distance to surface wth respect to its position and direction
//                 determines the mininum positive distance to intersection
//                 (returns a very large number if no intersection for ease of calculation down the line)
// reflective hit: Reflect angle corresponding to the grad of the S equation
//                 [ still for general plane only ]


///////////
// Plane //
///////////

// Plane-general
double Plane_Surface::eval( const Point_t& p )
{
	return a * p.x  +  b * p.y  +  c * p.z  - d;
}

double Plane_Surface::distance( const Particle_t& P )
{
	Point_t pos = P.pos();
	Point_t dir = P.dir();

	const double denom = a * dir.x  +  b * dir.y  +  c * dir.z;

	// Check if particle moves in a direction that is (or very close to) parallel to the surface
  	if ( std::fabs( denom ) > EPSILON ) 
	{
    		double dist = ( d - a * pos.x - b * pos.y - c * pos.z ) / denom;
    		// Check if particle moves away from the surface
		if ( dist > 0.0 ) { return dist; }
    		else { return MAX; }
  	}    	
	// It does! (Parallel)
  	else 
	{ return MAX; }
}

void Plane_Reflective::hit( Particle_t& P, const std::vector<std::shared_ptr<Region_t>>& Region  )
{
	// Reflect angle
	const double K = ( a * P.dir().x + b * P.dir().y + c * P.dir().z );
	Point_t q;
	
	q.x = P.dir().x - modx * K;
	q.y = P.dir().y - mody * K;
	q.z = P.dir().z - modz * K;

	P.setDirection(q);

	// Cross the surface (Epsilon kick)
	cross ( P );
	for ( auto e : estimators ) { e->score( P ); }
}



////////////
// Sphere //
////////////

double Sphere_Surface::eval( const Point_t& p ) 
{
  	const double x_t = p.x - x0;
  	const double y_t = p.y - y0;
  	const double z_t = p.z - z0;
	return x_t*x_t + y_t*y_t + z_t*z_t - rad_sq;
}


double Sphere_Surface::distance( const Particle_t& P ) 
{
  	Point_t p = P.pos();
  	Point_t u = P.dir();

  	// put into quadratic equation form: a*s^2 + b*s + c = 0, where a = 1
  	double b = 2.0 * ( ( p.x - x0 ) * u.x + ( p.y - y0 ) * u.y + ( p.z - z0 ) * u.z );
  	double c = eval( p );

  	return sphere_quad( b, c );
}



///////////////////////
// Infinite Cyliders //
///////////////////////

// Cylinder-X
double CylinderX_Surface::eval( const Point_t& p )
{
	const double y_t = p.y - y0;
	const double z_t = p.z - z0;
	return y_t*y_t + z_t*z_t - rad_sq;
}

double CylinderX_Surface::distance( const Particle_t& P )
{
  	Point_t p = P.pos();
  	Point_t u = P.dir();

  	// a*s^2 + b*s + c = 0
  	double a = 1.0 - u.x*u.x;
	double b = 2.0 * ( ( p.y - y0 ) * u.y + ( p.z - z0 ) * u.z );
  	double c = eval( p );

  	return solve_quad( a, b, c );
}


// Cylinder-Y
double CylinderY_Surface::eval( const Point_t& p )
{
	const double x_t = p.x - x0;
	const double z_t = p.y - z0;
	return x_t*x_t + z_t*z_t - rad_sq;
}

double CylinderY_Surface::distance( const Particle_t& P )
{
  	Point_t p = P.pos();
  	Point_t u = P.dir();

  	// a*s^2 + b*s + c = 0
  	double a = 1.0 - u.y*u.y;
	double b = 2.0 * ( ( p.x - x0 ) * u.x + ( p.z - z0 ) * u.z );
  	double c = eval( p );

  	return solve_quad( a, b, c );
}


// Cylinder-Z
double CylinderZ_Surface::eval( const Point_t& p )
{
	const double x_t = p.x - x0;
	const double y_t = p.y - y0;
	return x_t*x_t + y_t*y_t - rad_sq;
}


double CylinderZ_Surface::distance( const Particle_t& P )
{
  	Point_t p = P.pos();
  	Point_t u = P.dir();

  	// a*s^2 + b*s + c = 0
  	double a = 1.0 - u.z*u.z;
	double b = 2.0 * ( ( p.x - x0 ) * u.x + ( p.y - y0 ) * u.y );
  	double c = eval( p );

  	return solve_quad( a, b, c );
}


// Cone-X
double ConeX_Surface::eval( const Point_t& p )
{
	const double x_t = p.x - x0;
	const double y_t = p.y - y0;
	const double z_t = p.z - z0;
	return - rad_sq * x_t*x_t + y_t*y_t + z_t*z_t;
}


double ConeX_Surface::distance( const Particle_t& P )
{
  	Point_t p = P.pos();
  	Point_t u = P.dir();

  	// a*s^2 + b*s + c = 0
  	double a = 1.0 - ( rad_sq + 1.0 ) * u.x*u.x;
  	double b = 2.0 * ( - rad_sq * ( p.x - x0 ) * u.x + ( p.y - y0 ) * u.y + ( p.z - z0 ) * u.z );
  	double c = eval( p );

  	return solve_quad( a, b, c );
}

