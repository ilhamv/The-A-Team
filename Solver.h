#ifndef _SOLVER_HEADER_
#define _SOLVER_HEADER_

#include "Const.h" // MAX


// Quad solver for sphere surface
// return smallest positive real root if it exists; if it does not, return very big number
// *Seems redundantly hard-coded, yet it saves up the multiplication operation with a
double sphere_quad( const double b, const double c ) 
{
  	double d = b*b - 4.0 * c;
    	
	// roots are complex, no intersection, return huge number
	// or identical roots, tangent, return huge number
  	if ( d <= 0.0 ) { return MAX; }	
  	
	else 
	{
    		double sqrtd = std::sqrt(d);

    		double r1 = 0.5 * ( -1.0 * b - sqrtd );
    		double r2 = 0.5 * ( -1.0 * b + sqrtd );

		// Negative roots return huge number
		if ( r1 < 0 ) { r1 = MAX; }
		if ( r2 < 0 ) { r2 = MAX; }

    		return std::fmin( r1, r2 );
  	}
}


// Quad solver for infinite cylinder surface
// return smallest positive real root if it exists; if it does not, return very big number
double solve_quad( const double a, const double b, const double c ) 
{
  	double d = b*b - 4.0 * a * c;
    	
	// roots are complex, no intersection, return huge number
	// or identical roots, tangent, return huge number
  	if ( d <= 0.0 ) { return MAX; }
  	
	else 
	{
    		double sqrtd = std::sqrt(d);
		double ai = 0.5 / a;

    		double r1 = ai * ( -1.0 * b - sqrtd );
    		double r2 = ai * ( -1.0 * b + sqrtd );

		// Negative roots return huge number
		if ( r1 < 0 ) { r1 = MAX; }
		if ( r2 < 0 ) { r2 = MAX; }

    		return std::fmin( r1, r2 );
  	}
}


#endif
