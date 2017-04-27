#include <cmath>

#include "Point.h"

Point_t& Point_t::normalize() 
{
  	double norm = 1.0 / std::sqrt( x*x + y*y + z*z );
  	x *= norm; y *= norm; z *= norm;
	return *this;
}
