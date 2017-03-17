#ifndef _POINT_HEADER_
#define _POINT_HEADER_


// Point in 3d space
class Point_t
{
	public:
		double x, y, z;

		// Constructor: pass the point xyz value
		Point_t( double a = 0.0, double b = 0.0, double c = 0.0 ) : x(a), y(b), z(c) {};
		~Point_t() {};		

		void normalize();
};


#endif
