#ifndef _XSEC_HEADER_
#define _XSEC_HEADER_


#include <cmath> // sqrt


class XSec_t
{
	public:
     		 Xsec_t() {};
    		~Xsec_t() {};
		
		// Get cross section in energy
		virtual double xs( const double E ) const = 0;
};


// Constant cross section
class Constant_XSec
{
	const double val;

	public:
		 Constant_XSec( const double x ) : val(x) {};
		~Constant_XSec() {};

		double xs( const double E ) const { return val; }
}


// 1/v cross section
class OverV_XSec
{
	const double a, b; // a + b / sqrt(E)
	                   // E is in eV

	public:
		 OverV_XSec( const double p1, const double p2 ) : a(p1), b(p2) {};
		~OverV_XSec() {};

		double xs( const double E ) const { return a + b / std::sqrt(E); }
}


#endif
