#ifndef _XSEC_HEADER_
#define _XSEC_HEADER_


#include <cmath> // sqrt
#include <vector>
#include "Solver.h"


class XSec_t
{
	public:
     		 XSec_t() {};
    		~XSec_t() {};
		
		// Get cross section in energy
		virtual double xs( const double E ) const = 0;
};


// Constant cross section
class Constant_XSec : public XSec_t
{
	const double val;

	public:
		 Constant_XSec( const double x ) : val(x) {};
		~Constant_XSec() {};

		double xs( const double E ) const { return val; }
};


// 1/v cross section
class OverV_XSec : public XSec_t
{
	const double a, b; // a + b / sqrt(E)
	                   // E is in eV

	public:
		 OverV_XSec( const double p1, const double p2 ) : a(p1), b(p2) {};
		~OverV_XSec() {};

		double xs( const double E ) const { return a + b / std::sqrt(E); }
};

class lookup_XSec : public XSec_t
{
	std::vector<double> Edata, XSdata;

	public:
		 lookup_XSec( std::vector<double> p1, std::vector<double> p2 ) : Edata(p1), XSdata(p2) {};
		~lookup_XSec() {};

		double xs( const double E ) const;
        
};

#endif
