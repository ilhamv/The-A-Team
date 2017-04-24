#ifndef _DISTRIBUTION_HEADER_
#define _DISTRIBUTION_HEADER_

#include <vector>   // vector
#include <cmath>    // erf, floor
#include <string>
#include <memory>
#include <cassert>
#include <fstream>

#include "Const.h" // PI
#include "Random.h" // Urand
#include "Point.h"
#include "Particle.h"

// Sampling Distributions base class
template< class T >
class Distribution_t
{
	private:
		const std::string d_name;

	public:
     		 Distribution_t( const std::string label = "" ) : d_name(label) {};
    		~Distribution_t() {};
		
		virtual std::string name() final { return d_name; };
		
		// Get sample value
		virtual T sample() = 0;
    
        // This function is needed to pass the particle energy to sample fission
        virtual double sampleFis( const Particle_t &P ) = 0;

};


// Delta distribution
template <class T>
class Delta_Distribution : public Distribution_t<T> 
{
  	private:
    		T result;

  	public:
     		 Delta_Distribution( T val, const std::string label = "" ) : Distribution_t<T>(label), result(val) {};
    		~Delta_Distribution() {};

    		T sample() { return result; }
    
    //dummy function. it is needed to pass the particle energy to sample fission
    double sampleFis( const Particle_t &P ){return 1.0; };
};


// Discrete distribution base class
template <class T>
class Discrete_Distribution : public Distribution_t<T> 
{
	protected:
		std::vector< std::pair< T, double > > cdf;

	public:
		 Discrete_Distribution( const std::vector< std::pair< T, double > >& data, const std::string label = "" );
		~Discrete_Distribution() {};
		T sample();
    
    //dummy function. it is needed to pass the particle energy to sample fission
    double sampleFis( const Particle_t &P ){return 1.0; };
};


template < class T >
Discrete_Distribution<T>::Discrete_Distribution( const std::vector< std::pair< T, double > >& data, const std::string label /*= ""*/ ):
       Distribution_t<T>(label)	
{
	// convert pmf to cdf
	double c = 0.0;
	for ( auto& d : data ) 
	{
		// first is pointer to data type T and second is pmf input
		cdf.push_back( std::make_pair( d.first, d.second + c ) );
		c += d.second;
	}
    
}


template < class T >
T Discrete_Distribution<T>::sample() 
{
	double   r = Urand() * cdf.back().second;
	for ( auto& c : cdf ) 
	{
		// first is pointer to data type T and second is cdf
		if ( r < c.second ) { return c.first; };
	}
    assert( false ); // should never reach here
    
}


// Uniform distribution in [a,b]
class Uniform_Distribution : public Distribution_t<double>
{
  	private:
    		const double a, b, range;

	public:
     		 Uniform_Distribution( const double p1, const double p2, const std::string label = "" ) : a( p1 ), b( p2 ), range( p2 - p1 ), Distribution_t(label) {};
    		~Uniform_Distribution() {};
    		double sample();
    
            //dummy function. it is needed to pass the particle energy to sample fission
            double sampleFis( const Particle_t &P ){return 1.0; };
};


// Linear distribution
class Linear_Distribution : public Distribution_t<double> 
{
  	private:
    		double a, b, fa, fb;
  	
	public:
    		Linear_Distribution( double x1, double x2, double y1, double y2, const std::string label = "" )  
      			: Distribution_t(label), a(x1), b(x2), fa(y1), fb(y2) {};
    ~Linear_Distribution() {};
    double sample();
    
            //dummy function. it is needed to pass the particle energy to sample fission
            double sampleFis( const Particle_t &P ){return 1.0; };
};

/*
// Watt distribution
class Watt_Distribution : public Distribution_t <double>
{
private:
    std::vector<double> evChi;   // Energy grid
    std::vector<double> probChi; // PDF
    std::vector<double> cdfChi;  // CDF
    
public:
    Watt_Distribution( const std::string nameFile, const std::string label = "" ): Distribution_t(label) // Construct cdf from watt spectrum text file
    {
        std::ifstream inputFile(nameFile);
        std::string line;
        while(getline(inputFile, line))
        {
            if (!line.length() )
                continue;
            double x = 0.0, y = 0.0;
            sscanf(line.c_str(), "%lf %lf", &x, &y);
            evChi.push_back(x);
            probChi.push_back(y);
        }
        
        // Create the cdf vector
        double cdfNow = probChi[0];
        cdfChi.push_back(cdfNow);
        for(int a = 0 ; a<probChi.size()-1 ; a++)
        {
            cdfNow += 0.5 * ( evChi[a+1]-evChi[a] ) * ( probChi[a+1]+probChi[a] ) ;
            cdfChi.push_back(cdfNow);
        }
        
        // Normalize the distribution
        for(int b = 0 ; b<cdfChi.size() ; b++){
            cdfChi[b] = cdfChi[b] / cdfNow;
        }
    }
    ~Watt_Distribution() {};
    double sample();
};
 */

// Watt distribution - the hyperbolic sin model
class Watt_Distribution : public Distribution_t <double>
{
private:
    
    std::string nameNuc;
    
    double a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12;
    double b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12;
    double g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12;
    
    //default is U235 ( neutron energy 0 MeV - 1 MeV )
    double a = 0.988; //[MeV]
    double b = 2.249; //[MeV-1]
    double g = std::sqrt( pow (1.0+(0.988*2.249/8.0) ,2 ) - 1.0) + (1.0+(0.988*2.249/8.0) ) ;
    
public:
    Watt_Distribution( const std::string nuc, const std::string label = "" ): nameNuc(nuc), Distribution_t(label) {
    
        //Th232 ( 0 MeV - 0.025 MeV )
        a1 = 1.0888;
        b1 = 1.6871; 
        g1 = std::sqrt( pow (1.0+(1.0888*1.6871/8.0) ,2 ) - 1.0) + (1.0+(1.0888*1.6871/8.0) ) ;
        
        //Th232 ( 0.025 MeV - 1 MeV )
        a2 = 1.1096;
        b2 = 1.6316;
        g2 = std::sqrt( pow (1.0+(1.1096*1.6316/8.0) ,2 ) - 1.0) + (1.0+(1.1096*1.6316/8.0) ) ;
        
        //Th232 ( 1 MeV - 14 MeV    )
        a3 = 1.1700;
        b3 = 1.4610;
        g3 = std::sqrt( pow (1.0+(1.1700*1.4610/8.0) ,2 ) - 1.0) + (1.0+(1.1700*1.4610/8.0) ) ;
        
        //U233  ( 0 MeV - 0.025 MeV &&  0.025 MeV - 1 MeV )
        a4 = 0.977;
        b4 = 2.546;
        g4 = std::sqrt( pow (1.0+(0.977*2.546/8.0) ,2 ) - 1.0) + (1.0+(0.977*2.546/8.0) ) ;
        
        //U233  ( 1 MeV - 14 MeV )
        a5 = 1.0036;
        b5 = 2.6377;
        g5 = std::sqrt( pow (1.0+(1.0036*2.6377/8.0) ,2 ) - 1.0) + (1.0+(1.0036*2.6377/8.0) ) ;
        
        //U235  ( 1 MeV - 14 MeV )
        a7 = 1.028;
        b7 = 2.084;
        g7 = std::sqrt( pow (1.0+(1.028*2.084/8.0) ,2 ) - 1.0) + (1.0+(1.028*2.084/8.0) ) ;
        
        //U238  ( 0 MeV - 0.025 MeV )
        a8 = 0.88111;
        b8 = 3.4005;
        g8 = std::sqrt( pow (1.0+(0.88111*3.4005/8.0) ,2 ) - 1.0) + (1.0+(0.88111*3.4005/8.0) ) ;
        
        //U238  ( 0.025 MeV - 1 MeV )
        a9 = 0.89506;
        b9 = 3.2953;
        g9 = std::sqrt( pow (1.0+(0.89506*3.2953/8.0) ,2 ) - 1.0) + (1.0+(0.89506*3.2953/8.0) ) ;
        
        //U238  ( 1 MeV - 14 MeV    )
        a10 = 0.96534;
        b10 = 2.8330;
        g10 = std::sqrt( pow (1.0+(0.96534*2.8330/8.0) ,2 ) - 1.0) + (1.0+(0.96534*2.8330/8.0) ) ;
        
        //Pu239 ( 0 MeV - 0.025 MeV &&  0.025 MeV - 1 MeV )
        a11 = 0.966;
        b11 = 2.842;
        g11 = std::sqrt( pow (1.0+(0.966*2.842/8.0) ,2 ) - 1.0) + (1.0+(0.966*2.842/8.0) ) ;
        
        //Pu239 ( 1 MeV - 14 MeV )
        a12 = 1.055;
        b12 = 2.383; 
        g12 = std::sqrt( pow (1.0+(1.055*2.383/8.0) ,2 ) - 1.0) + (1.0+(1.055*2.383/8.0) ) ;
    
    };
    ~Watt_Distribution() {};
    
    //for the watt distribution, the sample() function is the dummy one
    double sample(){return 1.0; };
    
    //sample fission neutron energy distribution
    double sampleFis( const Particle_t &P );
};
 

// Cubic distribution [rejection sampling]
class Cubic_Distribution : public Distribution_t<double>
{
  	private:
    		const double a, b, c3, c2, c1, c0, fmax;
  		double f( const double x )
		{
			const double x_sq = x*x;
			return c3 * x_sq*x + c2 * x_sq + c1 * x + c0;
		}	

	public:
    		Cubic_Distribution( const double p1, const double p2, const double p3, double const p4,
			       const double p5, const double p6, const double p7, const std::string label = "" )  
      			: Distribution_t(label), a(p1), b(p2), c3(p3), c2(p4), c1(p5), c0(p6), fmax(p7) {};
   		~Cubic_Distribution() {};
   		double sample();
    
    //dummy function. it is needed to pass the particle energy to sample fission
    double sampleFis( const Particle_t &P ){return 1.0; };
};


// Rayleigh scattering distribution
class RayleighScatter_Distribution : public Distribution_t<double>
{
  	public:
     		 RayleighScatter_Distribution( const std::string label = "" ) : Distribution_t(label) {};
    		~RayleighScatter_Distribution() {};
    		double sample();
    
    //dummy function. it is needed to pass the particle energy to sample fission
    double sampleFis( const Particle_t &P ){return 1.0; };
};


// Normal distribution with mean and standard deviation sigma
class Normal_Distribution : public Distribution_t<double>
{
  	private:
    		const double twopi = 2.0 * PI;
    		const double       mean, sigma;
  	
	public:
     		 Normal_Distribution( const double p1, const double p2, const std::string label = "" ) : Distribution_t(label), mean(p1), sigma(p2) {};
    		~Normal_Distribution() {};
	    	double sample();
    
    //dummy function. it is needed to pass the particle energy to sample fission
    double sampleFis( const Particle_t &P ){return 1.0; };
};


// Henyey-Green scattering distribution with parameter g
class HGScatter_Distribution : public Distribution_t<double>
{
  	private:
    		const double g;
 		double A;
    		double B;
    		double C;
    		double D;
    		double E;

  	public:
     		HGScatter_Distribution( const double p, const std::string label = "" ) : g(p), Distribution_t(label)
		{
 			A = ( 1.0 + g*g ) / ( 2.0*g );
    			B = 2.0 * g / ( 1.0 - g*g );
    			C = 1.0 / sqrt( 2.0 * g + g*g + 1.0 );
			D = 1.0 / ( 2.0*g*B*B );
			E = C/B;
		};
    		~HGScatter_Distribution() {};
    		double sample();
    
    //dummy function. it is needed to pass the particle energy to sample fission
    double sampleFis( const Particle_t &P ){return 1.0; };

};


// Isotropic scattering distribution
class IsotropicScatter_Distribution : public Distribution_t<double>
{
	public:
		 IsotropicScatter_Distribution( const std::string label = "" ) : Distribution_t(label) {};
		~IsotropicScatter_Distribution() {};

    double sample(){ return 2.0 * Urand() - 1.0; }
    
    //dummy function. it is needed to pass the particle energy to sample fission
    double sampleFis( const Particle_t &P ){return 1.0; };
};


// Linear scattering distribution
// using linear decomposition
class LinearScatter_Distribution : public Distribution_t<double>
{
	private:
		const double prob; // Probability for the first pdf

	public:
		 LinearScatter_Distribution( const double mubar, const std::string label = "" ) : prob( 3.0 * mubar ), Distribution_t(label) {};
		~LinearScatter_Distribution() {};

		double sample();
    
        //dummy function. it is needed to pass the particle energy to sample fission
        double sampleFis( const Particle_t &P ){return 1.0; };
};


// Isotropic direction distribution
class IsotropicDirection_Distribution : public Distribution_t<Point_t>
{
public:
	IsotropicDirection_Distribution( const std::string label = "" ) : Distribution_t(label) {};
    ~IsotropicDirection_Distribution() {};
    Point_t sample();
    
    //dummy function. it is needed to pass the particle energy to sample fission
    double sampleFis( const Particle_t &P ) {return 1.0; };
    
};


// Average multiplicity distribution
class Average_Multiplicity_Distribution : public Distribution_t<int>
{
	private:
		const double nubar;

	public:
		 Average_Multiplicity_Distribution( const double p1, const std::string label = "" ) : nubar(p1), Distribution_t(label) {};
		~Average_Multiplicity_Distribution() {};
		int sample() { return std::floor( nubar + Urand() ); }
    
        //dummy function. it is needed to pass the particle energy to sample fission
        double sampleFis( const Particle_t &P ){return 1.0; };
};

// Terrel multiplicity distribution
class Terrel_Multiplicity_Distribution : public Discrete_Distribution<int>
{
	private:
		const double nubar, gamma, b, n;

	public:
		Terrel_Multiplicity_Distribution( const double p1, const double p2, const double p3, const int p4, const std::vector< std::pair< int, double > >& data, const std::string label = "" ) :
			nubar(p1), gamma(p2), b(p3), n(p4+1), Discrete_Distribution(data,label)
		{
			for ( int nu = 0; nu < n; nu++ )
			{
				double prob = 0.5 * ( erf( ( nu - nubar + 0.5 + b ) / ( gamma * sqrt(2.0) ) ) + 1.0 );
				cdf.push_back( std::make_pair( nu, prob ) );
			}
		}		
		~Terrel_Multiplicity_Distribution() {};
};


// Independent 3point distribution
class IndependentXYZ_Distribution : public Distribution_t<Point_t>
{
    private:
        std::shared_ptr<Distribution_t<double>> dist_x, dist_y, dist_z;
    
    public:
        IndependentXYZ_Distribution( std::shared_ptr<Distribution_t<double>> dx,
                                std::shared_ptr<Distribution_t<double>> dy, std::shared_ptr<Distribution_t<double>> dz, const std::string label = "" )
        : Distribution_t(label), dist_x(dx), dist_y(dy), dist_z(dz) {};
        ~IndependentXYZ_Distribution() {};
    
        Point_t sample();
    
        //dummy function. it is needed to pass the particle energy to sample fission
        double sampleFis( const Particle_t &P ){return 1.0; };
};


// Anisotropic direction distribution
class AnisotropicDirection_Distribution : public Distribution_t<Point_t>
{
    private:
        double sin_t;
    
        Point_t axis;
        std::shared_ptr<Distribution_t<double>> dist_mu;
    
    public:
        AnisotropicDirection_Distribution( Point_t p, std::shared_ptr<Distribution_t<double>> dmu, const std::string label = "" )
        : Distribution_t(label), axis(p), dist_mu(dmu)
        { axis.normalize(); sin_t = std::sqrt( 1.0 - axis.z * axis.z ); };
        ~AnisotropicDirection_Distribution() {};
    
        Point_t sample();
    
        //dummy function. it is needed to pass the particle energy to sample fission
        double sampleFis( const Particle_t &P ){return 1.0; };
    
};


#endif
