#ifndef _ESTIMATOR_HEADER_
#define _ESTIMATOR_HEADER_

#include <cmath>    // sqrt
#include <iostream> // cout
#include <vector>   // vector
#include <fstream>  // ofstream
#include <memory>   // shared_ptr
#include <cstring>  // string

#include "Particle.h"
#include "Region.h"
#include "Surface.h"

// Forward declaration
class Surface_t;


// Estimator base class
class Estimator_t 
{
	protected:
		std::vector<std::string> o_names;     // Estimated object names
		const std::string        e_name;      // Estimator name
		unsigned long long       nhist   = 0; // # of histories estimated
		virtual void             stats() = 0; // Compute mean and variance
	public:
		// Constructor: pass the estimator name
		 Estimator_t( const std::string n ) : e_name(n) {};
		~Estimator_t() {};

		// Add estimated object name, for reference in the report
		virtual void addName( const std::string s ) final 
		{ o_names.push_back( s );}
		
		// Score at events
		virtual void score( const Particle_t& P, const double p = 0.0 ) = 0; 
		
		// Closeout history
		virtual void endHistory() = 0;              
		
		// Report results
		virtual void report( const double tTime ) = 0;
};


// Single valued estimator base class
class Single_Valued_Estimator : public Estimator_t 
{
	protected:
		// History sum, sum over all histories, and sum of history sum squared over all histories
		double tally_hist    = 0.0; // History sum
		double tally_sum     = 0.0; // Sum over all histories
	       	double tally_squared = 0.0; // Sum of history sum squared over all histories
		double mean;  		    // Estimated mean
		double var;       	    // Estimated variance of the population
		double meanUncer; 	    // Uncertainty of the estimated mean
		                            // sqrt(variance of the estimated mean)

		// Compute mean, variance and mean statistical uncertainty
		virtual void stats() final 
		{
			// Estimated mean of the population
			mean      = tally_sum / nhist;
			// Estimated variance of the population
			var       = ( tally_squared - nhist * mean*mean ) / ( nhist - 1.0 );
			// Uncertainty of the estimated mean
			meanUncer = sqrt( ( tally_squared / nhist - mean*mean ) / ( nhist - 1.0 ) );
		};
	
	public:
		// Constructor: pass the name
		Single_Valued_Estimator( const std::string n ) : Estimator_t(n) {};
		~Single_Valued_Estimator() {};

		// Closeout history
		// Update the sum and sum of squared, and restart history sum
		virtual void endHistory() final 
		{ 
			nhist++;
			tally_sum     += tally_hist;
			tally_squared += tally_hist * tally_hist;
			tally_hist     = 0.0; 
		}

		// Score at events
		virtual void score( const Particle_t& P, const double p = 0.0 ) = 0;
		
		// Report results
		virtual void report( const double tTime ) = 0;                   
};


// Crossing surface estimator
class Surface_Current_Estimator : public Single_Valued_Estimator 
{
	public:
		// Constructor: pass the name
		 Surface_Current_Estimator( const std::string n ) : Single_Valued_Estimator( n ) {};
		~Surface_Current_Estimator() {};

		// Score at events
		void score( const Particle_t& P, const double null = 0.0 ); 

		// Report results
		void report( const double tTime );
};


// Volume-averaged scalar flux in a cell
class Region_Flux_Estimator : public Single_Valued_Estimator 
{
  	private:
		std::shared_ptr<Region_t> R;                  // Estimated region
		const bool capture, scatter, fission, absorp; // Reaction rate estimator flag
  	
	public:
		// Constructor: pass the name, region pointer, and reaction flags
     		 Region_Flux_Estimator( const std::string n, const std::shared_ptr<Region_t>& reg, const bool p1, const bool p2, const bool p3, const bool p4 ) : 
       			Single_Valued_Estimator(n), R(reg), scatter(p1), capture(p2), fission(p3), absorp(p4) {};
    		~Region_Flux_Estimator() {};

    		// Score at events
		void score( const Particle_t& P, const double path_length = 0.0 );

    		// Report results
		void report( const double tTime ); 
};


// Unsigned integer PMF estimator base class
class UInteger_PMF_Estimator : public Estimator_t 
{
	protected:
		int tally_hist = 0;           // History sum
		double mean, var;             // Esimated mean and variance of the population
		std::vector<double> pmf;      // Recorded PMF, it expands as new valid unsigned integer arises
		                              // the vector index means the unsigned integer random variable
		std::vector<double> pmfUncer; // Recorded PMF uncertainty

		// Compute mean, variance and PMF statistical uncertainty
		virtual void stats() final 
		{
			// Create PMF uncertainty vector with appropriate size
			pmfUncer.resize( pmf.size() );
			
			double s1=0,s2=0; // first and second moment accumulator
			for ( int i = 0; i < pmf.size(); i++ )
			{
				s1 += i * pmf[i];
				s2 += i*i *pmf[i];
				pmfUncer[i] = sqrt( pmf[i] / nhist );
			}
			
			mean      = s1;
			var       = s2 - s1*s1; 
		};
	
	public:
		// Constructor: pass the name
		UInteger_PMF_Estimator( const std::string n ) : Estimator_t(n) {};
		~UInteger_PMF_Estimator() {};

		// Closeout History
		// Update the PMF (not normalized yet)
		virtual void endHistory() final 
		{ 
			nhist++;
			
			// Check if the sampled unsigned integer (tally_hist) is inside the PMF already
			if ( pmf.size() > tally_hist )
			{
				// It is: increment its probability bin
				pmf[tally_hist]++;	
			}
			// It isn't yet: resize the PMF vector up to the sampled unsigned integer index
			// vector size becomes tally_hist+1, then increment the sampled unsigned integer bin
			else { pmf.resize( tally_hist + 1 ); pmf.back()++; }
			tally_hist     = 0; 
		};

		// Score at events
		virtual void score( const Particle_t& P, const double p = 0.0 ) = 0; 
	
		// Report results
		virtual void report( const double tTime ) = 0; 

		// Normalize the recorded PMF
		virtual void normalize() final
		{ for ( int i = 0; i < pmf.size(); i++ ) { pmf[i] = pmf[i]/nhist; } }
};


// Crossing surface PMF estimator
class Surface_PMF_Estimator : public UInteger_PMF_Estimator 
{
	public:
		// Initialize variables at construction 
		 Surface_PMF_Estimator( const std::string n ) : UInteger_PMF_Estimator(n) {};
		~Surface_PMF_Estimator() {};

		// Score at events
		void score( const Particle_t& P, const double null = 0.0 );

		// Report results
		// output.txt file providing the PMF is created (or overwritten)
		void report( const double tTime );
};

#endif

