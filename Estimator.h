#ifndef _ESTIMATOR_HEADER_
#define _ESTIMATOR_HEADER_

#include <cmath>    // sqrt
#include <iostream> // cout
#include <iomanip>
#include <vector>   // vector
#include <fstream>  // ofstream
#include <memory>   // shared_ptr
#include <cstring>  // string
#include <sstream>  // ostringstream

#include "Particle.h"
#include "Geometry.h"


///////////////
/// Scoring ///
///////////////

// Score base class
class Score_t
{
	private:
		const std::string s_name;
	public:
		 Score_t( const std::string n ) : s_name(n) {};
		~Score_t() {};

		// Get score type name
		virtual std::string name() final { return s_name; }

		// Get score to be added at event
		virtual double add_score( const Particle_t& P, const double l = 0.0 ) = 0;
};

// Current (crossing)
class Current_Score : public Score_t
{
	public:
		 Current_Score() : Score_t( "Current" ) {};
		~Current_Score() {};

		double add_score( const Particle_t& P, const double l = 0.0 );
};

// Flux (path length)
class Flux_Score : public Score_t
{
	public:
		 Flux_Score() : Score_t( "Flux" ) {};
		~Flux_Score() {};

		double add_score( const Particle_t& P, const double track = 0.0 );
};

// Absorption (path length)
class Absorption_Score : public Score_t
{
	public:
		 Absorption_Score() : Score_t( "Abs. Rate" ) {};
		~Absorption_Score() {};

		double add_score( const Particle_t& P, const double track = 0.0 );
};

// Scatter (path length)
class Scatter_Score : public Score_t
{
	public:
		 Scatter_Score() : Score_t( "Scat. Rate" ) {};
		~Scatter_Score() {};

		double add_score( const Particle_t& P, const double track = 0.0 );
};

// Capture (path length)
class Capture_Score : public Score_t
{
	public:
		 Capture_Score() : Score_t( "Capt. Rate" ) {};
		~Capture_Score() {};

		double add_score( const Particle_t& P, const double track = 0.0 );
};

// Fission (path length)
class Fission_Score : public Score_t
{
	public:
		 Fission_Score() : Score_t( "Fis. Rate" ) {};
		~Fission_Score() {};

		double add_score( const Particle_t& P, const double track = 0.0 );
};

// Total (path length)
class Total_Score : public Score_t
{
	public:
		 Total_Score() : Score_t( "Tot. Rate" ) {};
		~Total_Score() {};

		double add_score( const Particle_t& P, const double track = 0.0 );
};



/////////////
/// Tally ///
/////////////

class Tally_t
{
	public:
		double hist    = 0.0; // History sum
		double sum     = 0.0; // Sum over all histories
	       	double squared = 0.0; // Sum of history sum squared over all histories
		double mean;          // Estimated mean
		double var;           // Estimated variance of the population
		double meanUncer;     // Uncertainty of the estimated mean
		                            // sqrt(variance of the estimated mean)
		double relUncer;      // Relative uncertainty of estimated mean
		double FOM;           // Figure of merit		
	
	public:
		 Tally_t() {};
		~Tally_t() {};
};



///////////
/// Bin ///
///////////

// Bin base class 
class Bin_t
{
	public:
		std::vector<double>                   grid;   // Bin grid
		std::vector<std::vector<Tally_t>>     tally;  // Bin tallies ( indexing --> [bin#][score#] )
		const int                             Nbin;   // # of bins
		std::vector<std::shared_ptr<Score_t>> scores; // Things to be scored
		const int                             Nscore; // # of scores
		const std::string                     unit;   // result unit, for report
	
	public:
		// Constructor: pass grid points and construct bin and their tallies
		 Bin_t( const std::vector<double> gr, const std::vector<Tally_t> total_tally, std::vector<std::shared_ptr<Score_t>>& s, 
				 const std::string str ) : 
			 grid(gr), Nbin(gr.size() - 1), unit(str), Nscore(s.size())
		{ 
			// Store scores
			scores = s;
			
			// Each bin is set with the same # of score tallies as the total tally of the estimator
			tally.resize( Nbin, total_tally );
		}
		~Bin_t() {};

		// Score bin
		virtual void score( const Particle_t& Pold, const Particle_t& P, const double track = 0.0 ) = 0;
};

// Energy bin
class Energy_Bin : public Bin_t
{
	public:
		// Constructor: pass grid points
		 Energy_Bin( const std::vector<double> grid, const std::vector<Tally_t> total_tally, std::vector<std::shared_ptr<Score_t>>& s ) : 
			 Bin_t(grid,total_tally,s,"eV") {};
		~Energy_Bin() {};

		void score( const Particle_t& Pold, const Particle_t& P, const double track = 0.0 );
};

// Time bin
class Time_Bin : public Bin_t
{
	public:
		// Constructor: pass grid points
		 Time_Bin( const std::vector<double> grid, const std::vector<Tally_t> total_tally, std::vector<std::shared_ptr<Score_t>>& s ) : 
			 Bin_t(grid,total_tally,s,"sec") {};
		~Time_Bin() {};

		void score( const Particle_t& Pold, const Particle_t& P, const double track = 0.0 );
};



/////////////////
/// Estimator ///
/////////////////

// Estimator base class
class Estimator_t
{
	protected:
		const std::string    e_name;     // Estimator name
		unsigned long long   nhist = 0;  // # of histories estimated

	public:
		// Constructor: pass the estimator name and bin grid
		 Estimator_t( const std::string n ) : e_name(n) {};
		~Estimator_t() {};

		// Add thing to be scored
		virtual void addScore( const std::shared_ptr<Score_t>& S ) = 0;

		// Set bin grid and corresponding tallies
		virtual void setBin( const std::string type, std::vector<double> bin ) = 0;
		
		// Score at events
		virtual void score( const Particle_t& Pold, const Particle_t& P, const double track = 0.0 ) = 0;
		
		// Closeout history
		virtual void endHistory() = 0;              
		
		// Report results
		virtual void report( std::ostringstream& output, const double tTime ) = 0;
};


// Generic estimator
////////////////////

class Generic_Estimator : public Estimator_t
{
	protected:
		std::vector<std::shared_ptr<Score_t>> scores;      // Things to be scored
		int                                   Nscore = 0;  // # of scores (things to be scored)
		std::vector<Tally_t>                  total_tally; // Total tallies [Nscore]
		
		std::shared_ptr<Bin_t>                bin;         // Estimator bin
		int                                   Nbin = 0;    // # of bins

	public:
		// Constructor: pass the estimator name
		 Generic_Estimator( const std::string n ) : Estimator_t(n) {};
		~Generic_Estimator() {};

		// Add thing to be scored and push new total tally
		void addScore( const std::shared_ptr<Score_t>& S )
		{ 
			// Push new score
			scores.push_back( S );
			Nscore++;
			
			// Push new total tally
			Tally_t T;
			total_tally.push_back( T );
			// Note: Index in total_tally corresponds to its score
		}

		// Set bin
		void setBin( const std::string type, const std::vector<double> grid )
		{
			Nbin = grid.size() - 1;
			if      ( type == "energy" ) { bin = std::make_shared<Energy_Bin> ( grid, total_tally, scores ); }
			else if ( type == "time" )   { bin = std::make_shared<Time_Bin>   ( grid, total_tally, scores ); }
		}
		
		// Update the sum and sum of squared, and restart history sum of a tally
		void tally_endHistory( Tally_t& T )
		{ 
			T.sum     += T.hist;
			T.squared += T.hist * T.hist;
			T.hist     = 0.0; 
		}
		
		// Compute mean, variance and mean statistical uncertainty of a tally
		void tally_stats( Tally_t& T, const double trackTime )
		{
			// Estimated mean of the population
			T.mean      = T.sum / nhist;
			// Estimated variance of the population
			T.var       = ( T.squared - nhist * T.mean*T.mean ) / ( nhist - 1.0 );
			// Uncertainty of the estimated mean
			T.meanUncer = sqrt( ( T.squared / nhist - T.mean*T.mean ) / ( nhist - 1.0 ) );
			// Relative uncertainty of estimated mean
			T.relUncer  = T.meanUncer / T.mean;
			// Figure of merit
			T.FOM       = 1.0 / ( T.relUncer*T.relUncer * trackTime );
		};
		
		// Score at events
		void score( const Particle_t& Pold, const Particle_t& P, const double track = 0.0 );

		// Closeout history
		// Update the sum and sum of squared, and restart history sum of all tallies
		void endHistory()
		{
			nhist++;
			for ( int i = 0 ; i < Nscore ; i++ )
			{
				// Total tally
				tally_endHistory( total_tally[i] );
				
				// Bin tally
				for ( int j = 0 ; j < Nbin ; j++ )
				{
					tally_endHistory( bin->tally[j][i] );
				}
			}
		}

		// Report results
		void report( std::ostringstream& output, const double trackTime )
		{
			// Compute mean, variance, and uncertainty of all tallies
			for ( int i = 0 ; i < Nscore ; i++ )
			{
				// Total tally
				tally_stats( total_tally[i], trackTime );
				
				// Bin tally
				for ( int j = 0 ; j < Nbin ; j++ )
				{
					tally_stats( bin->tally[j][i], trackTime );
				}	
			}
			
			output << "\n\n";
			output << "Estimator report: " + e_name + "\n";
			for ( int i = 0 ; i < e_name.length()+18 ; i++ ) { output << "="; }
			output << "\n";
			
			// Total tallies
			output << "Total tallies,\n";
			for ( int i = 0 ; i < Nscore ; i++ )
			{
				output << "  " + scores[i]->name() + ":\n";
				output << "  -> Mean     = " << std::scientific << total_tally[i].mean;
			       	output << "  +/-  " << std::scientific << total_tally[i].meanUncer;
				output << "  (" << std::scientific << total_tally[i].relUncer * 100.0;
				output << "%)\n";
				output << "  -> Variance = " << std::scientific << total_tally[i].var;
			       	output << "\n";
				output << "  [F.O.M.: " << std::scientific << total_tally[i].FOM;
			        output << "]\n\n\n";
			}

			// Bin tallies
			if ( Nbin != 0 )
			{
				output << "Bin tallies," << std::endl;
		
				output << "bin#\t";
				output << "lower(" + bin->unit + ")\tupper(" + bin->unit + ")\t";
				for ( int i = 0 ; i < Nscore ; i++ ) 
				{ 
					output << std::setw(12) << std::left << scores[i]->name() << "\t"; 
					output << std::setw(12) << std::left << "_uncertainty\t"; 
				}
				output << "\n";
			
				output << "----\t" << "------------\t" << "------------\t";
				for ( int i = 0 ; i < Nscore ; i++ ) 
				{ 
					output << "------------\t" << "------------\t"; 
				}
				output << "\n";
			
				for ( int j = 0 ; j < Nbin ; j++ )
				{
					output << j+1 << "\t";
					output << std::scientific << bin->grid[j];
				        output << "\t" << std::scientific << bin->grid[j+1]; 
					output << "\t";

					for ( int i = 0 ; i < Nscore ; i++ )
					{
						output << std::scientific << bin->tally[j][i].mean;
					        output << "\t" << std::scientific << bin->tally[j][i].meanUncer;
					        output << "\t";
					}
					output << "\n";
				}
			}
		}
};



/// Miscellaneous Estimator
///////////////////////////

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
		virtual void score( const Particle_t& Pold, const Particle_t& P, const double track = 0.0 ) = 0;
	
		// Report results
		virtual void report( std::ostringstream& output, const double tTime ) = 0; 

		// Normalize the recorded PMF
		virtual void normalize() final
		{ for ( int i = 0; i < pmf.size(); i++ ) { pmf[i] = pmf[i]/nhist; } }
		
		// Add thing to be scored [unsuported]
		void addScore( const std::shared_ptr<Score_t>& S ) { return; }

		// Set bin grid and corresponding tallies [unsuported]
		void setBin( const std::string type, std::vector<double> bin ) { return; };
};


// Crossing surface PMF estimator (a.k.a. surface counting)
class Surface_PMF_Estimator : public UInteger_PMF_Estimator 
{
	public:
		// Initialize variables at construction 
		 Surface_PMF_Estimator( const std::string n ) : UInteger_PMF_Estimator(n) {};
		~Surface_PMF_Estimator() {};

		// Score at events
		void score( const Particle_t& Pold, const Particle_t& P, const double track = 0.0 );

		// Report results
		// output.txt file providing the PMF is created (or overwritten)
		void report( std::ostringstream& output, const double tTime );
};

#endif

