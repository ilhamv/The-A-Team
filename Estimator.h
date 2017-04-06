#ifndef _ESTIMATOR_HEADER_
#define _ESTIMATOR_HEADER_

#include <cmath>    // sqrt
#include <iostream> // cout
#include <vector>   // vector
#include <fstream>  // ofstream
#include <memory>   // shared_ptr
#include <cstring>  // string

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

		virtual std::string name() final { return s_name; }
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
	
	public:
		 Tally_t() {};
		~Tally_t() {};
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

		// Score at events
		virtual void score( const Particle_t& P, const double p = 0.0, bool reg_flag = false, const double t_old = 0.0 ) = 0; 
		
		// Add thing to be scored
		virtual void addScore( const std::shared_ptr<Score_t>& S ) = 0;

		// Set bin grid and corresponding tallies
		virtual void setBin( const std::string type, std::vector<double> bin, const bool br ) = 0;
		
		// Closeout history
		virtual void endHistory() = 0;              
		
		// Report results
		virtual void report( const std::string simName, const double tTime ) = 0;
};


// Generic estimator
class Generic_Estimator : public Estimator_t
{
	protected:
		std::vector<std::shared_ptr<Score_t>> scores;             // Things to be scored
		std::vector<Tally_t>                  total_tally;        // Total tallies [Nscore]
		bool                                  bin_active = false; // Bin active flag
		std::vector<double>                   bin_grid;           // Bin grid
		std::vector<std::vector<Tally_t>>     bin_tally;          // Bin tallies   [(Ngrid-1) x Nscore]
		std::string                           bin_type;           // Bin type (time or energy)
		int                                   Nscore = 0;         // Score size
		int                                   Nbin = 0;           // Bin size
		bool                                  bin_report;         // Bin report on monitor flag

	public:
		// Constructor: pass the estimator name and bin grid
		 Generic_Estimator( const std::string n ) : Estimator_t(n) {};
		~Generic_Estimator() {};

		// Add thing to be scored and push new total tally
		void addScore( const std::shared_ptr<Score_t>& S )
		{ 
			scores.push_back( S );
			Tally_t T;
			total_tally.push_back( T );
			Nscore++;
		}

		// Set bin grid and corresponding tallies
		void setBin( const std::string type, std::vector<double> bin, const bool br )
		{
			bin_active = true;
			bin_type   = type;
			bin_report = br;
			bin_grid   = bin;
			Nbin       = bin.size() - 1;
			bin_tally.resize( Nbin, total_tally );
		}
		
		// Update the sum and sum of squared, and restart history sum of a tally
		void tally_endHistory( Tally_t& T )
		{ 
			T.sum     += T.hist;
			T.squared += T.hist * T.hist;
			T.hist     = 0.0; 
		}
		
		// Compute mean, variance and mean statistical uncertainty of a tally
		void tally_stats( Tally_t& T )
		{
			// Estimated mean of the population
			T.mean      = T.sum / nhist;
			// Estimated variance of the population
			T.var       = ( T.squared - nhist * T.mean*T.mean ) / ( nhist - 1.0 );
			// Uncertainty of the estimated mean
			T.meanUncer = sqrt( ( T.squared / nhist - T.mean*T.mean ) / ( nhist - 1.0 ) );
			// Relative uncertainty of estimated mean
			T.relUncer  = T.meanUncer / T.mean;
		};
		
		// Score at events
		void score( const Particle_t& P, const double track = 0.0, bool reg_flag = false, const double t_old = 0.0 );

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
					tally_endHistory( bin_tally[j][i] );
				}
			}
		}

		// Report results
		void report( const std::string simName, const double tTime )
		{
			// Compute mean, variance, and uncertainty of all tallies
			for ( int i = 0 ; i < Nscore ; i++ )
			{
				// Total tally
				tally_stats( total_tally[i] );
				// Bin tally
				for ( int j = 0 ; j < Nbin ; j++ )
				{
					tally_stats( bin_tally[j][i] );
				}	
			}
			
			std::cout<< std::endl << std::endl;
			std::cout<< "Estimator report: " << e_name << std::endl;
			for ( int i = 0 ; i < e_name.length()+18 ; i++ ) { std::cout<< "="; }
			std::cout<< std::endl;
			
			// Total tallies
			std::cout<< "Total tallies," << std::endl;
			for ( int i = 0 ; i < Nscore ; i++ )
			{
				std::cout<< "  " << scores[i]->name() << ":" << std::endl;
				std::cout<< "  -> Mean     = " << total_tally[i].mean << "  +/-  " << total_tally[i].meanUncer
					              << "  (" << total_tally[i].relUncer * 100.0 << "%)" << std::endl;
				std::cout<< "  -> Variance = " << total_tally[i].var << std::endl;
				std::cout<< "  [F.O.M.: " << 1.0 / ( total_tally[i].relUncer*total_tally[i].relUncer * tTime ) << "]" << std::endl << std::endl;
			}

			// Bin tallies
			if ( bin_report )
			{
	
				// Printouts (Uncertainty is only printed out in the output file)
				std::cout << std::endl;
				std::cout << "Bin tallies," << std::endl;
			
				std::cout<<"bin#\t";
				if      ( bin_type == "energy" ) { std::cout<<"lower(eV)\tupper(eV)\t"; }
				else if ( bin_type == "time" )   { std::cout<<"lower(sec)\tupper(sec)\t"; }
				for ( int i = 0 ; i < Nscore ; i++ )
				{
					std::cout << scores[i]->name() << "\t\t\t";
				}
				std::cout << std::endl;
				
				std::cout << "----\t";
				std::cout << "----------\t";
				std::cout << "----------\t";
				for ( int i = 0 ; i < Nscore ; i++ )
				{
					std::cout << "------------\t\t";
				}
				std::cout << std::endl;
				
				for ( int j = 0 ; j < Nbin ; j++ )
				{
					std::cout<< j+1 << "\t";
					std::cout<< bin_grid[j] << "\t\t" << bin_grid[j+1] << "\t\t";
					for ( int i = 0 ; i < Nscore ; i++ )
					{
						std::cout << bin_tally[j][i].mean << "\t\t";
					}
					std::cout << std::endl;
				}
			}

			if ( bin_active )
			{
				std::ofstream file( simName + " - " + e_name + ".txt" ); // Create .txt file
				// Printouts (Uncertainty is only printed out in the output file)
				file << std::endl;
				file << "Bin tallies," << std::endl;
		
				file <<"bin#\t";
				if      ( bin_type == "energy" ) { file<<"lower(eV)\tupper(eV)\t"; }
				else if ( bin_type == "time" )   { file<<"lower(sec)\tupper(sec)\t"; }
				for ( int i = 0 ; i < Nscore ; i++ )
				{
					file<< scores[i]->name() << "\t";
				}
				file << std::endl;
			
				file << "----\t";
				file << "----------\t";
				file << "----------\t";
				for ( int i = 0 ; i < Nscore ; i++ )
				{
					file << "------------\t";
				}
				file << std::endl;
			
				for ( int j = 0 ; j < Nbin ; j++ )
				{
					file<< j+1 << "\t";
					file<< bin_grid[j] << "\t" << bin_grid[j+1] << "\t";
					for ( int i = 0 ; i < Nscore ; i++ )
					{
						file << bin_tally[j][i].mean << "\t";
					}
					file << std::endl;
				}
			}
		}
};



/////////////////////
/// Miscellaneous ///
////////////////////


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
		virtual void score( const Particle_t& P, const double p = 0.0 , bool reg_flag = false, const double t_old = 0.0 ) = 0; 
	
		// Report results
		virtual void report( const std::string simName, const double tTime ) = 0; 

		// Normalize the recorded PMF
		virtual void normalize() final
		{ for ( int i = 0; i < pmf.size(); i++ ) { pmf[i] = pmf[i]/nhist; } }
		
		// Add thing to be scored [unsuported]
		void addScore( const std::shared_ptr<Score_t>& S ) { return; }

		// Set bin grid and corresponding tallies [unsuported]
		void setBin( const std::string type, std::vector<double> bin, const bool br ) { return; };
};


// Crossing surface PMF estimator (a.k.a. surface counting)
class Surface_PMF_Estimator : public UInteger_PMF_Estimator 
{
	public:
		// Initialize variables at construction 
		 Surface_PMF_Estimator( const std::string n ) : UInteger_PMF_Estimator(n) {};
		~Surface_PMF_Estimator() {};

		// Score at events
		void score( const Particle_t& P, const double null = 0.0, bool reg_flag = false, const double t_old = 0.0 );

		// Report results
		// output.txt file providing the PMF is created (or overwritten)
		void report( const std::string simName, const double tTime );
};

#endif

