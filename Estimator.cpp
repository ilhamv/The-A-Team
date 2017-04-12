#include <cmath>
#include <iostream>

#include "Estimator.h"
#include "Geometry.h"
#include "Solver.h"  // Binary_Search


///////////////
/// Scoring ///
///////////////

// Current (crossing)
double Current_Score::add_score( const Particle_t& P, const double l /*= 0.0*/ )
{ return P.weight(); }

// Flux (path length)
double Flux_Score::add_score( const Particle_t& P, const double track /*= 0.0*/ )
{ return track * P.weight(); }

// Absorption (path length)
double Absorption_Score::add_score( const Particle_t& P, const double track /*= 0.0*/ )
{ return ( P.region()->SigmaC( P.energy() ) + P.region()->SigmaF( P.energy() ) ) * track * P.weight(); }

// Scatter (path length)
double Scatter_Score::add_score( const Particle_t& P, const double track /*= 0.0*/ )
{ return P.region()->SigmaS( P.energy() ) * track * P.weight(); }

// Capture (path length)
double Capture_Score::add_score( const Particle_t& P, const double track /*= 0.0*/ )
{ return P.region()->SigmaC( P.energy() ) * track * P.weight(); }

// Fission (path length)
double Fission_Score::add_score( const Particle_t& P, const double track /*= 0.0*/ )
{ return P.region()->SigmaF( P.energy() ) * track * P.weight(); }

// Total (path length)
double Total_Score::add_score( const Particle_t& P, const double track /*= 0.0*/ )
{ return P.region()->SigmaT( P.energy() ) * track * P.weight(); }



///////////
/// Bin ///
///////////

// Energy Bin
void Energy_Bin::score( const Particle_t& Pold, const Particle_t& P, const double track /*= 0.0*/ )
{
	// Search bin location	
	const int loc = Binary_Search( P.energy(), grid );

	// Check if inside the grids
	if ( loc >= 0 && loc < Nbin )
	{
		// Iterate over scores
		for ( int i = 0 ; i < Nscore ; i++ )
		{
			tally[loc][i].hist += scores[i]->add_score( P, track );
		}
		// Note: In case of using cross section table, 
		// might want to pass an index pointing to the XSec table location
		// to avoid XS search for each reaction score!
	}
	// Note: value exactly at the grid point constributes to bin whose upper bound is the value
}

// Time Bin
void Time_Bin::score( const Particle_t& Pold, const Particle_t& P, const double track /*= 0.0*/ )
{
	// Since score migth be distributed across bins...
	
	// Pair of bin location and corresponding track to be scored
	std::vector< std::pair<int,double> > loc_track;

	// Search bin location
	int loc1    = Binary_Search( Pold.time(), grid );    // before track generation
	int loc2    = Binary_Search( P.time()   , grid ); // after
	
	// 1 bin spanned [need to consider if it is outside the time grid]
	if ( loc1 == loc2 ) 
	{
		if ( loc1 >= 0 && loc1 < Nbin )
		{ loc_track.push_back( std::make_pair( loc1, track ) ); }
	}
	
	// >1 bin spanned
	else
	{				
		int    num_bin = loc2 - loc1 - 1; // # of full bins spanned
		double new_track;                 // to hold bin track
		
		// First partial bin [need to consider if it's outside]
		if ( loc1 >= 0 )
		{
			new_track = ( grid[loc1+1] - Pold.time() ) * P.speed();
			loc_track.push_back( std::make_pair( loc1, new_track ) );
		}
		
		// Intermediate full bin
		for ( int i = 1 ; i <= num_bin ; i++ )
		{
			new_track = ( grid[loc1+i+1] - grid[loc1+i] ) * P.speed();
			loc_track.push_back( std::make_pair( loc1+i, new_track ) );
		}
		
		// Last partial bin [consider if it's outside]
		if ( loc2 < Nbin )
		{
			new_track = ( P.time() - grid[loc2] ) * P.speed();
			loc_track.push_back( std::make_pair( loc2, new_track ) );
		}
		// Note: this algorithm covers the following "extreme" cases
		// 	loc1 : <lowest_grid  or at grid point
		// 	loc2 : >highest_grid or at grid point
	}

	// Iterate over scores
	for ( int i = 0 ; i < Nscore ; i++ )
	{
		// Iterate over scored bins
		for ( auto& LnT : loc_track )
		{
			tally[LnT.first][i].hist += scores[i]->add_score( P, LnT.second );
		}	
		// Note: In case of using cross section table, 
		// might want to pass an index pointing to the location in XSec table
		// to avoid XS search for each reaction score!
	}
}




/////////////////
/// Estimator ///
/////////////////

// Generic estimator
////////////////////

// Score at events
void Generic_Estimator::score( const Particle_t& Pold, const Particle_t& P, const double track /*= 0.0*/ )
{
        // Total tallies
        for ( int i = 0 ; i < Nscore ; i++ )
        {
                total_tally[i].hist += scores[i]->add_score( P, track );
        }

	// Bin tallies, if any
	if ( Nbin != 0 )
	{ bin->score( Pold, P, track ); }
}



// Miscellaneous Estimator
//////////////////////////

// Surface PMF estimator (surface counting)
void Surface_PMF_Estimator::score( const Particle_t& Pold, const Particle_t& P, const double track /*= 0.0*/ )
{ tally_hist++;}

void Surface_PMF_Estimator::report( const std::string simName, const double tTime ) 
{
	normalize();                        // Normalize the recorded PMF
	stats();                            // Compute mean, variance and PMF statistical uncertainty
	std::ofstream file( simName + " - " + e_name + ".txt" ); // Create .txt file
	
	// Printouts (Uncertainty is only printed out in the output.txt file)
	std::cout<< std::endl << std::endl;
	std::cout<< "Estimator report: " << e_name << std::endl;
	for ( int i = 0 ; i < e_name.length()+18 ; i++ ) { std::cout<< "="; }
	std::cout<< std::endl;
	std::cout << "PMF of # of particle crossing surface ";
	     file << "PMF of # of particle crossing surface ";
	std::cout << ":\nx\tP(x)\n---\t-------"                        << std::endl;
	     file << ":\nx\tP(x)\t+/-\n---\t-------\t------------" << std::endl;
	for ( int i = 0 ; i < pmf.size() ; i++ )
	{ 
		std::cout<< i << "\t" << pmf[i]                          << std::endl; 
		     file<< i << "\t" << pmf[i] << "\t" << pmfUncer[i] << std::endl; 
	}
	
	std::cout << "\nMean    : " << mean     << std::endl;
	std::cout << "Variance: "   << var << std::endl;	
	     file << "\nMean\t"     << mean     << std::endl;
	     file << "Variance\t"   << var << std::endl;	
}		
