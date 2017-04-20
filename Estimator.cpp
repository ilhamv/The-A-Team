#include <cmath>
#include <iostream>
#include <sstream>  // ostringstream

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

// Production (path length)
//double Production_Score::add_score( const Particle_t& P, const double track /*= 0.0*/ )
//{ return P.region()->nuSigmaF( P.energy() ) * track * P.weight(); }

// Total (path length)
double Total_Score::add_score( const Particle_t& P, const double track /*= 0.0*/ )
{ return P.region()->SigmaT( P.energy() ) * track * P.weight(); }



///////////
/// Bin ///
///////////

// Energy Bin (multi-score)
void Energy_Bin::score( const Particle_t& Pold, const Particle_t& P, const std::vector<double>& grid, const double track /*= 0.0*/ )
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

// Time Bin (multi-score)
void Time_Bin::score( const Particle_t& Pold, const Particle_t& P, const std::vector<double>& grid, const double track /*= 0.0*/ )
{
	// Since score migth be distributed across bins...
	
	// Pair of bin location and corresponding track to be scored
	std::vector< std::pair<int,double> > loc_track;

	// Search bin location
	int loc1    = Binary_Search( Pold.time(), grid ); // before track generation
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
		
// Add thing to be scored and push new total tally
void Generic_Estimator::addScore( const std::shared_ptr<Score_t>& S )
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
void Generic_Estimator::setBin( const std::string type, const std::vector<double> gr )
{
	grid = gr;
	Nbin = grid.size() - 1;
	if      ( type == "energy" ) { bin = std::make_shared<Energy_Bin> ( grid, total_tally, scores ); }
	else if ( type == "time" )   { bin = std::make_shared<Time_Bin>   ( grid, total_tally, scores ); }
}
	
// Update the sum and sum of squared, and restart history sum of a tally
void Generic_Estimator::tally_endHistory( Tally_t& T )
{ 
	T.sum     += T.hist;
	T.squared += T.hist * T.hist;
	T.hist     = 0.0; 
}
	
// Compute mean, variance and mean statistical uncertainty of a tally
void Generic_Estimator::tally_stats( Tally_t& T, const double trackTime )
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
void Generic_Estimator::score( const Particle_t& Pold, const Particle_t& P, const double track /*= 0.0*/ )
{
        // Total tallies
        for ( int i = 0 ; i < Nscore ; i++ )
        {
                total_tally[i].hist += scores[i]->add_score( P, track );
        }

	// Bin tallies, if any
	if ( Nbin != 0 )
	{ bin->score( Pold, P, grid, track ); }
}

// Closeout history
// Update the sum and sum of squared, and restart history sum of all tallies
void Generic_Estimator::endHistory()
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
void Generic_Estimator::report( std::ostringstream& output, const double trackTime )
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
		output << "  -> Mean     = " << total_tally[i].mean;
	       	output << "  +/-  " << total_tally[i].meanUncer;
		output << "  (" << std::defaultfloat << total_tally[i].relUncer * 100.0;
		output << "%)\n";
		output << "  -> Variance = " << total_tally[i].var;
	       	output << "\n";
		output << "  [F.O.M.: " << total_tally[i].FOM;
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
			output << std::setw(12) << std::left << "uncertainty\t"; 
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
			output << j+1;
			output << "\t" << std::setw(12) << std::left << grid[j];
		        output << "\t" << std::setw(12) << std::left << grid[j+1]; 
			output << "\t";

			for ( int i = 0 ; i < Nscore ; i++ )
			{
				output << std::setw(12) << bin->tally[j][i].mean << "\t";
			        output << std::setw(12) << bin->tally[j][i].meanUncer << "\t";
			}
			output << "\n";
		}
	}
}



// Homogenized MG Constant Generator
////////////////////////////////////

// Score at events
void MGXS_Estimator::score( const Particle_t& Pold, const Particle_t& P, const double track /*= 0.0*/ )
{
	// Flux, Capture, Fission, nuFission, Total, Scater Tallies
	bin->score( Pold, P, grid, track );
	
	// Scattering matrix Tallies
	// First, simulate scattering reaction to get post-scattering particle information
	Particle_t P_final = P;
	P.region()->simulate_scatter( P_final );
	
	// Then, find the final energy bin location
	const int loc = Binary_Search( P_final.energy(), grid );
	
	// Final energy bin location is set, now score into the initial energy bin location
	matrix_bin[loc]->score( Pold, P, grid, track );
}


// Closeout history
// Update the sum and sum of squared, and restart history sum of all tallies
void MGXS_Estimator::endHistory()
{
	nhist++;
	// Flux, Capture, Fission, nuFission, Total, Scater Tallies
	for ( int i = 0 ; i < Nscore ; i++ )
	{
		for ( int j = 0 ; j < Nbin ; j++ )
		{
			tally_endHistory( bin->tally[j][i] );
		}
	}
	
	// Scattering matrix Tallies
	for ( int i = 0 ; i < Nbin ; i++ )
	{
		for ( int j = 0 ; j < Nbin ; j++ )
		{
			tally_endHistory( matrix_bin[i]->tally[j][0] );
		}
	}
}


// Report results
void MGXS_Estimator::report( std::ostringstream& output, const double trackTime )
{
	// Compute mean, variance, uncertainty, and F.O.M of all tallies
	// Flux, Capture, Fission, nuFission, Total, Scater Tallies
	for ( int i = 0 ; i < Nscore ; i++ )
	{
		for ( int j = 0 ; j < Nbin ; j++ )
		{
			tally_stats( bin->tally[j][i], trackTime );
		}	
	}	
	// Scattering matrix Tallies
	for ( int i = 0 ; i < Nbin ; i++ )
	{
		for ( int j = 0 ; j < Nbin ; j++ )
		{
			tally_stats( matrix_bin[i]->tally[j][0], trackTime );
		}	
	}	
	

	// Convert tally (except flux) into homogenized multigroup constant [mean and uncertainty]
	for ( int j = 0 ; j < Nbin ; j++ )
	{
		const double B = bin->tally[j][0].meanUncer / bin->tally[j][0].mean;
		for ( int i = 1 ; i < Nscore ; i++ )
		{
			const double A = bin->tally[j][i].meanUncer / bin->tally[j][i].mean;
			bin->tally[j][i].meanUncer = std::sqrt( A*A + B*B );

			bin->tally[j][i].mean = bin->tally[j][i].mean / bin->tally[j][0].mean;
			bin->tally[j][i].meanUncer = bin->tally[j][i].meanUncer * bin->tally[j][i].mean;
		}	
		for ( int i = 0 ; i < Nbin ; i++ )
		{
			const double A = matrix_bin[i]->tally[j][0].meanUncer / matrix_bin[i]->tally[j][0].mean;
			matrix_bin[i]->tally[j][0].meanUncer = std::sqrt( A*A + B*B );

			matrix_bin[i]->tally[j][0].mean = matrix_bin[i]->tally[j][0].mean / bin->tally[j][0].mean;
			matrix_bin[i]->tally[j][0].meanUncer = matrix_bin[i]->tally[j][0].meanUncer * matrix_bin[i]->tally[j][0].mean;
		}	
	}	
	

	output << "\n\n";
	output << "Estimator report: " + e_name + "\n";
	for ( int i = 0 ; i < e_name.length()+18 ; i++ ) { output << "="; }
	output << "\n";
	
	// Bin tallies
	if ( Nbin != 0 )
	{
		output << "Homogenized MG Constants," << std::endl;

		output << "g\t";
		output << "lower(" + bin->unit + ")\tupper(" + bin->unit + ")\t";
		for ( int i = 1 ; i < Nscore ; i++ ) 
		{ 
			output << std::setw(12) << std::left << scores[i]->name() << "\t"; 
			output << std::setw(12) << std::left << "uncertainty\t"; 
		}
		for ( int j = 0 ; j < Nbin ; j++ )
		{
			const std::string g_prime = std::to_string(j+1);
			output << std::setw(12) << std::left << "g -> " + g_prime << "\t"; 
			output << std::setw(12) << std::left << "uncertainty\t"; 
		}
		output << "\n";
	
		output << "----\t" << "------------\t" << "------------\t";
		for ( int i = 1 ; i < Nscore ; i++ ) 
		{ 
			output << "------------\t" << "------------\t"; 
		}
		for ( int j = 0 ; j < Nbin ; j++ )
		{ 
			output << "------------\t" << "------------\t"; 
		}
		output << "\n";
	
		for ( int j = 0 ; j < Nbin ; j++ )
		{
			output << j+1;
			output << "\t" << std::setw(12) << std::left << grid[j];
		        output << "\t" << std::setw(12) << std::left << grid[j+1]; 
			output << "\t";

			for ( int i = 1 ; i < Nscore ; i++ )
			{
				output << std::setw(12) << bin->tally[j][i].mean << "\t";
			        output << std::setw(12) << bin->tally[j][i].meanUncer << "\t";
			}
			for ( int i = 0 ; i < Nbin ; i++ )
			{
				output << std::setw(12) << matrix_bin[i]->tally[j][0].mean << "\t"; 
				output << std::setw(12) << matrix_bin[i]->tally[j][0].meanUncer << "\t"; 
			}
			output << "\n";
		}
	}
}


// Miscellaneous Estimator
//////////////////////////

// Surface PMF estimator (surface counting)
void Surface_PMF_Estimator::score( const Particle_t& Pold, const Particle_t& P, const double track /*= 0.0*/ )
{ tally_hist++;}

void Surface_PMF_Estimator::report( std::ostringstream& output, const double tTime ) 
{
	normalize();                        // Normalize the recorded PMF
	stats();                            // Compute mean, variance and PMF statistical uncertainty
	
	// Printouts (Uncertainty is only printed out in the output.txt file)
	output << "\n\n";
	output << "Estimator report: " << e_name << "\n";
	for ( int i = 0 ; i < e_name.length()+18 ; i++ ) { output << "="; }
	output << "\n";
	output << "PMF of # of particle crossing surface ";
	output << ":\nx\t" << std::setw(12) << std::left << "P(x)\t" << "uncertainty\n";
	output << "---\t------------\t------------\n";
	for ( int i = 0 ; i < pmf.size() ; i++ )
	{ output << i << "\t" << std::setw(12) << std::left << pmf[i] << "\t" << pmfUncer[i] << std::endl;  }
	
	output << "\nMean    : " << mean << std::endl;
	output << "Variance: "   << var  << std::endl;	
}		
