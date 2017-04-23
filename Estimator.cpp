#include <cmath>
#include <iostream>
#include <sstream>  // ostringstream

#include "Estimator.h"
#include "Geometry.h"
#include "Solver.h"  // Binary_Search, Linterpolate


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
void Energy_Bin::score( const Particle_t& P, const std::vector<double>& grid, const double told, const double track /*= 0.0*/ )
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
void Time_Bin::score( const Particle_t& P, const std::vector<double>& grid, const double told, const double track /*= 0.0*/ )
{
	// Since score migth be distributed across bins...
	
	// Pair of bin location and corresponding track to be scored
	std::vector< std::pair<int,double> > loc_track;

	// Search bin location
	int loc1    = Binary_Search( told, grid ); // before track generation
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
			new_track = ( grid[loc1+1] - told ) * P.speed();
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
void Generic_Estimator::score( const Particle_t& P, const double told, const double track /*= 0.0*/ )
{
        // Total tallies
        for ( int i = 0 ; i < Nscore ; i++ )
        {
                total_tally[i].hist += scores[i]->add_score( P, track );
        }

	// Bin tallies, if any
	if ( Nbin != 0 )
	{ bin->score( P, grid, told, track ); }
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
		output << "  (" << total_tally[i].relUncer * 100.0;
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

// Set bin (or group structure), and calculate Chi group constants
void MGXS_Estimator::setBin( const std::string type, const std::vector<double> gr )
{
	// Set energy grid points
	grid = gr; 

	// Set scores = Flux, Capture, Fission, nuFission, Total, Scatter
	scores.push_back( std::make_shared<Flux_Score>()      );
	scores.push_back( std::make_shared<Capture_Score>()   );
	scores.push_back( std::make_shared<Fission_Score>()   );
	//scores.push_back( std::make_shared<nuFission_Score>() );
	scores.push_back( std::make_shared<Total_Score>()     );
	scores.push_back( std::make_shared<Scatter_Score>()   );
	Nscore = scores.size();

	// Set energy bin for scores
	std::vector<Tally_t> Tvec; // Vector of tallies corresponding to each score
	Tally_t T;
	Tvec.resize( Nscore, T );
	
	bin = std::make_shared<Energy_Bin> ( grid, Tvec, scores );
	Nbin = grid.size() - 1.0;

	// Set vector of energy bins scoring scattering only --> scattering matrix
	std::vector<std::shared_ptr<Score_t>> temp_scores;          // Scattering score
	temp_scores.push_back( std::make_shared<Scatter_Score>() );
	
	Tally_t Tsingle; // Single tally
	Tvec.resize( 1.0, Tsingle );

	for ( int i = 0 ; i < Nbin ; i++ )
	{
		std::shared_ptr<Bin_t> temp_bin = std::make_shared<Energy_Bin> ( grid, Tvec, temp_scores ); // Bin of initial energy with only one tally for scattering score
		matrix_bin.push_back( temp_bin ); // Generate the matrix whose column vectors are the initial energy bins
	}

	// Calculate Chi group constants
	calculateChi();
}


// Score at events
void MGXS_Estimator::score( const Particle_t& P, const double told, const double track /*= 0.0*/ )
{
	// Flux, Capture, Fission, nuFission, Total, Scater Tallies
	bin->score( P, grid, told, track );
	
	// Scattering matrix Tallies
	// First, simulate scattering reaction to get post-scattering particle information
	Particle_t P_final = P;
	P.region()->simulate_scatter( P_final );
	
	// Then, find the final energy bin location
	const int loc = Binary_Search( P_final.energy(), grid );
	
	// Final energy bin location is set, now score into the initial energy bin location
	matrix_bin[loc]->score( P, grid, told, track );
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


// Calculate Chi group constant from the universal Chi table
void MGXS_Estimator::calculateChi()
{
	// Read watt spectrum text file
	std::vector<double> evChi;   // Energy grid 
	std::vector<double> probChi; // PDF
    	std::ifstream inputFile("Chi.txt");
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
	
	// Integrate Chi spectrum for each group (marching up scheme)
	double Elow = 0.0; // Lower bound at each step of numerical integration
	double ylow = 0.0; // Its value
	double Eup  = 0.0; // Upper bound at each step of numerical integration (T.B.D. at each step)
	double yup  = 0.0; // Its value
	int    idx  = 1;   // Index in Chi spectrum data to be compared to decide the upper bound
	
	for ( int i = 0 ; i < Nbin ; i++ )
	{
		double sum = 0.0;       // Numerical integration sum
		
		while ( Eup < grid[i+1] )
		{
			if ( evChi[idx] > grid[i+1] )
			{
				// Upper bound is energy group grid
				Eup = grid[i+1];
				yup = Linterpolate( Eup, Elow, evChi[idx], ylow, probChi[idx] );
			}
			else
			{
				// Upper bound is chi spectrum data grid
				Eup = evChi[idx];
				yup = probChi[idx];
				idx++;
			}
			sum += ( yup + ylow ) * ( Eup - Elow ) * 0.5;
			
			Elow = Eup;
			ylow = yup;
		}
		Chi.push_back( sum );

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
	
	// Note energy group number is reversed for convenience
	if ( Nbin != 0 )
	{
		output << "Homogenized MG Constants," << std::endl;

		output << "g\t";
		output << "upper(" + bin->unit + ")\tlower(" + bin->unit + ")\t";
		output << std::setw(12) << std::left << "Chi" << "\t"; 
		for ( int i = 1 ; i < Nscore ; i++ ) 
		{ 
			output << std::setw(12) << std::left << scores[i]->name() << "\t"; 
			output << std::setw(12) << std::left << "uncertainty\t"; 
		}
		output << "\n";
	
		output << "----\t" << "------------\t" << "------------\t" << "------------\t";
		for ( int i = 1 ; i < Nscore ; i++ ) 
		{ 
			output << "------------\t" << "------------\t"; 
		}
		output << "\n";
	
		for ( int j = 0 ; j < Nbin ; j++ )
		{
			output << j+1;
			output << "\t" << std::setw(12) << std::left << grid[Nbin-j];
		        output << "\t" << std::setw(12) << std::left << grid[Nbin-j-1]; 
		        output << "\t" << std::setw(12) << std::left << Chi[Nbin-j-1]; 
			output << "\t";

			for ( int i = 1 ; i < Nscore ; i++ )
			{
				output << std::setw(12) << bin->tally[Nbin-j-1][i].mean << "\t";
			        output << std::setw(12) << bin->tally[Nbin-j-1][i].meanUncer << "\t";
			}
			output << "\n";
		}
		
		output << "\nScattering matrix," << std::endl;

		output << "g\t";
		output << "upper(" + bin->unit + ")\tlower(" + bin->unit + ")\t";
		for ( int j = 0 ; j < Nbin ; j++ )
		{
			const std::string g_prime = std::to_string(j+1);
			output << std::setw(12) << std::left << "g -> " + g_prime << "\t"; 
			output << std::setw(12) << std::left << "uncertainty\t"; 
		}
		output << "\n";
	
		output << "----\t" << "------------\t" << "------------\t";
		for ( int j = 0 ; j < Nbin ; j++ )
		{ 
			output << "------------\t" << "------------\t"; 
		}
		output << "\n";
	
		for ( int j = 0 ; j < Nbin ; j++ )
		{
			output << j+1;
			output << "\t" << std::setw(12) << std::left << grid[Nbin-j];
		        output << "\t" << std::setw(12) << std::left << grid[Nbin-j-1]; 
			output << "\t";
			for ( int i = 0 ; i < Nbin ; i++ )
			{
				output << std::setw(12) << matrix_bin[Nbin-i-1]->tally[Nbin-j-1][0].mean << "\t"; 
				output << std::setw(12) << matrix_bin[Nbin-i-1]->tally[Nbin-j-1][0].meanUncer << "\t"; 
			}
			output << "\n";
		}
	}
}


// Miscellaneous Estimator
//////////////////////////

// Surface PMF estimator (surface counting)
void Surface_PMF_Estimator::score( const Particle_t& P, const double told, const double track /*= 0.0*/ )
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
