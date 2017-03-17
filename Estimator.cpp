#include <cmath>

#include "Estimator.h"
#include "Surface.h"
#include "Region.h"


// Surface current estimator //
// Score at events
void Surface_Current_Estimator::score( const double w, const double null /*= 0.0*/ ) 
{ tally_hist += w; }

// Report results
void Surface_Current_Estimator::report( const double tTime ) 
{
	stats();                           // Compute mean, variance, and uncertainty
	const double r = meanUncer / mean; // relative uncertaianty of mean

	std::cout<< std::endl;
	std::cout<< e_name << " report" << std::endl;
	for ( int i = 0 ; i < e_name.length()+7 ; i++ ) { std::cout<< "="; }
	std::cout<< std::endl;
	std::cout<< "Number of particles crossing surface ";
	for ( int i = 0 ; i < o_names.size() ; i++ ) { std::cout<< o_names[i] << " "; }
	std::cout<< ":" << std::endl;
	std::cout<< "Mean     = " << mean << "  +/-  " << meanUncer << "  (" << r * 100.0 << "%)" << std::endl;
	std::cout<< "Variance = " << var << std::endl;
	std::cout<< "[F.O.M.: " << 1.0 / ( r*r * tTime ) << "]" << std::endl;
}



// Region flux estimator //
// Score at events
void Region_Flux_Estimator::score( const double w, const double path_length /*= 0.0*/ )
{ tally_hist += w * path_length; }

// Report results
void Region_Flux_Estimator::report( const double tTime )
{
      	stats();
			
	double volume;                      // to hold region volume
	const double vol = R->volume();  // to hold dummy volume (if vol = 0 --> give note!)
	if ( vol == 0.0 ) { volume = 1.0; }
	else { volume == vol; }
	
	const double r = meanUncer / mean;        // relative uncertaianty of mean
	
	std::cout<< std::endl;
	std::cout<< e_name << " report" << std::endl;
	for ( int i = 0 ; i < e_name.length()+7 ; i++ ) { std::cout<< "="; }
	std::cout<< std::endl;
	std::cout << "Volume-averaged particle flux in region " << o_names[0] << ":" << std::endl;
	std::cout<< "Mean     = " << mean / volume << "  +/-  " << meanUncer / volume << "  (" << r * 100.0 << "%)" << std::endl;
	std::cout<< "Variance = " << var / volume / volume  << std::endl;
	if ( vol == 0.0 ) { std::cout<<"..Note that region volume is not given, value above is not volume-averaged.." << std::endl; }
	std::cout<< "[F.O.M.: " << 1.0 / ( r*r * tTime ) << "]" << std::endl;

	if (absorp)
	{
		std::cout<< std::endl;
		std::cout<< "Absorption rate in region " << o_names[0] << ":" << std::endl;
		std::cout<< mean * R->SigmaA() << "  +/-  " << meanUncer * R->SigmaA() << std::endl;
	}
	if (scatter)
	{
		std::cout<< std::endl;
		std::cout<< "Scattering rate in region " << o_names[0] << ":" << std::endl;
		std::cout<< mean * R->SigmaS() << "  +/-  " << meanUncer * R->SigmaS() << std::endl;
	}
	if (capture)
	{
		std::cout<< std::endl;
		std::cout<< "Capture rate in region " << o_names[0] << ":" << std::endl;
		std::cout<< mean * R->SigmaC() << "  +/-  " << meanUncer * R->SigmaC() << std::endl;
	}
	if (fission)
	{
		std::cout<< std::endl;
		std::cout<< "Fission rate in region " << o_names[0] << ":" << std::endl;
		std::cout<< mean * R->SigmaF() << "  +/-  " << meanUncer * R->SigmaF() << std::endl;
	}
}



// Surface PMF estimator
// Score at events
void Surface_PMF_Estimator::score( const double w, const double null /*= 0.0*/ ) 
{ tally_hist++;}

// Report results
// output.txt file providing the PMF is created (or overwritten)
void Surface_PMF_Estimator::report( const double tTime ) 
{
	normalize();                        // Normalize the recorded PMF
	stats();                            // Compute mean, variance and PMF statistical uncertainty
	std::ofstream file( "output.txt" ); // Create output.txt file
	
	// Printouts (Uncertainty is only printed out in the output.txt file)
	std::cout<<std::endl;
	std::cout << e_name << " report" << std::endl;
	for ( int i = 0 ; i < e_name.length()+7 ; i++ ) { std::cout<< "="; }
	std::cout<< std::endl;
	std::cout << "PMF of # of particle crossing surface ";
	     file << "PMF of # of particle crossing surface ";
	for ( int i = 0 ; i < o_names.size() ; i++ ) { std::cout<< o_names[i] << " "; }
	for ( int i = 0 ; i < o_names.size() ; i++ ) {      file<< o_names[i] << " "; }
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
