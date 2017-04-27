#include <vector>       // vector
#include <iostream>     // cout
#include <cstring>      // strcmp
#include <memory>       // shared_ptr, make_shared
#include <stack>        // stack
#include <cmath>        // exp
#include <sstream>      // ostringstream

#include "VReduction.h" // Split_Roulette
#include "Const.h"      // MAX
#include "pugixml.hpp"
#include "Geometry.h"
#include "Particle.h"
#include "Distribution.h"
#include "Source.h"
#include "Nuclide.h"
#include "Material.h"
#include "Reaction.h"
#include "Estimator.h"
#include "Shannon_Entropy_Mesh.h"
#include "XSec.h"
#include "XMLparser.h"


int main()
{
	// Variable declarations
	std::string                                   simName;           // Simulation name
	unsigned long long                            nhist;             // Number of particle samples
	double                                        Ecut_off  = 0.0;   // Energy cut-off
	double                                        tcut_off  = MAX;   // Time cut-off
	bool 										  eigenvalue = 0;	 // toggle for fixed-source vs k-eigenvalue simulation
	double                                        k = 1;			 // k-eigenvalue; initial guess of 1
	int										      ncycles = 1;		 // total number of cycles to run for k-eigenvalue simulation
	int										      npassive = 0;      // number of cycles to skip for source convergence for k-eigenvalue simulation
	std::shared_ptr <K_Eigenvalue_Estimator>      k_est;		     // estimator to compute new k after each cycle, and sample mean of k after source has converged
	std::shared_ptr <Shannon_Entropy_Mesh> 		  shannon_mesh;      // Shannon entropy mesh for eigenfunction convergence
	unsigned long long                            trackTime = 0;     // "Computation time" (particle track) for variance reduction
	Source_Bank                                   Sbank;             // Source bank
	std::stack  < Particle_t >                    Pbank;             // Particle bank
	std::stack  < Particle_t >                    Fbank;			 // Fission bank to update source for k-eigenvalue simulations
	std::vector < std::shared_ptr<Surface_t>    > Surface;           // Surfaces
	std::vector < std::shared_ptr<Region_t>     > Region;            // Regions
	std::vector < std::shared_ptr<Nuclide_t>    > Nuclide;           // Nuclides
	std::vector < std::shared_ptr<Material_t>   > Material;          // Materials
	std::vector < std::shared_ptr<Estimator_t>  > Estimator;         // Estimators  	
	// User-defined distributions
	std::vector < std::shared_ptr<Distribution_t<double>> > double_distributions;
  	std::vector < std::shared_ptr<Distribution_t<int>>    > int_distributions;
  	std::vector < std::shared_ptr<Distribution_t<Point_t>>> point_distributions;
    double                                           transportMethod = 0.0;   // standard is 0, anything else is delta
	
	// XML input parser	
	XML_input
	( simName, nhist, Ecut_off, tcut_off, eigenvalue, ncycles, npassive, k_est, shannon_mesh, Sbank, Surface, Region, Nuclide,
      Material, Estimator, double_distributions, int_distributions, point_distributions, transportMethod );
	
	std::cout<<"\nSimulation setup done,\nNow running the simulation...\n\n";
	std::cout.flush();

	// Simulation loop
	if (transportMethod != 0) // delta tracking
    {   std::cout << "RUNNING DELTA TRACKING MODE!"<<std::endl;

		eigenvalue = 0;	// delta tracking overrides k-eigenvalue

        for ( unsigned int isample = 0 ; isample < nhist ; isample++ )
	    {
	    	std::shared_ptr <Source_t> Source = Sbank.getSource();
			Particle_t source_particle = Source->getSource();
			source_particle.searchRegion( Region );
			Pbank.push( source_particle );

	    	// History loop
	    	while ( !Pbank.empty() )
	    	{
	    		Particle_t                P = Pbank.top(); // Working particle
	    		std::shared_ptr<Region_t> R = P.region();  // Working region
	    		Pbank.pop();
	    		// Particle loop
	    		while ( P.alive() )
	    		{
	    			double XS_max = 0.0;
                    for (const auto& R : Region)
                    {
                        if (R->SigmaT(P.energy()) > XS_max)
                        {XS_max = R->SigmaT(P.energy());}

                    }
                    
        	    			
	    			// Determine collision distance with max cross section
	    			double dcol = - std::log( Urand() ) / XS_max;;
	    
	    		    R->moveParticle( P, dcol, eigenvalue, npassive );
                    
                    P.searchRegion(Region);
                    
	    			// test collision
	    			if (Urand() < (R->SigmaT(P.energy()))/XS_max)
	    			{
	    				// Sample the collision and tally if there is any region tally
	    				R->collision( P, eigenvalue, k, Pbank, Fbank, shannon_mesh );
	    				
	    			}
        
	    				// Splitting & Roulette Variance Reduction Technique
	    				// 	More important : split
	    				// 	Less important : roulette
	    				// 	Same importance: do nothing
	    				// 	Note: old working region has the previous region importance and will be updated
	    				Split_Roulette( R, P, Pbank );
                        //Point_t test_point;
                        //test_point = P.pos();
                        //std::cout<<test_point.x<<std::endl;
                        //std::cout<<R->importance()<<std::endl;
                        
                    
                
                    // Accumulate "computation time"
	    				trackTime++;
	    		
	    		// Cut-off working particle?
	    		if ( P.energy() < Ecut_off || P.time() > tcut_off  || R->importance() ==0) { P.kill();}

	    		} // Particle is dead, end of particle loop		
	    		// Transport next Particle in the bank
        
	    	} // Particle bank is empty, end of history loop

	    	// Estimator history closeout
	    	for ( auto& E : Estimator ) { E->endHistory(); }
	    	// Start next history

	    } // All histories are done, end of simulation loop
    }
    else // standard transport
    {
        // Cycle loop
		for ( unsigned int icycle = 0; icycle < ncycles; icycle++ )
		{
			//std::cout << "cycle: " << icycle << std::endl;

			// Start timing for progress updating later
			std::clock_t start = std::clock();

			// Simulation loop
			for ( unsigned int isample = 0 ; isample < nhist ; isample++ )
			{

				//std::cout << "history: " << isample << std::endl;

				std::shared_ptr <Source_t> Source = Sbank.getSource();
				Particle_t source_particle = Source->getSource();
				source_particle.searchRegion( Region );
				Pbank.push( source_particle );

				//std::cout << "GOT TO HISTORY!" << std::endl;

				// History loop
				while ( !Pbank.empty() )
				{
					Particle_t                P = Pbank.top(); // Working particle
					std::shared_ptr<Region_t> R = P.region();  // Working region
					Pbank.pop();

					//particle_counter++;
					//std::cout << particle_counter << std::endl;

					//std::cout << "GOT TO PARTICLE!" << std::endl;

					// Particle loop
					while ( P.alive() )
					{

						//std::cout << "should this print: " << P.alive() << std::endl;

						std::pair< std::shared_ptr< Surface_t >, double > SnD; // To hold nearest surface and its distance
				
						// Determine nearest surface and its distance
						SnD = R->surface_intersect( P );

						// Determine collision distance
						double dcol = R->collision_distance( P.energy() );

						//std::cout << "COLLISION OR SURFACE CROSS?" << std::endl;
				
						// Hit surface?
						if ( dcol > SnD.second )
						{	
							// Surface hit! Move particle to surface, tally if there is any Region Tally
							R->moveParticle( P, SnD.second, eigenvalue, npassive, icycle );

							//std::cout << "PARTICLE MOVED!" << std::endl;

							// Implement surface hit:
							// 	Reflect angle for reflective surface
							// 	Cross the surface (move epsilon distance)
							// 	Search new particle region for transmission surface
							// 	Tally if there is any Surface Tally
							// 	Note: particle weight and working region are not updated yet
							SnD.first->hit( P, Region, eigenvalue, npassive, icycle );

							//std::cout << "SURFACE CROSS!" << std::endl;

							// std::cout << P.region()->name() << std::endl;

							// Splitting & Roulette Variance Reduction Technique
							// 	More important : split
							// 	Less important : roulette
							// 	Same importance: do nothing
							// 	Note: old working region has the previous region importance and will be updated
							Split_Roulette( R, P, Pbank );
					
							// Accumulate "computation time"
							trackTime++;
						}
				
						// Collide!!
						else
						{
							// Move particle to collision site and sample the collision and tally if there is any region tally
							R->moveParticle( P, dcol, eigenvalue, npassive, icycle );

							//std::cout << "PARTICLE MOVED!" << std::endl;

							R->collision( P, eigenvalue, k, Pbank, Fbank, shannon_mesh );
					
							//std::cout << "COLLISION!" << std::endl;

							// Accumulate "computation time"
							trackTime++;
						}
			
						// Cut-off working particle?
						if ( P.energy() < Ecut_off || P.time() > tcut_off ) { P.kill();}

					} // Particle is dead, end of particle loop		
				// Transport next Particle in the bank

				} // Particle bank is empty, end of history loop

				// Estimator history closeout
				if ( eigenvalue )
				{
					if ( icycle >= npassive )
						for ( auto& E : Estimator ) { E->endHistory(); }
				}
				else
					for ( auto& E : Estimator ) { E->endHistory(); }

				// Print progress report to terminal
				if ( ( !eigenvalue & ( fmod( std::log10( isample + 1 ), 1 ) == 0 ) ) || ( isample + 1 == nhist ) )
				{
		  			double duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
		  			if ( duration != 0 )
					{
		    			// to print scientific notation
		    			double sci1 = ( isample + 1 ) / std::pow( 10, std::floor( std::log10( isample + 1 ) ) );
		    			double sci2 = std::floor( std::log10( isample + 1 ) );
		    			std::cout << "  " << sci1 << "E" << sci2 << " histories took " << duration << " seconds.";
		    			if ( isample + 1 == nhist ) {std::cout << " Finished on "; }
		    			else {std::cout << ".. should finish on "; }
		    			// predict time left
		    			double timeLeft = duration / (isample + 1) * ( nhist - ( isample + 1 ) );
		    			// print local time + time left
		    			time_t rawtime;
		    			struct tm * timeinfo;
		    			time(&rawtime);
		    			rawtime += timeLeft;
		    			timeinfo = localtime (&rawtime);
		    			std::cout << asctime(timeinfo);
		  			}
				}

				// Start next history

			} // All histories are done, end of simulation loop

			if ( eigenvalue )
			{
				k = k_est->new_k( npassive, icycle );					// get new estimate of k for next iteration
				std::shared_ptr <Source_t> Source = Sbank.getSource();  // Update source with fission bank for next iteration

				//std::cout << "old nhist: " << nhist << std::endl;

				nhist = Fbank.size();

				//std::cout << "new nhist: " << nhist << std::endl;

				Source->update( Fbank );

				//std::cout << "new fission bank size: " << Fbank.size() <<std::endl;

				std::cout << icycle + 1 << " cycles completed -- k = " << k << " -- Shannon Entropy = " << shannon_mesh->entropy( nhist ) << std::endl;
				shannon_mesh->clear();
				if ( icycle == ( npassive - 1 ) ) { std::cout << std::endl << "All passive cycles completed; estimators activated." << std::endl << std::endl; }
				if ( icycle >= npassive ) { k_est->endCycle(); }		// if source has converged, start sampling mean of k
			}

		} // All cycles are done, end of cycle loop
	}

	// Generate outputs
	std::ostringstream output;                       // Output text
	std::ofstream file( simName + " - output.txt" ); // .txt file

	// Header
	output << "\n";
       	for ( int i = 0 ; i < simName.length()+6 ; i++ ) { output << "="; }
	output << "\n";
	output << "== " + simName + " ==\n";
       	for ( int i = 0 ; i < simName.length()+6 ; i++ ) { output << "="; }
	output << "\n";
	output << "Number of histories: " << nhist << "\n";
	output << "Track time: " << trackTime << "\n";

	// Report tallies

	if ( eigenvalue && (transportMethod == 0) ) { Estimator.push_back( k_est ); }
	for ( auto& E : Estimator ) { E->report( output, trackTime ); } // TrackTime is passed for F.O.M.
	
	// Print on monitor and file
	std::cout<< output.str();
	     file<< output.str();

	return 0;
}
