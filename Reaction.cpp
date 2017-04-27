#include <stack> // stack

#include "Reaction.h"
#include "Particle.h"
#include "Distribution.h"


// Capture kills the working particle
void Capture_Reaction::sample( Particle_t& P, std::stack< Particle_t >& Pbank )
{ P.kill(); }


// Scatter the working particle
void Scatter_Reaction::sample( Particle_t& P, std::stack< Particle_t >& Pbank )
{ P.scatter( scatter_dist->sample(), A ); }


// Fission reaction sample
void Fission_Reaction::sample( Particle_t& P, std::stack< Particle_t >& Pbank )
{
	// create random number of secondaries from multiplicity distributon nu_dist and
	// push all but one of them into the Particle bank, and reset the top particle 
	// if no secondaries, kill the particle

	int n = nu_dist->sample( r_nu->xs( P.energy() ) ); // sampled multiplicity
    
	if ( n != 0 ) 
    {
        // bank all but last particle (skips if n = 1)
	    for ( int i = 0 ; i < n - 1 ; i++ )
        {
	            Particle_t p( P.pos(), isotropic.sample(), Chi_dist->sample( P.energy() ), P.time(), P.weight() );
        	    p.setRegion( P.region() );
            	Pbank.push( p );
        }

		// reset the top particle
		P.setDirection( isotropic.sample() );
		P.setEnergy( Chi_dist->sample( P.energy() ) );
	}
	else
	{ P.kill(); }
}

// Implicit fission reaction sample
void Implicit_Fission_Reaction::sample( Particle_t& P, std::stack< Particle_t >& Pbank )
{ P.kill(); }
   

// Photoelectric is like capture
void Photoelectric_Reaction::sample( Particle_t& P, std::stack< Particle_t >& Pbank )
{ P.kill(); }

// Compton scatter the working particle
void ComptonScatter_Reaction::sample( Particle_t& P, std::stack< Particle_t >& Pbank )
{  
  Point_t dir1 = P.dir();
  Point_t dir2 = KleinNishina.sample();
  double dir1mag = std::sqrt( dir1.x*dir1.x + dir1.y*dir1.y + dir1.z*dir1.z );
  double dir2mag = std::sqrt( dir2.x*dir2.x + dir2.y*dir2.y + dir2.z*dir2.z );
  double mu0 = ( dir1.x*dir2.x + dir1.y*dir2.y + dir1.z*dir2.z ) / ( dir1mag * dir2mag );
  P.setEnergy( P.energy() / (1 + P.energy() / ME * ( 1 - mu0 ) ) ); // formula 2.17 in Knoll 4th ed.
  P.setDirection( dir2 );
}

// Pair produce reaction sample
void PairProduction_Reaction::sample( Particle_t& P, std::stack< Particle_t >& Pbank )
{
// create a secondary and push it into the Particle bank; reset the top particle
    Point_t dir1 = isotropic.sample();
    Point_t dir2( -dir1.x, -dir1.y, -dir1.z ); // 180 degrees opposite
    Particle_t p( P.pos(), dir1, ME, P.time(), P.weight() ); // 511 keV, assume all KE is deposited
    std::string photon = "photon";
    p.type = photon;
    p.setRegion( P.region() );
    Pbank.push( p );

    // reset the top particle
    P.setDirection( dir2 );
    P.setEnergy( ME ); // 511 keV
}


// Check reaction type
bool Capture_Reaction::type( const std::string s ) { if( s.compare("capture") == 0 ) { return true; } return false; }
bool Scatter_Reaction::type( const std::string s ) { if( s.compare("scatter") == 0 ) { return true; } return false; }
bool Fission_Reaction::type( const std::string s ) { if( s.compare("fission") == 0 ) { return true; } return false; }
bool Photoelectric_Reaction::type( const std::string s ) { if( s.compare("photoelectric") == 0 ) { return true; } return false; }
bool ComptonScatter_Reaction::type( const std::string s ) { if( s.compare("compton_scatter") == 0 ) { return true; } return false; }
bool PairProduction_Reaction::type( const std::string s ) { if( s.compare("pair_production") == 0 ) { return true; } return false; }

