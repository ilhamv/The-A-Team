#include <vector>       // vector
#include <iostream>     // cout
#include <cstring>      // strcmp
#include <memory>       // shared_ptr, make_shared
#include <stack>        // stack
#include <cmath>        // exp

#include "VReduction.h" // Split_Roulette
#include "pugixml.hpp"
#include "Surface.h"
#include "Particle.h"
#include "Distribution.h"
#include "Region.h"
#include "Source.h"
#include "Nuclide.h"
#include "Material.h"
#include "Reaction.h"
#include "Estimator.h"


// Function that returns an item from a vector of objects of type T by name provided
// the object class has a string and a method called name() allowing for it to be returned
template< typename T >
std::shared_ptr< T > findByName( const std::vector< std::shared_ptr< T > >& vec, const std::string name ) 
{
	for ( auto v : vec ) 
	{
		if ( v->name() == name ) { return v; }
	}
	return nullptr;
}


int main()
{
	// Variable declarations
	std::string                                   simName;    // Simulation name
	unsigned long long                            nhist;      // Number of particle samples
	unsigned long long                            trackTime;  // "Computation time" (particle track) for variance reduction
	Source_Bank                                   Sbank;      // Source bank
	std::stack  < Particle_t >                    Pbank;      // Particle bank
	std::vector < std::shared_ptr<Surface_t>    > Surface;    // Surfaces
	std::vector < std::shared_ptr<Region_t>     > Region;     // Regions
	std::vector < std::shared_ptr<Nuclide_t>    > Nuclide;    // Nuclides
	std::vector < std::shared_ptr<Material_t>   > Material;   // Materials
	std::vector < std::shared_ptr<Estimator_t>  > Estimator;  // Estimators  	
	// User-defined distributions
	std::vector < std::shared_ptr<Distribution_t<double>> > double_distributions;
  	std::vector < std::shared_ptr<Distribution_t<int>>    > int_distributions;
  	std::vector < std::shared_ptr<Distribution_t<Point_t>>> point_distributions;
	
	
	// XML input treatment //
	
	// User enters the XML file name and pugixml will attempt to load
	std::string input_file_name;
  	std::cout << "\nEnter XML input file name:\n";
  	std::cin  >> input_file_name;
	
	// XML input file
	pugi::xml_document input_file;
	pugi::xml_parse_result load_result = input_file.load_file( input_file_name.c_str() );

	// Check to see if result failed and throw an exception if it did
	if ( ! load_result ) 
	{
		std::cout << load_result.description() << std::endl;
		throw;
	}
	
	// Set simulation description: name and # of histories
  	pugi::xml_node input_simulation = input_file.child("simulation");
	simName = input_simulation.first_child().attribute("name").value();          // simulation name
	nhist   = input_simulation.first_child().attribute("histories").as_double(); // # of histories
	
	// Set user distributuions
  	pugi::xml_node input_distributions = input_file.child("distributions");

  	// Find total number of distributions
  	unsigned int num_distributions = 0;
  	for ( const auto& d : input_distributions ) { num_distributions++; }

  	// Since distributions may depend on other distributions, need to iterate
  	unsigned int set_distributions = 0;
  	while ( set_distributions < num_distributions ) 
	{
    		int previous_set_distributions = set_distributions;

    		for ( const auto& d : input_distributions ) 
		{
      			const std::string type = d.name();
      			const std::string name = d.attribute("name").value();
      			const std::string data = d.attribute("datatype").value();

      			// Double
			if ( data == "double" ) 
			{
        			// Skip rest of loop if distribution already done
        			if ( findByName( double_distributions, name ) ) { continue; }
	        		std::shared_ptr<Distribution_t<double>> Dist;
        			
				// Delta-double
				if ( type == "delta" ) 
				{
          				const double a = d.attribute("a").as_double();
          				Dist = std::make_shared< Delta_Distribution< double > > ( a, name );
	        		}

        			// Uniform-double
				else if ( type == "uniform" ) 
				{
          				const double a = d.attribute("a").as_double();
          				const double b = d.attribute("b").as_double();
	          			Dist = std::make_shared< Uniform_Distribution > ( a, b, name );
        			}
        			
				// Linear-double
				else if ( type == "linear" ) 
				{
          				const double a  = d.attribute("a").as_double();
	          			const double b  = d.attribute("b").as_double();
        	  			const double fa = d.attribute("fa").as_double();
          				const double fb = d.attribute("fb").as_double();
          				Dist = std::make_shared< Linear_Distribution > ( a, b, fa, fb, name );
        			}
	        		
				// Unknown
				else 
				{
          				std::cout << "unsupported distribution with data type " << data << std::endl;
          				throw;
        			}
				
				// Push new double-distribution
	        		double_distributions.push_back( Dist );
      			}

			// 3D point
			else if ( data == "point" ) 
			{
        			// Skip rest of loop if distribution already done
        			if ( findByName( point_distributions, name ) ) { continue; }
        			std::shared_ptr< Distribution_t< Point_t > > Dist;
		        	
				// Delta-point
				if ( type == "delta" ) 
				{
	          			const double x = d.attribute("x").as_double(); 
	          			const double y = d.attribute("y").as_double(); 
        	  			const double z = d.attribute("z").as_double();         
	          			Dist = std::make_shared< Delta_Distribution< Point_t > > ( Point_t( x, y, z ), name );
	        		}
        			
				// Isotropic-point
				else if ( type == "isotropic" ) 
				{
          				Dist = std::make_shared< IsotropicDirection_Distribution > ( name );
        			}
        			
				// Anisotropic-point
				else if ( type == "anisotropic" ) 
				{
          				const double u = d.attribute("u").as_double(); 
          				const double v = d.attribute("v").as_double(); 
          				const double w = d.attribute("w").as_double();        
          				std::shared_ptr< Distribution_t<double> > angDist = 
            				findByName( double_distributions, d.attribute("distribution").value() );
      
          				// in the angular distribution does not yet, skip to the end of the loop
          				if ( ! angDist ) { continue; }

          				Dist = std::make_shared< AnisotropicDirection_Distribution > ( Point_t( u, v, w ), angDist, name );
        			}

				// XYZ-point
        			else if ( type == "independentXYZ" ) 
				{
		          		std::shared_ptr< Distribution_t<double> > distX = findByName( double_distributions, d.attribute("x").value() ); 
        		  		std::shared_ptr< Distribution_t<double> > distY = findByName( double_distributions, d.attribute("y").value() ); 
          				std::shared_ptr< Distribution_t<double> > distZ = findByName( double_distributions, d.attribute("z").value() ); 

          				// if any of these distributions have not yet been resolved, skip to the end of the loop
	          			if ( !distX || !distY || !distZ ) { continue; }

		          		Dist = std::make_shared< IndependentXYZ_Distribution > ( distX, distY, distZ, name );
		        	}
		        	
				// Unknown
				else 
				{
          				std::cout << "unsupported " << data << " distribution of type " << type << std::endl;
		          		throw;
		        	}
				
				// Push new point-distribution
        			point_distributions.push_back( Dist );
      			}
      			
			// Unknown datatype
			else 
			{
        			std::cout << "unsupported distribution with data type " << data << std::endl;
        			throw;
      			}
      			
			// if we reach here, assume distribution has been set
      			set_distributions++;
		}
    		
		// check to see if number of distributions has increased, if not, caught in an infinite loop
    		if ( previous_set_distributions == set_distributions ) 
		{ 
      			std::cout << "distributions could not be resolved. " << std::endl;
      			throw;
    		}
  	}	

	// Set nuclides
	pugi::xml_node input_nuclides = input_file.child("nuclides");
  	for ( const auto& n : input_nuclides )
	{
    		const std::string          name = n.attribute("name").value();
		std::shared_ptr<Nuclide_t> Nuc  = std::make_shared<Nuclide_t> ( name );

    		// Add nuclide reactions
    		for ( const auto& r : n.children() ) 
		{
      			const std::string           rxn_type = r.name();
      			const double                xs       = r.attribute("xs").as_double();
      			
			// Capture
			if ( rxn_type == "capture" )
			{
        			Nuc->addReaction( std::make_shared<Capture_Reaction> ( xs ) );
      			}      
			
			// Scatter
			else if ( rxn_type == "scatter" )
			{
				if ( !r.attribute("distribution") ) 
				{ 
					std::cout << "distribution is required for scattering reaction" << std::endl;
					throw;
				}
				
				const std::string dist_name = r.attribute("distribution").value();
        			
				// Isotropic
				if ( dist_name == "isotropic" ) 
				{
					Nuc->addReaction( std::make_shared< Scatter_Reaction > ( xs, std::make_shared< IsotropicScatter_Distribution > () ) );
        			}
				
				// Henyey-Greenstein
				else if ( dist_name == "henyey-greenstein" ) 
				{
					if ( !r.attribute("g") ) 
					{ 
						std::cout << "parameter g is required for henyey-greenstein scattering distribution" << std::endl;
						throw;
					}
					const double g = r.attribute("g").as_double();
					Nuc->addReaction( std::make_shared< Scatter_Reaction > ( xs, std::make_shared< HGScatter_Distribution > ( g ) ) );
        			}
				
				// Linearly anisotropic
				else if ( dist_name == "linear" )
				{
					const double mubar = r.attribute("mubar").as_double();
					Nuc->addReaction( std::make_shared< Scatter_Reaction > ( xs, std::make_shared< HGScatter_Distribution > ( mubar ) ) );
        			}
        			
				// Unknown scattering distribution
				else 
				{
          				std::cout << "unknown scattering distribution " << dist_name << " in nuclide " << name << std::endl;
          				throw;
        			}
      			}

			// Fission
			else if ( rxn_type == "fission" )
			{
				if ( !r.attribute("multiplicity") ) 
				{ 
					std::cout << "multiplicity is required for fission reaction" << std::endl;
					throw;
				}

				const std::string mult_dist_name = r.attribute("multiplicity").value();
				
				// Average
				if ( mult_dist_name == "average" )
				{
					if ( !r.attribute("nubar") ) 
					{ 
						std::cout << "parameter nubar is required for average fission multiplicity" << std::endl;
						throw;
					}
					const double nubar = r.attribute("nubar").as_double();
					Nuc->addReaction( std::make_shared< Fission_Reaction > ( xs, std::make_shared< Average_Multiplicity_Distribution > ( nubar ) ) );					
				}

				// Terrel
				else if ( mult_dist_name == "terrel" )
				{
					if ( !r.attribute("nubar") || !r.attribute("gamma") || !r.attribute("b") || !r.attribute("nmax") ) 
					{ 
						std::cout << "parameter nubar, gamma, b, and nmax are required for terrel fission multiplicity" << std::endl;
						throw;
					}

					const double nubar = r.attribute("nubar").as_double();
					const double gamma = r.attribute("gamma").as_double();
					const double b     = r.attribute("b").as_double();
					const int    nmax  = r.attribute("nmax").as_int();
					const std::vector< std::pair< int, double > > v;     // a dummy, as it is required for discrete distribution base class
					Nuc->addReaction( std::make_shared< Fission_Reaction > ( xs, std::make_shared< Terrel_Multiplicity_Distribution > ( nubar, gamma, b, nmax, v ) ) );
				}
				
				// Unknown multiplicity distribution
				else 
				{
          				std::cout << "unknown fission multiplicity distribution " << mult_dist_name << " in nuclide " << name << std::endl;
          				throw;
        			}
			}
      			
			// Unknown reaction
			else 
			{
        			std::cout << "unknown reaction type " << rxn_type << std::endl;
        			throw;
      			}
    		}
		
		// Push new nuclide
		Nuclide.push_back( Nuc );
 	}

  	// Set materials
  	pugi::xml_node input_materials = input_file.child("materials");
  	for ( const auto& m : input_materials )
	{
    		const           std::string name = m.attribute("name").value();
    		std::shared_ptr<Material_t>  Mat = std::make_shared<Material_t> ( name );

    		// Add material nuclides
    		for ( const auto& n : m.children() )
		{
			if ( (std::string) n.name() == "nuclide" ) 
			{
        			const std::string                   nuclide_name = n.attribute("name").value();
        			const double                        density      = n.attribute("density").as_double();
				const std::shared_ptr<Nuclide_t>    nucPtr       = findByName( Nuclide, nuclide_name );
				
      				if ( nucPtr ) 
				{
        				Mat->addNuclide( nucPtr, density );
      				}
      				else
		       		{
        				std::cout << "unknown nuclide " << nuclide_name << " in material " << name << std::endl;
        				throw;
      				}
      			}
    		}
    		
		// Push new material
		Material.push_back( Mat );
  	}
  
	// Set surfaces
  	pugi::xml_node input_surfaces = input_file.child("surfaces");
  	for ( const auto& s : input_surfaces )
	{
    		const std::string type = s.name();
      		const std::string name = s.attribute("name").value();
		const std::string bc   = s.attribute("bc").value();
    		std::shared_ptr< Surface_t > S;
    		
		// Plane
		if ( type == "plane" ) 
		{
      			const double a = s.attribute("a").as_double();
      			const double b = s.attribute("b").as_double();
      			const double c = s.attribute("c").as_double();
      			const double d = s.attribute("d").as_double();
      			if ( bc == "reflect" )
			{
				S = std::make_shared< Plane_Reflective > ( name, a, b, c, d );
			}
			else
			{
				S = std::make_shared< Plane_Surface > ( name, a, b, c, d );
			}
    		}
    		
		// Sphere
		else if ( type == "sphere" )
		{
      			const double x = s.attribute("x").as_double();
      			const double y = s.attribute("y").as_double();
      			const double z = s.attribute("z").as_double();
      			const double r  = s.attribute("r").as_double();
      			S = std::make_shared< Sphere_Surface > ( name, x, y, z, r );
		}
		
		// Cylinder-x
		else if ( type == "cylinder_x" )
		{
      			const double y = s.attribute("y").as_double();
      			const double z = s.attribute("z").as_double();
      			const double r  = s.attribute("r").as_double();
      			S = std::make_shared< CylinderX_Surface > ( name, y, z, r );
		}
		
		// Cylinder-z
		else if ( type == "cylinder_z" )
		{
      			const double x = s.attribute("x").as_double();
      			const double y = s.attribute("y").as_double();
      			const double r  = s.attribute("r").as_double();
      			S = std::make_shared< CylinderZ_Surface > ( name, x, y, r );
		}
		
		// Unknown surface typw
		else 
		{
      			std::cout << " unkown surface type " << type << std::endl;
      			throw;
    		}
    	
		// Push new surface
		Surface.push_back( S );
  	}
  
	// Set regions
  	pugi::xml_node input_regions = input_file.child("regions");
  	for ( const auto& r : input_regions ) 
	{
    		const std::string name       = r.attribute("name").value();
		double            importance = 1.0; // default
		double            volume     = 0.0; // default
    		std::shared_ptr<Region_t> Reg;

    		// Modify region importance
    		if ( r.attribute("importance") ) 
    		{
	    		importance = r.attribute("importance").as_double();
    		}
    		
		// Modify region volume
    		if ( r.attribute("volume") ) 
    		{
	    		volume = r.attribute("volume").as_double();
    		}

    		// Set region material
    		if ( r.attribute("material") ) 
		{
      			const std::shared_ptr<Material_t> matPtr = findByName( Material, r.attribute("material").value() );
      			if ( matPtr ) 
			{
    				Reg  = std::make_shared<Region_t> ( name, importance, volume );
        			Reg->setMaterial( matPtr );
      			}
			else
		       	{
        			std::cout << "unknown material " << r.attribute("material").value() << " in region " << name << std::endl;
        			throw;
      			} 
   		}
		
		// Vacuum (no material)
		else
		{
    			Reg  = std::make_shared<Region_Vacuum> ( name, importance, volume );
		}

   
    		// Set region bounding surfaces
    		for ( const auto& s : r.children() ) 
		{
      			if ( (std::string) s.name() == "surface" ) 
			{
        			std::string name  = s.attribute("name").value();
        			const int   sense = s.attribute("sense").as_int();

        			std::shared_ptr<Surface_t> SurfPtr = findByName( Surface, name );
        			
				if ( SurfPtr ) 
				{
          				Reg->addSurface( findByName( Surface, name ), sense );
        			}
				else 
				{
          				std::cout << "unknown surface with name " << name << std::endl;
          				throw;
        			}
      			}
			else 
			{
        			std::cout << "unknown data type " << s.name() << " in region " << name << std::endl;
        			throw;
      			}
    		} 
    		
		// Push new region
		Region.push_back( Reg );
  	}
  	
	// Set estimators
  	pugi::xml_node input_estimators = input_file.child("estimators");
  	for ( const auto& e : input_estimators )
	{
    		const std::string type = e.name();
    		const std::string name = e.attribute("name").value();
    		std::shared_ptr<Estimator_t> Est;
    		
		// Surface current estimator
		if ( type == "surface_current" ) 
		{
      			Est = std::make_shared<Surface_Current_Estimator> ( name );

      			// Get the surfaces
      			for ( const auto& s : e.children() )
			{
        			if ( (std::string) s.name() == "surface" )
				{
          				const std::string                name    = s.attribute("name").value();
          				const std::shared_ptr<Surface_t> SurfPtr = findByName( Surface, name );
          				
					// Add estimator to surface
					if ( SurfPtr ) 
					{
            					SurfPtr->addEstimator( Est );
						Est->addName( SurfPtr->name() );
          				}
          				else 
					{
            					std::cout << "unknown surface label " << name << " in estimator " << e.attribute("name").value() << std::endl;
						throw;
          				}
        			}
      			}
    		}
    		
		// Counting surface estimator (results Probability Mass Function)
		else if ( type == "surface_counting" ) 
		{
      			Est = std::make_shared<Surface_PMF_Estimator > ( name );

      			// Get the surfaces
      			for ( const auto& s : e.children() ) 
			{
        			if ( (std::string) s.name() == "surface" ) 
				{
          				const std::string                 name    = s.attribute("name").value();
          				const std::shared_ptr<Surface_t> SurfPtr = findByName( Surface, name );
          				
					// Add estimator to surface and vice versa
					if ( SurfPtr ) 
					{
            					SurfPtr->addEstimator( Est );
						Est->addName( SurfPtr->name() );
          				}
          				else 
					{
            					std::cout << "unknown surface label " << name << " in estimator " << e.attribute("name").value() << std::endl;
						throw;
          				}
        			}
      			}			 
    		}
		
		// Region flux estimator
		else if ( type == "region_flux" ) 
		{
      			bool scatter    = false; // defaults
      			bool capture    = false;
      			bool fission    = false;
      			bool absorption = false;
			
			// Modify reaction rate flags
    			if ( e.attribute("scatter") ) 
    			{
	    			scatter = e.attribute("scatter").as_bool();
    			}
    			if ( e.attribute("capture") ) 
    			{
	    			capture = e.attribute("capture").as_bool();
    			}
    			if ( e.attribute("fission") ) 
    			{
	    			fission = e.attribute("fission").as_bool();
    			}
    			if ( e.attribute("absorption") ) 
    			{
	    			absorption = e.attribute("absorption").as_bool();
    			}
			
			// Get the region
      			for ( const auto& r : e.children() )
			{
        			if ( (std::string) r.name() == "region" )
				{
          				const std::string          r_name  = r.attribute("name").value();
          				std::shared_ptr<Region_t>  RegPtr  = findByName( Region, r_name );
          				
					// Add estimator to the region
					if ( RegPtr ) 
					{
						Est = std::make_shared<Region_Flux_Estimator> ( name, RegPtr, scatter, capture, fission, absorption );
            					RegPtr->addEstimator( Est );
						Est->addName( RegPtr->name() );
          				}
          				else 
					{
            					std::cout << "unknown region label " << r_name << " in estimator " << name << std::endl;
						throw;
          				}
        			}
      			}
    		}
    		
		// Unknown estimator type
		else 
		{
      			std::cout << "unknown estimator type " << name << std::endl;
      			throw;
    		}
    
		// Push new estimator
    		Estimator.push_back( Est );
  	}

	// Set source bank
  	pugi::xml_node input_sources = input_file.child("sources");
  	for ( const auto& s : input_sources )
	{
		double prob = 1.0; // default
    		std::shared_ptr<Source_t> Src;

		// Modify source probability
    		if ( s.attribute("probability") )
    		{
	    		prob = s.attribute("probability").as_double();
    		}

		// 1D normal plane-x source
		if ( (std::string) s.name() == "normal_plane_x" )
		{
			const double x     = s.attribute("x").as_double();
			const int    sense = s.attribute("sense").as_double();
			Src = std::make_shared<Normal_HPlanex_Source> ( x, sense );
    		}
		
		// 1D isotropic plane-x source
		else if ( (std::string) s.name() == "isotropic_plane_x" )
		{
			const double x     = s.attribute("x").as_double();
			const int    sense = s.attribute("sense").as_double();
			Src = std::make_shared<Isotropic_HPlanex_Source> ( x, sense );
    		}
    		
		// 1D slab-x source
		else if ( (std::string) s.name() == "isotropic_slab_x" )
		{
			const double a = s.attribute("a").as_double();
			const double b = s.attribute("b").as_double();
			Src = std::make_shared<Isotropic_Slabx_Source> ( a, b );
    		}
		
		// Isotropic point
		else if ( (std::string) s.name() == "isotropic_point" )
		{
			const double x = s.attribute("x").as_double();
			const double y = s.attribute("y").as_double();
			const double z = s.attribute("z").as_double();
			Src = std::make_shared<Isotropic_Point_Source> ( x, y, z );
    		}
		
		// Isotropic spherical shell
		else if ( (std::string) s.name() == "sphere_shell" )
		{
			const double x  = s.attribute("x").as_double();
			const double y  = s.attribute("y").as_double();
			const double z  = s.attribute("z").as_double();
			const double ri = s.attribute("ri").as_double();
			const double ro = s.attribute("ro").as_double();
			Src = std::make_shared<Sphere_Shell_Source> ( x, y, z, ri, ro );
    		}

		// Normal disk-z
		else if ( (std::string) s.name() == "disk_z" )
		{
			const double x     = s.attribute("x").as_double();
			const double y     = s.attribute("y").as_double();
			const double z     = s.attribute("z").as_double();
			const double r     = s.attribute("r").as_double();
			const int    sense = s.attribute("sense").as_double();
			Src = std::make_shared<DiskZ_Source> ( x, y, z, r, sense );
    		}
		
		// User defined
		else if ( (std::string) s.name() == "usource" )
		{
		
  			std::string pos_dist_name = s.attribute("position").value();
  			std::string dir_dist_name = s.attribute("direction").value();

  			std::shared_ptr< Distribution_t<Point_t> > posDist = findByName( point_distributions, pos_dist_name );
  			std::shared_ptr< Distribution_t<Point_t> > dirDist = findByName( point_distributions, dir_dist_name );

            		if ( posDist && dirDist ) 
			{
				Src = std::make_shared<User_Source> ( posDist, dirDist );
  			}
  			else 
			{
    				if ( ! posDist ) { std::cout << " unknown position distribution "  << pos_dist_name << " in user-defined source " << std::endl; }
    				if ( ! dirDist ) { std::cout << " unknown direction distribution " << dir_dist_name << " in user-defined source " << std::endl; }
    				throw;
  			}
		}

		// Unknown source type
		else 
		{
            		std::cout << "unknown source type " << (std::string) s.name() << std::endl;
			throw;
          	}

		// Push new source
		Sbank.addSource( Src, prob );
  	}

	// XML input treatment done //
	
	std::cout<<"\nSimulation setup done,\nNow running the simulation...\n\n";
	std::cout.flush();


	// Simulation loop
	for ( unsigned int isample = 0 ; isample < nhist ; isample++ )
	{
		Pbank.push( Sbank.getSource( Region ) );

		// History loop
		while ( !Pbank.empty() )
		{
			Particle_t                P = Pbank.top(); // Working particle
			std::shared_ptr<Region_t> R = P.region();  // Working region
			Pbank.pop();
			
			// Particle loop
			while ( P.alive() )
			{
				std::pair< std::shared_ptr< Surface_t >, double > SnD; // To hold nearest surface and its distance
				
				// determine nearest surface and its distance
				SnD = R->surface_intersect( P );

				// determine collision distance
				double dcol = R->collision_distance();
				
				// hit surface?
				if ( dcol > SnD.second )
				{	
					// Surface hit! Move to surface, tally if there is any Region Tally
					R->moveParticle( P, SnD.second );
					
					// Accumulate "computation time" for variance reduction
					trackTime++;

					// Implement surface hit:
					//  non-Reflective -> move epsilon distance, search new region
					// 	Reflective -> reflect direction, move epsilon distance
					// 	tally if there is any Surface Tally
					// 	particle weight and working region are not updated yet
					SnD.first->hit( P, Region );

					// Splitting & Roulette Variance Reduction Technique
					// More important : split
					// Less important : roulette
					// Same importance: do nothing
					// old working region has the previous region importance and will be updated
					Split_Roulette( R, P, Pbank );
				}
					
				else
				{
					// Collide!! Move to collision site and sample the collision
					// tally if there is any region tally
					R->moveParticle( P, dcol );
					R->collision( P, Pbank );
					
					// Accumulate "computation time" for variance reduction
					trackTime++;
				}	

			} // Particle is dead, end of particle loop
			
			// Transport next Particle in the bank queue

		} // Particle bank is empty, end of history loop

		// Estimator history closeout
		for ( auto E : Estimator ) { E->endHistory(); }
		// Start next history

	} // All histories are done, end of simulation loop

	// Output header
	std::cout<< std::endl;
       	for ( int i = 0 ; i < simName.length()+6 ; i++ ) { std::cout<< "="; }
	std::cout<< std::endl;
	std::cout<< "== "<< simName << " ==" << std::endl;
       	for ( int i = 0 ; i < simName.length()+6 ; i++ ) { std::cout<< "="; }
	std::cout<< std::endl;
	std::cout<< "Number of histories: " << nhist <<std::endl;
	std::cout<< "Track time: " << trackTime << std::endl;

	// Report tallies
	for ( auto E : Estimator ) { E->report( trackTime ); }
	
	return 0;
}
