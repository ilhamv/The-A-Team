#define CATCH_CONFIG_MAIN // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"

#include <vector>       // vector
#include <iostream>     // cout
#include <cstring>      // strcmp
#include <memory>       // shared_ptr, make_shared
#include <cmath>        // exp
#include <sstream>      // ostringstream
#include <limits>
#include <string>

#include "Estimator.h"
#include "Particle.h"
#include "Point.h"
#include "Geometry.h"


TEST_CASE( "Pulse height", "[PulseHeight_Estimator]" ) {
    
    SECTION ( "Pulse height unit test" ) {
        
        //create a defined detector region - let's call its ptr regCell
        std::string aregionname = "a region name";
        std::shared_ptr<Region_t> regCell = std::make_shared<Region_t> ( aregionname, 1.0 );
        
        //create a pulse height tally estimator and assign the cell of the detector
        const std::string n = "name of pulse height";
        PulseHeight_Estimator est( n , regCell );
        
        //set a Bin grid from 0 to 10 MeV
        std::vector <double> grid;
        for(int i = 0; i < 11 ; i++){grid.push_back(i*1.0e6);}
        est.setBin("bin", grid);
        
        //create a set of random energy for particles
        //the nullptr ones are the particles that go outside of the detector
        //the regCell ones are those that go into the detector
        //create dummy point for position and direction distribution
        Point_t point;
        //history 1
        Particle_t p1(point, point, 9.4e6, 0.0, 1.0);
                    p1.setRegion( regCell );
        Particle_t p2(point, point, 4.5e6, 0.0, 1.0);
                    p2.setRegion( nullptr );
        //history 2
        Particle_t p3(point, point, 5.2e6, 0.0, 1.0);
                    p3.setRegion( nullptr );
        Particle_t p4(point, point, 6.7e6, 0.0, 1.0);
                    p4.setRegion( regCell );
        //history 3
        Particle_t p5(point, point, 2.9e6, 0.0, 1.0);
                    p5.setRegion( regCell );
        Particle_t p6(point, point, 4.3e6, 0.0, 1.0);
                    p6.setRegion( nullptr );
        Particle_t p7(point, point, 6.4e6, 0.0, 1.0);
                    p7.setRegion( regCell );
        Particle_t p8(point, point, 2.1e6, 0.0, 1.0);
                    p8.setRegion( regCell );
        
        //score the particles for each history
        //history 1
        est.score( p1, 0.0, 0.0);
        est.score( p2, 0.0, 0.0);
        est.endHistory();
        //history 2
        est.score( p3, 0.0, 0.0);
        est.score( p4, 0.0, 0.0);
        est.endHistory();
        //hstory 3
        est.score( p5, 0.0, 0.0);
        est.score( p6, 0.0, 0.0);
        est.score( p7, 0.0, 0.0);
        est.score( p8, 0.0, 0.0);
        est.endHistory();
        
        //set the expected result
        std::vector <double> expectedResult;
        for(int j = 0 ; j < 10 ; j++ ){
            expectedResult.push_back(0.0);
        }
        expectedResult[4] = 1.0;
        expectedResult[1] = 1.0;
        expectedResult[7] = 1.0;

        //get the result from the code
        std::vector <double> codeResult = est.getTallySumVec();
        
        //then compare the results
        for (int k = 0; k < 10; k++ ){
            REQUIRE(codeResult[k] == Approx(  expectedResult[k]  )  );
        }
        
    }
    
}
