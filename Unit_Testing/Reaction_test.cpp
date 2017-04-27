#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include <cmath>
#include <vector>  // vector
#include <memory>  // shared_ptr
#include <stack>   // stack
#include <cstring> // string

#include "Reaction.h"
#include "Particle.h"
#include "Distribution.h"
#include "Point.h"
#include "XSec.h"

TEST_CASE( "Reaction_t method", "[Reaction_t]" ) {

  SECTION ( "xs" ) {
    Constant_XSec aConstant_XSec( 5.0 );
    std::shared_ptr< Constant_XSec > xspointer = std::make_shared< Constant_XSec >( aConstant_XSec );
    Capture_Reaction aCapture_Reaction( xspointer );
    REQUIRE( aCapture_Reaction.xs( 1E6 ) == 5.0 );
  }

}

TEST_CASE( "capture reaction", "[Capture_Reaction]" ) {

  SECTION ( "sample" ) {
    Point_t aposition( 1.0, 2.0, -5.0 );
    Point_t adirection( -1.0, -2.0, 5.0 );
    adirection.normalize();
    Particle_t aParticle( aposition, adirection );
    std::stack< Particle_t > aBank;
    Constant_XSec aConstant_XSec( 5.0 );
    std::shared_ptr< Constant_XSec > xspointer = std::make_shared< Constant_XSec >( aConstant_XSec );
    Capture_Reaction aCapture_Reaction( xspointer );
    aCapture_Reaction.sample( aParticle, aBank );
    REQUIRE( !aParticle.alive() );
  }

  SECTION ( "type" ) {
    Constant_XSec aConstant_XSec( 5.0 );
    std::shared_ptr< Constant_XSec > xspointer = std::make_shared< Constant_XSec >( aConstant_XSec );
    Capture_Reaction aCapture_Reaction( xspointer );
    REQUIRE( aCapture_Reaction.type( "capture" ) );
    REQUIRE( !aCapture_Reaction.type( "scatter" ) );
    REQUIRE( !aCapture_Reaction.type( "fission" ) );
  }

}

TEST_CASE( "scatter reaction", "[Scatter_Reaction]" ) {

  SECTION ( "sample" ) {
    Point_t aposition( 1.0, 2.0, -5.0 );
    Point_t adirection( -1.0, -2.0, 5.0 );
    adirection.normalize();
    Particle_t aParticle( aposition, adirection );
    std::stack< Particle_t > aBank;
    Constant_XSec aConstant_XSec( 5.0 );
    std::shared_ptr< Constant_XSec > xspointer = std::make_shared< Constant_XSec >( aConstant_XSec );
    Capture_Reaction aCapture_Reaction( xspointer );
    aCapture_Reaction.sample( aParticle, aBank );
    REQUIRE( !aParticle.alive() );
  }

  SECTION ( "type" ) {
    Constant_XSec aConstant_XSec( 5.0 );
    std::shared_ptr< Constant_XSec > xspointer = std::make_shared< Constant_XSec >( aConstant_XSec );
    Capture_Reaction aCapture_Reaction( xspointer );
    REQUIRE( aCapture_Reaction.type( "capture" ) );
    REQUIRE( !aCapture_Reaction.type( "scatter" ) );
    REQUIRE( !aCapture_Reaction.type( "fission" ) );
  }

}

TEST_CASE( "fission reaction", "[Fission_Reaction]" ) {

  SECTION ( "sample" ) {
    Point_t aposition( 1.0, 2.0, -5.0 );
    Point_t adirection( -1.0, -2.0, 5.0 );
    adirection.normalize();
    Particle_t aParticle( aposition, adirection );
    std::stack< Particle_t > aBank;
    Constant_XSec aConstant_XSec( 5.0 );
    std::shared_ptr< Constant_XSec > xspointer = std::make_shared< Constant_XSec >( aConstant_XSec );
    Capture_Reaction aCapture_Reaction( xspointer );
    aCapture_Reaction.sample( aParticle, aBank );
    REQUIRE( !aParticle.alive() );
  }

  SECTION ( "type" ) {
    Constant_XSec aConstant_XSec( 5.0 );
    std::shared_ptr< Constant_XSec > xspointer = std::make_shared< Constant_XSec >( aConstant_XSec );
    Capture_Reaction aCapture_Reaction( xspointer );
    REQUIRE( aCapture_Reaction.type( "capture" ) );
    REQUIRE( !aCapture_Reaction.type( "scatter" ) );
    REQUIRE( !aCapture_Reaction.type( "fission" ) );
  }

}
