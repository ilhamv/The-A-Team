#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include "XSec.h"


TEST_CASE( "Constant cross section", "[Constant_XSec]" ) {

  double aXS = 5.0;
  Constant_XSec aConstant_XSec( aXS );

  SECTION ( "xs" ) {

    double aE = 1.0;
    REQUIRE( aConstant_XSec.xs( aE ) == aXS );

  }
}

TEST_CASE( "1/V cross section", "[OverV_XSec]" ) {

  double a = 5.0;
  double b = 0.7;
  OverV_XSec aOverV_XSec( a, b );

  SECTION ( "xs" ) {

    double aE = 1E6; //eV
    REQUIRE( aOverV_XSec.xs( aE ) == 5.7 );

  }
}
