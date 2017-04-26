#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include "Geometry.h"

#include <memory>      // shared_ptr
#include <vector>      // vector
#include <cstring>     // string

#include "Const.h"     // EPSILON
#include "Particle.h"
#include "Estimator.h"
#include "Point.h"
#include "Material.h"

TEST_CASE( "testing geometry template and surface template methods", "[Geometry_t and Surface_t]" ) {

  SECTION ( "name" ) {
    std::string asurfacename = "a surface name";
    std::string aboundarycondition = "transmission";
    double alocation = 3.0;
    PlaneX_Surface aPlaneX_Surface( asurfacename, aboundarycondition, alocation );
    REQUIRE( aPlaneX_Surface.name() == asurfacename );
  }

  SECTION ( "addEstimator" ) {
/* not finished
  std::string asurfacename = "a surface name";
  std::string aboundarycondition = "transmission";
  double alocation = 3.0;
  PlaneX_Surface aPlaneX_Surface( asurfacename, aboundarycondition, alocation );
    Tally_t aTally;
    aPlaneX_Surface.addEstimator( aTally** );
    std::string aregionname = "a region name";
    double aimportance = 1.0;
    Region_t aRegion( aregionname, aimportance );
    aRegion.addSurface( aPlaneX_Surface*, 1 );
    Point_t aposition( -1.0, 0.0, 0.0 );
    Point_t adirection( 1.0, 0.0, 0.0 );
    Particle_t aParticle( aposition, adirection );
    double adistance = 100;
    aRegion.moveParticle( aParticle*, adistance);
    REQUIRE( false );
*/  }

  SECTION ( "cross" ) {
    // protected class
  }

  SECTION ( "hit and cross" ) {
    std::string asurfacename = "a surface name";
    std::string aboundarycondition = "transmission";
    double alocation = 3.0;
    PlaneX_Surface aPlaneX_Surface( asurfacename, aboundarycondition, alocation );
    Point_t aposition( 3.0, 0.0, 0.0);
    Point_t adirection( 1.0, 0.0, 0.0 );
    Particle_t aParticle( aposition, adirection );
    std::vector < std::shared_ptr<Region_t> > regions;
    std::string aregionname = "a region name";
    std::shared_ptr<Region_t> aRegion = std::make_shared<Region_t> ( aregionname, 1.0 );
    std::shared_ptr<PlaneX_Surface> surfacepointer = std::make_shared<PlaneX_Surface>( aPlaneX_Surface );
    aRegion->addSurface( surfacepointer, 1 );
    regions.push_back( aRegion );
    aPlaneX_Surface.hit( aParticle, regions );
    REQUIRE( aParticle.pos().x == Approx( 3.0+EPSILON ) );
    REQUIRE( aParticle.region() == aRegion );
  }

}

TEST_CASE( "infinite plane perpendicular to x axis", "[PlaneX_Surface]" ) {

  SECTION ( "reflect" ) {
    std::string asurfacename = "a surface name";
    std::string aboundarycondition = "reflective";
    double alocation = 3.0;
    PlaneX_Surface aPlaneX_Surface( asurfacename, aboundarycondition, alocation );
    Point_t aposition( 3.0, 0.0, 0.0);
    Point_t adirection( 1.0, -2.0, 30.0 );
    Particle_t aParticle( aposition, adirection );
    std::vector < std::shared_ptr<Region_t> > regions;
    std::string aregionname = "a region name";
    std::shared_ptr<Region_t> aRegion = std::make_shared<Region_t> ( aregionname, 1.0 );
    std::shared_ptr<PlaneX_Surface> surfacepointer = std::make_shared<PlaneX_Surface>( aPlaneX_Surface );
    aRegion->addSurface( surfacepointer, -1 );
    regions.push_back( aRegion );
    aPlaneX_Surface.hit( aParticle, regions );
    REQUIRE( aParticle.pos().x == Approx( 3.0-EPSILON ) );
    REQUIRE( aParticle.dir().x == Approx( -0.0332411 ) );
    REQUIRE( aParticle.dir().y == Approx( -0.0664822 ) );
    REQUIRE( aParticle.dir().z == Approx( 0.99723374 ) );
  }

  SECTION ( "eval" ) {
    std::string asurfacename = "a surface name";
    std::string aboundarycondition = "reflective";
    double alocation = 3.0;
    PlaneX_Surface aPlaneX_Surface( asurfacename, aboundarycondition, alocation );
    Point_t aposition( -7.0, -2.0, 30.0);
    REQUIRE( aPlaneX_Surface.eval(aposition) == -10.0 );
  }

  SECTION ( "distance" ) {
    std::string asurfacename = "a surface name";
    std::string aboundarycondition = "reflective";
    double alocation = 3.0;
    PlaneX_Surface aPlaneX_Surface( asurfacename, aboundarycondition, alocation );
    Point_t aposition( 0.0, 0.0, 0.0);
    Point_t adirection( 1.0, 1.0, 0.0 );
    Particle_t aParticle( aposition, adirection );
    REQUIRE( aPlaneX_Surface.distance(aParticle) == Approx( std::sqrt(18) ) );
  }

}

TEST_CASE( "infinite plane perpendicular to y axis", "[PlaneY_Surface]" ) {

  SECTION ( "reflect" ) {
    std::string asurfacename = "a surface name";
    std::string aboundarycondition = "reflective";
    double alocation = 3.0;
    PlaneY_Surface aPlaneY_Surface( asurfacename, aboundarycondition, alocation );
    Point_t aposition( 0.0, 3.0, 0.0);
    Point_t adirection( -1.0, 2.0, 30.0 );
    Particle_t aParticle( aposition, adirection );
    std::vector < std::shared_ptr<Region_t> > regions;
    std::string aregionname = "a region name";
    std::shared_ptr<Region_t> aRegion = std::make_shared<Region_t> ( aregionname, 1.0 );
    std::shared_ptr<PlaneY_Surface> surfacepointer = std::make_shared<PlaneY_Surface>( aPlaneY_Surface );
    aRegion->addSurface( surfacepointer, -1 );
    regions.push_back( aRegion );
    aPlaneY_Surface.hit( aParticle, regions );
    REQUIRE( aParticle.pos().y == Approx( 3.0-EPSILON ) );
    REQUIRE( aParticle.dir().x == Approx( -0.0332411 ) );
    REQUIRE( aParticle.dir().y == Approx( -0.0664822 ) );
    REQUIRE( aParticle.dir().z == Approx( 0.99723374 ) );
  }

  SECTION ( "eval" ) {
    std::string asurfacename = "a surface name";
    std::string aboundarycondition = "reflective";
    double alocation = 3.0;
    PlaneY_Surface aPlaneY_Surface( asurfacename, aboundarycondition, alocation );
    Point_t aposition( -7.0, -2.0, 30.0);
    REQUIRE( aPlaneY_Surface.eval(aposition) == -5.0 );
  }

  SECTION ( "distance" ) {
    std::string asurfacename = "a surface name";
    std::string aboundarycondition = "reflective";
    double alocation = 3.0;
    PlaneY_Surface aPlaneY_Surface( asurfacename, aboundarycondition, alocation );
    Point_t aposition( 0.0, 0.0, 0.0);
    Point_t adirection( 1.0, 1.0, 0.0 );
    Particle_t aParticle( aposition, adirection );
    REQUIRE( aPlaneY_Surface.distance(aParticle) == Approx( std::sqrt(18) ) );
  }

}

TEST_CASE( "infinite plane perpendicular to z axis", "[PlaneZ_Surface]" ) {

  SECTION ( "reflect" ) {
    std::string asurfacename = "a surface name";
    std::string aboundarycondition = "reflective";
    double alocation = 3.0;
    PlaneZ_Surface aPlaneZ_Surface( asurfacename, aboundarycondition, alocation );
    Point_t aposition( 0.0, 0.0, 3.0);
    Point_t adirection( -1.0, -2.0, 30.0 );
    Particle_t aParticle( aposition, adirection );
    std::vector < std::shared_ptr<Region_t> > regions;
    std::string aregionname = "a region name";
    std::shared_ptr<Region_t> aRegion = std::make_shared<Region_t> ( aregionname, 1.0 );
    std::shared_ptr<PlaneZ_Surface> surfacepointer = std::make_shared<PlaneZ_Surface>( aPlaneZ_Surface );
    aRegion->addSurface( surfacepointer, -1 );
    regions.push_back( aRegion );
    aPlaneZ_Surface.hit( aParticle, regions );
    REQUIRE( aParticle.pos().z == Approx( 3.0-EPSILON ) );
    REQUIRE( aParticle.dir().x == Approx( -0.0332411 ) );
    REQUIRE( aParticle.dir().y == Approx( -0.0664822 ) );
    REQUIRE( aParticle.dir().z == Approx( -0.99723374 ) );
  }

  SECTION ( "eval" ) {
    std::string asurfacename = "a surface name";
    std::string aboundarycondition = "reflective";
    double alocation = 3.0;
    PlaneZ_Surface aPlaneZ_Surface( asurfacename, aboundarycondition, alocation );
    Point_t aposition( -7.0, -2.0, 30.0);
    REQUIRE( aPlaneZ_Surface.eval(aposition) == 27.0 );
  }

  SECTION ( "distance" ) {
    std::string asurfacename = "a surface name";
    std::string aboundarycondition = "reflective";
    double alocation = 3.0;
    PlaneZ_Surface aPlaneZ_Surface( asurfacename, aboundarycondition, alocation );
    Point_t aposition( 0.0, 0.0, 0.0);
    Point_t adirection( 0.0, 1.0, 1.0 );
    Particle_t aParticle( aposition, adirection );
    REQUIRE( aPlaneZ_Surface.distance(aParticle) == Approx( std::sqrt(18) ) );
  }

}
/*
TEST_CASE( "infinite arbitrary plane", "[Plane_Surface]" ) {

  SECTION ( "reflect" ) {
    std::string asurfacename = "a surface name";
    std::string aboundarycondition = "reflective";
    Plane_Surface aPlane_Surface( asurfacename, aboundarycondition, 1.0, -2.0, 3.0, 0.5 );
    Point_t aposition( 0.5, 0.0, 0.0);
    Point_t adirection( -2.0, -1.0, 3.0 );
    adirection.normalize();
    Particle_t aParticle( aposition, adirection );
    std::vector < std::shared_ptr<Region_t> > regions;
    std::string aregionname = "a region name";
    std::shared_ptr<Region_t> aRegion = std::make_shared<Region_t> ( aregionname, 1.0 );
    std::shared_ptr<Plane_Surface> surfacepointer = std::make_shared<Plane_Surface>( aPlane_Surface );
    aRegion->addSurface( surfacepointer, -1 );
    regions.push_back( aRegion );
    aPlane_Surface.hit( aParticle, regions );
    //REQUIRE( aParticle.dir().x == 0 );
    //REQUIRE( aParticle.dir().y == 0 );
    //REQUIRE( aParticle.dir().z == 0 );
  }

  SECTION ( "eval" ) {
    std::string asurfacename = "a surface name";
    std::string aboundarycondition = "reflective";
    Plane_Surface aPlane_Surface( asurfacename, aboundarycondition, 1.0, -2.0, 3.0, 0.5 );
    Point_t aposition( -7.0, -2.0, 30.0);
    REQUIRE( aPlane_Surface.eval(aposition) == 86.5 );
  }

  SECTION ( "distance" ) {
    std::string asurfacename = "a surface name";
    std::string aboundarycondition = "reflective";
    Plane_Surface aPlane_Surface( asurfacename, aboundarycondition, 0.0, 0.0, 1.0, 3.0 );
    Point_t aposition( 0.0, 0.0, 0.0);
    Point_t adirection( 0.0, 1.0, 1.0 );
    Particle_t aParticle( aposition, adirection );
    REQUIRE( aPlane_Surface.distance(aParticle) == std::sqrt(18) );
  }

}

TEST_CASE( "arbitrary sphere", "[Sphere_Surface]" ) {

  SECTION ( "" ) {
    
  }

}

TEST_CASE( "infinite cylinder || x axis", "[CylinderX_Surface]" ) {

  SECTION ( "" ) {
    
  }

}

TEST_CASE( "infinite cylinder|| y axis", "[CylinderY_Surface]" ) {

  SECTION ( "" ) {
    
  }

}

TEST_CASE( "infinite cylinder || z axis", "[CylinderZ_Surface]" ) {

  SECTION ( "" ) {
    
  }

}

TEST_CASE( "infinite cone || x axis", "[ConeX_Surface]" ) {

  SECTION ( "" ) {
    
  }

}

TEST_CASE( "infinite cone || y axis", "[ConeY_Surface]" ) {

  SECTION ( "" ) {
    
  }

}

TEST_CASE( "infinite cone || z axis", "[ConeZ_Surface]" ) {

  SECTION ( "" ) {
    
  }

}
*/
TEST_CASE( "a cell", "[Region_t]" ) {

  SECTION ( "importance" ) {
    std::string aregionname = "a region name";
    double aimportance = 1.4;
    Region_t aRegion( aregionname, aimportance );
    REQUIRE( aRegion.importance() == 1.4 );
  }

  SECTION ( "SigmaT" ) {
    
  }

  SECTION ( "SigmaS" ) {
    
  }

  SECTION ( "SigmaC" ) {
    
  }

  SECTION ( "SigmaF" ) {
    
  }

  SECTION ( "setMaterial" ) {
    
  }

  SECTION ( "addSurface" ) {
    
  }

  SECTION ( "testPoint" ) {
    
  }

  SECTION ( "moveParticle" ) {
    
  }

  SECTION ( "surface_intersect" ) {
    
  }

  SECTION ( "collision_distance" ) {
    
  }

  SECTION ( "collision" ) {
    
  }

  SECTION ( "simulate_scatter" ) {
    
  }

}
