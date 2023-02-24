#include "../MonteCarlo/simulation.h"
#include <iostream>
#include <stdexcept>
#include <catch2/catch_test_macros.hpp>

bool check_seeding(int seed1,double min_x,double max_x, double min_y, double max_y, double min_z, double max_z){

    walkers particles(100000, seed1,"constant");
    simulation rw_particles(particles);

    Eigen::VectorXd box(6);
    box << min_x,max_x,min_y,max_y,min_z,max_z;

    auto isValid = rw_particles.seedParticlesInBox(box);
    if (!isValid){
        throw std::runtime_error("error");
    }

    auto positions = rw_particles.getPositions();

    for (int i=0; i<positions.cols(); i++){
        auto column = positions.col(i);
        double column_min = column.minCoeff();
        double column_max = column.maxCoeff();
        if (column_min < box(2*i) or column_max > box(2*i+1)){
            return false;
        }
        else{
            continue;
        }
    }
    return true;
}

bool unique_seeding(int seed1, int seed2,  double min_x,double max_x, double min_y, double max_y, double min_z, double max_z){
    walkers particles_1(100000, seed1,"constant");
    walkers particles_2(100000, seed2,"constant");
    simulation rw_particles1(particles_1);
    simulation rw_particles2(particles_2);

    Eigen::VectorXd box(6);
    box << min_x,max_x,min_y,max_y,min_z,max_z;

    rw_particles1.seedParticlesInBox(box);
    rw_particles2.seedParticlesInBox(box);

    return rw_particles1.getPositions()==rw_particles2.getPositions();
}

TEST_CASE( "Check particles are bounded within box", "[check_seeding]" ) {
    REQUIRE(check_seeding(13412,0.,100.,0.,100.,0.,700.) == true );
    REQUIRE(check_seeding(13412,235.,2345.,212.,1512.,1.,5.) == true );
    SECTION("error condition"){
    REQUIRE_THROWS_AS(check_seeding(1,235.,2345.,2000.,1512.,1.,5.),std::runtime_error);
    REQUIRE_THROWS_AS(check_seeding(1800,235.,2345.,2000.,1512.,1.,5.),std::runtime_error);
    }
}

TEST_CASE( "Check unique seeding","[unique_seeding]"){
    REQUIRE(unique_seeding(1230,1234,0,100,0,100,0,600)==false);
    REQUIRE(unique_seeding(1230,1230,0,100,0,100,0,600)==true);
}