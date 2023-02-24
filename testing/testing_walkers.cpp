#include "../MonteCarlo/walkers.h"
#include <iostream>
#include <catch2/catch_test_macros.hpp>

bool unique_seeds(unsigned int seed1, unsigned int seed2){
    walkers first_set_of_particles(10000, seed1,"constant");
    walkers second_set_of_particles(10000,seed2,"constant");

    return first_set_of_particles.get_seeds()==second_set_of_particles.get_seeds();
}

TEST_CASE( "Unique seeds is checked", "[unique seeds]" ) {
    REQUIRE( unique_seeds(1,2) == false );
    REQUIRE( unique_seeds(4,7456) == false );
    REQUIRE( unique_seeds(3,3) == true );
    REQUIRE( unique_seeds(1234,1234) == true );
}
