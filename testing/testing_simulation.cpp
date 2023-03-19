#include "../MonteCarlo/simulation.h"
#include <iostream>
#include <stdexcept>
#include <catch2/catch_test_macros.hpp>
#include <boost/filesystem.hpp>

/*
TEST_CASE("Seeding particles in bounding box", "[simulation]")
{

    int seed = 1;
    walkers particles(1000, seed);
    simulation sim(particles);

    SECTION("Valid bounding box")
    {
        Eigen::VectorXd bounding_box(6);
        bounding_box << 0.0, 1.0, 0.0, 1.0, 0.0, 1.0;
        bool isSeedingSuccessful = sim.seedParticlesInBox(bounding_box);
        REQUIRE(isSeedingSuccessful == true);
    }

    SECTION("Invalid bounding box (min > max)")
    {
        Eigen::VectorXd bounding_box(6);
        bounding_box << 1.0, 0.0, 1.0, 0.0, 1.0, 0.0;
        REQUIRE(sim.seedParticlesInBox(bounding_box) == false);
    }

    SECTION("Check particles in bounding box")
    {
        Eigen::VectorXd bounding_box(6);
        bounding_box << 0.0, 250., 0.0, 250., 0.0, 800.;
        if (sim.seedParticlesInBox(bounding_box))
        {
            for (int i = 0; i < particles.get_number_of_particles(); i++)
            {
                Particle p = particles.get_particle(i);
                REQUIRE(p.position(0) >= bounding_box(0));
                REQUIRE(p.position(0) <= bounding_box(1));
                REQUIRE(p.position(1) >= bounding_box(2));
                REQUIRE(p.position(1) <= bounding_box(3));
                REQUIRE(p.position(2) >= bounding_box(4));
                REQUIRE(p.position(2) <= bounding_box(5));
            }
        }
    }
}
*/
TEST_CASE("one_dt", "[one_dt]")
{

    oneGenerator fixed_rng_engine;

    auto cwd = boost::filesystem::current_path();
    auto parent_path = cwd.parent_path();
    std::string parent_path_str = parent_path.string();
    std::string grandparent_path_str = parent_path_str.substr(0, parent_path_str.size() - 5);

    std::string sub_path = "testing/data/geometry_1.mat";
    std::string sub_full_path = grandparent_path_str + sub_path;

    simulation sim(1,1);
    std::string step_type = "constant";
    std::string transit_model = "whatever";
    sim.set_parameters(2, 3, step_type, transit_model, 2.5, 1., 0.05);
    auto myo = utility_substrate::read_mat_file(sub_full_path);
    substrate sub(myo.first, myo.second);

    transform t(0.01, true, false);
    t.set_block(495.3992, 392.3432, 126.5612);
    sub.setTransform(t);

    sim.get_particle(0).position << 1433.5427938, 910.0219397, 3206.792634463;
    int myo_index = sub.searchPolygon(sim.get_particle(0).position, "global");
    sim.get_particle(0).myocyte_index = myo_index;
    sim.one_dt<oneGenerator>(sim.get_particle(0), sub, fixed_rng_engine, 1);
    Eigen::Vector3d expected_position(1434.9570073623731, 911.4361532623731, 3208.2068479267462);
    Eigen::Vector3d diff = expected_position - sim.get_particle(0).position;
    REQUIRE(std::abs(diff.norm()) < 1e-6);


    //This case checks intersection with block, intersection with cardiomyoycte (reflection) and  no intersection
    sim.get_particle(0).position << 905.37515026, 1542.7427839, 586.83012998;
    myo_index = sub.searchPolygon(sim.get_particle(0).position, "global");
    sim.get_particle(0).myocyte_index = myo_index;
    sim.one_dt<oneGenerator>(sim.get_particle(0), sub, fixed_rng_engine, 0.4886875);
    Eigen::Vector3d expected_position_2(906.13982179, 1544.271804, 585.858874785);
    Eigen::Vector3d diff_2 = expected_position_2 - sim.get_particle(0).position;
    REQUIRE(std::abs(diff_2.norm()) < 1e-6);


     //This case checks intersection with block, intersection with cardiomyoycte (reflection) and  no intersection
    sim.get_particle(0).position << 979.858916463534,614.924794934872,7289.33203526935;
    myo_index = sub.searchPolygon(sim.get_particle(0).position, "global");
    sim.get_particle(0).myocyte_index = myo_index;
    sim.one_dt<oneGenerator>(sim.get_particle(0), sub, fixed_rng_engine, 0.4886875);
    Eigen::Vector3d expected_position_3(981.422066328390,616.487944799729,7290.89518513421);
    Eigen::Vector3d diff_3 = expected_position_3 - sim.get_particle(0).position;
    REQUIRE(std::abs(diff_3.norm()) < 1e-6);
    


}
