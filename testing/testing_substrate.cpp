#include <matioCpp/matioCpp.h>
#include <catch2/catch_test_macros.hpp>
#include "../Substrate/substrate.h"
#include <boost/filesystem.hpp>
#include <iostream>
#include <string>

TEST_CASE("Set substrate", "[substrate]")
{
    std::string file = "/Users/ia4118/CLionProjects/3D_RandomWalk/testing/data/geometry_1.mat";
    std::pair<std::vector<Eigen::MatrixXd>, std::vector<Eigen::MatrixXd>> myo = utility_substrate::read_mat_file(file);
    substrate sub(myo.first, myo.second);
    Eigen::Vector3d pos1(250., 250., 250.);
    int myo_index1 = sub.searchPolygon(pos1);
    REQUIRE(myo_index1 == -1);
}
