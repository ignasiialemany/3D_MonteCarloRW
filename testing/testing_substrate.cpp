#include <matioCpp/matioCpp.h>
#include <catch2/catch_test_macros.hpp>
#include "../Substrate/substrate.h"
#include <boost/filesystem.hpp>
#include <iostream>
#include <string>
#include <random>
#include <Eigen/Dense>
#include <string>
#include <iomanip>
#include <boost/variant.hpp>
TEST_CASE("Set substrate", "[substrate]")
{
    std::string file = "/Users/ia4118/CLionProjects/3D_RandomWalk/testing/data/geometry_1.mat";
    std::pair<std::vector<Eigen::MatrixXd>, std::vector<Eigen::MatrixXd>> myo = utility_substrate::read_mat_file(file);
    substrate sub(myo.first, myo.second);

    Eigen::MatrixXd mat(20, 3);
    mat << 0, 51, 95,
        225, 208, 27,
        23, 265, 85,
        458, 149, 65,
        408, 13, 6,
        260, 262, 0,
        188, 26, 52,
        337, 230, 117,
        415, 206, 11,
        321, 162, 88,
        446, 298, 33,
        23, 287, 41,
        310, 295, 124,
        179, 96, 123,
        354, 294, 82,
        35, 246, 111,
        133, 170, 96,
        234, 92, 34,
        176, 65, 61,
        440, 355, 7;

    // From matlab code
    Eigen::VectorXd expected_indices(20);
    expected_indices << -1, 491, 51, 1105, 977, -1, 393, 797, 995, 767, 1055, 57, 731, -1, -1, 59, -1, -1, 395, 1053;

    Eigen::VectorXd myo_indices(mat.rows());

    for (int i = 0; i < mat.rows(); i++)
    {
        Eigen::Vector3d pos(mat(i, 0), mat(i, 1), mat(i, 2));
        myo_indices(i) = sub.searchPolygon(pos);
    }

    REQUIRE(myo_indices == expected_indices);

    // Check intersections
    Eigen::MatrixXi mat_intersection(50, 3);
    std::default_random_engine generator(1201231);

    std::uniform_int_distribution<int> dist1(20, 400);
    std::uniform_int_distribution<int> dist2(20, 360);
    std::uniform_int_distribution<int> dist3(20, 110);

    for (int i = 0; i < mat_intersection.rows(); i++)
    {
        mat_intersection(i, 0) = (int)dist1(generator);
        mat_intersection(i, 1) = (int)dist2(generator);
        mat_intersection(i, 2) = (int)dist3(generator);
    }

    std::vector<int> intersection_indices;
    std::vector<std::string> intersection_distances;

    for (int i = 0; i < mat_intersection.rows(); i++)
    {
        Eigen::Vector3d pos(mat_intersection(i, 0), mat_intersection(i, 1), mat_intersection(i, 2));
        Eigen::Vector3d step(10.51, 0.12, 9.12);
        boost::variant<bool, std::pair<int, double>> result = sub.intersection(pos, step);
        if (result.which() == 0)
        {
            intersection_indices.push_back(-1);
            intersection_distances.push_back("-1");
        }
        else if (result.which() == 1)
        {
            intersection_indices.push_back(boost::get<std::pair<int, double>>(result).first);
            std::ostringstream stream;
            stream << std::fixed << std::setprecision(10) << boost::get<std::pair<int, double>>(result).second;
            double rounded_val = std::round(stod(stream.str()) * 100000) / 100000.0;
            stream.str("");
            stream << std::fixed << std::setprecision(5) << rounded_val;
            intersection_distances.push_back(stream.str());
        }
    }

    // Matlab values
    std::vector<int> expected_inter_ind{-1, 721, 446, -1, 308, -1, 718, 636, -1, -1, 562, -1, -1, -1, -1, 511, 612, -1, 242, 546, 352, -1, -1, -1, 244, 696, 500, 731, 331, 347, 869, 285, 191, 826, 390, 459, 689, -1, 308, 689, 396, 851, -1, 343, -1, -1, 734, 656, 511, 787};
    std::vector<std::string> expected_inter_distances =  {"-1", "11.99485", "1.64477", "-1", "8.03472", "-1", "3.53643", "0.29116",
  "-1", "-1", "1.43554", "-1", "-1", "-1", "-1", "5.61721", "7.39377", "-1",
  "10.64748", "4.02009", "10.67212", "-1", "-1", "-1", "9.84459", "13.43410",
  "7.77907", "2.50008", "12.20142", "1.75464", "11.45332", "9.59238", "1.30711", "0.81742", "11.31343", "7.66653", "8.44054", "-1", "2.10813", "7.82980", "7.55792", "2.02607", "-1", "6.83735", "-1", "-1", "8.16040", "4.82247", "1.78377", "12.19611"};

    REQUIRE(intersection_indices == expected_inter_ind);
    REQUIRE(intersection_distances == expected_inter_distances);
}
