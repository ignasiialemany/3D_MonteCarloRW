//
// Created by Alemany Juvanteny, Ignasi on 02/03/2023.
//
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "../Geometry/polygon.h"
#include <Eigen/Dense>

TEST_CASE("Compute volume and surface of a polygon", "[polygon]")
{
    // Define vertices and faces of a tetrahedron
    Eigen::MatrixXd vertices(4, 3);
    Eigen::MatrixXd faces(4, 3);
    vertices << 0, 0, 0,
        1, 0, 0,
        0, 1, 0,
        0, 0, 1;
    faces << 1, 3, 2,
        1, 2, 4,
        1, 4, 3,
        2, 3, 4;

    polygon poly(vertices, faces);

    SECTION("Compute volume")
    {
        double expected_volume = 1.0 / 6.0;
        double computed_volume = poly.computeVolume();
        REQUIRE(std::abs(computed_volume - expected_volume) < 1e-6);
    }

    SECTION("Compute surface")
    {
        double expected_surface = 1.5 + std::sqrt(3) / 2;
        double computed_surface = poly.computeSurface();
        REQUIRE(std::abs(computed_surface - expected_surface) < 1e-6);
    }
}

TEST_CASE("Compute volume and surface area of a regular tetrahedron", "[polygon]")
{
    // Define vertices and faces of a regular tetrahedron
    Eigen::MatrixXd vertices(4, 3);
    Eigen::MatrixXd faces(4, 3);
    vertices << 0, 0, 0,
        1, 0, 0,
        0.5, std::sqrt(3) / 2, 0,
        0.5, std::sqrt(3) / 6, std::sqrt(6) / 3;
    faces << 1, 3, 2,
        1, 2, 4,
        2, 3, 4,
        3, 1, 4;

    polygon poly(vertices, faces);

    SECTION("Compute volume")
    {
        double expected_volume = std::sqrt(2.0) / 12.0;
        double computed_volume = poly.computeVolume();
        REQUIRE(std::abs(computed_volume - expected_volume) < 1e-6);
    }

    SECTION("Compute surface area")
    {
        double expected_surface = std::sqrt(3.0);
        double computed_surface = poly.computeSurface();
        REQUIRE(std::abs(computed_surface - expected_surface) < 1e-6);
    }
}

TEST_CASE("Test intersection of point and step with polygon face that is too close to an edge", "[polygon]")
{
    Eigen::MatrixXd vertices(4, 3);
    Eigen::MatrixXd faces(4, 3);
    vertices << 0, 0, 0,
        1, 0, 0,
        0, 1, 0,
        0, 0, 1;
    faces << 1, 3, 2,
        1, 2, 4,
        1, 4, 3,
        2, 3, 4;

    polygon poly(vertices, faces);

    // Define a point and step that intersect with the forth face but is too close to an edge
    Eigen::Vector3d point(0.5, 0.5, 0);
    Eigen::Vector3d step(0, 0, 1);

    // Check if the point and step intersect with any face
    std::pair<int, double> intersection_info = poly.intersection(point, step);

    // Check that there is no intersection
    REQUIRE(intersection_info.second == -1);
}

TEST_CASE("Test intersection of point and step with polygon face", "[polygon]")
{
    Eigen::MatrixXd vertices(4, 3);
    Eigen::MatrixXd faces(4, 3);
    vertices << 0, 0, 0,
        1, 0, 0,
        0, 1, 0,
        0, 0, 1;
    faces << 1, 3, 2,
        1, 2, 4,
        1, 4, 3,
        2, 3, 4;

    polygon poly(vertices, faces);

    // Define a point and step that intersect with the fourth face
    Eigen::Vector3d point(0.2, 0.2, 0);
    Eigen::Vector3d step(0, 0, 1);

    // Check if the point and step intersect with the fourth face
    std::pair<int, double> intersection_info = poly.intersection(point, step);

    // The expected intersection info is that the point and step intersect with the fourth face
    std::pair<int, double> expected_intersection_info(4, 0.6);

    // Check that the computed intersection info matches the expected intersection info
    REQUIRE(intersection_info.first == expected_intersection_info.first);
    REQUIRE(std::abs(intersection_info.second - expected_intersection_info.second) < 1e-6);
}
