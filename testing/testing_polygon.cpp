//
// Created by Alemany Juvanteny, Ignasi on 02/03/2023.
//
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "../Geometry/polygon.h"
#include <matioCpp/matioCpp.h>


TEST_CASE("Compute surface of polygon", "[polygon]")
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

    SECTION("Check if polyhedron contains a point")
    {
        Eigen::Vector3d point_inside(0.25, 0.25, 0.25);
        Eigen::Vector3d point_outside(2, 2, 2);

        REQUIRE(poly.containsPoint(point_inside) == true);
        REQUIRE(poly.containsPoint(point_outside) == false);
    }

    SECTION("Check if polyhedron contains a point")
    {
        // Check point containment
        Eigen::Vector3d inside_point(0.2, 0.2, 0.2);
        Eigen::Vector3d outside_point(2, 2, 2);
        Eigen::Vector3d on_face_point(0.5, std::sqrt(3) / 6, std::sqrt(6) / 6);
        Eigen::Vector3d on_vertex_point(1, 0, 0);
        Eigen::Vector3d on_edge_point(0.25, 0, 0);
        Eigen::Vector3d outside_eps(0.8, 0.05, -1e-8);

        REQUIRE(poly.containsPoint(inside_point) == true);
        REQUIRE(poly.containsPoint(outside_point) == false);
        REQUIRE(poly.containsPoint(on_face_point) == true);
        REQUIRE(poly.containsPoint(on_vertex_point) == true);
        REQUIRE(poly.containsPoint(on_edge_point) == true);
        REQUIRE(poly.containsPoint(outside_eps) == false);
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
    
    boost::optional<std::pair<int,double>> intersection = poly.intersection(point, step);

    if (!intersection)
    {
        FAIL("The point and step do not intersect with the polygon");
    }
    else{
        REQUIRE((*intersection).first == 0);
        REQUIRE((*intersection).second == 0);
    }
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
    boost::optional<std::pair<int, double>> intersection_info = poly.intersection(point, step);

    // The expected intersection info is that the point and step intersect with the fourth face
    std::pair<int, double> expected_intersection_info(0, 0);

    if (!intersection_info)
    {
        FAIL("The point and step do not intersect with the polygon");
    }
    else{
        REQUIRE((*intersection_info).first == expected_intersection_info.first);
        REQUIRE((*intersection_info).second == expected_intersection_info.second);
    }
}

TEST_CASE("Test intersection of point and step with polygon face with distance", "[polygon]")
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
    Eigen::Vector3d point(0.2, 0.2, 0.1);
    Eigen::Vector3d step(0, 0, 1);

    // Check if the point and step intersect with the fourth face
    boost::optional<std::pair<int, double>> intersection_info = poly.intersection(point, step);

    // The expected intersection info is that the point and step intersect with the fourth face
    std::pair<int, double> expected_intersection_info(3, 0.5);

    if (!intersection_info)
    {
        FAIL("The point and step do not intersect with the polygon");
    }
    else{
        REQUIRE((*intersection_info).first == expected_intersection_info.first);
        REQUIRE((*intersection_info).second == expected_intersection_info.second);
    }
}