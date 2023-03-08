//
// Created by Alemany Juvanteny, Ignasi on 27/02/2023.
//

#ifndef INC_3DRANDOMWALK_UTILITY_H
#define INC_3DRANDOMWALK_UTILITY_H

#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>

class utility
{
    //using Point3 = boost::geometry::model::point<double, 3, bg::cs::cartesian>;
    //using Polygon3 = boost::geometry::model::polygon<Point3>;
    //using Box3 = boost::geometry::model::box<Point3>;

public:

    int x;
    //static std::vector<std::tuple<int, double, bool>> rayIntersectsPolygon(const Eigen::Vector3d &point, const Eigen::Vector3d &direction, const Eigen::MatrixXd &V1, const Eigen::MatrixXd &V2, const Eigen::MatrixXd &V3);

    //static std::pair<int, double> intersection(const Eigen::Vector3d &point, const Eigen::Vector3d &direction, const Eigen::MatrixXd &V1, const Eigen::MatrixXd &V2, const Eigen::MatrixXd &V3);

    //static Eigen::MatrixXd getOrderedVertices(const Eigen::MatrixXd &vertices, const Eigen::MatrixXd &column, int index);

    //static Eigen::MatrixXd crossMat(Eigen::MatrixXd &a, Eigen::MatrixXd &b);

    //static bool polygonContainsPoint(const Eigen::MatrixXd &vertices, const Eigen::MatrixXd &faces, const Eigen::Vector3d &point,const Eigen::MatrixXd &V1, const Eigen::MatrixXd &V2, const Eigen::MatrixXd &V3);

    //static Eigen::Vector3d computeBarycentricCoordinates(const Polygon3 &face, const Point3D &point);

};

#endif // INC_3DRANDOMWALK_UTILITY_H
