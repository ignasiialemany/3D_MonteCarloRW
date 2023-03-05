//
// Created by Alemany Juvanteny, Ignasi on 27/02/2023.
//

#ifndef INC_3DRANDOMWALK_UTILITY_H
#define INC_3DRANDOMWALK_UTILITY_H

#include<Eigen/Dense>
#include<cmath>
#include<iostream>

class utility {

public:

    static std::vector<std::pair<int,double>> rayIntersectsPolygon(const Eigen::Vector3d &point, const Eigen::Vector3d &direction, const Eigen::MatrixXd &V1, const Eigen::MatrixXd &V2,const Eigen::MatrixXd &V3);

    static std::pair<int,double> intersection(const Eigen::Vector3d &point, const Eigen::Vector3d &direction, const Eigen::MatrixXd &V1, const Eigen::MatrixXd &V2, const Eigen::MatrixXd &V3);

    static Eigen::MatrixXd getOrderedVertices(const Eigen::MatrixXd &vertices, const Eigen::MatrixXd &column, int index);

    static Eigen::MatrixXd crossMat(Eigen::MatrixXd &a, Eigen::MatrixXd &b);

};


#endif //INC_3DRANDOMWALK_UTILITY_H
