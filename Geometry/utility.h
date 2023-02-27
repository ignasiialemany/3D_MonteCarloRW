//
// Created by Alemany Juvanteny, Ignasi on 27/02/2023.
//

#ifndef INC_3DRANDOMWALK_UTILITY_H
#define INC_3DRANDOMWALK_UTILITY_H

#include<Eigen/Dense>
#include<cmath>

struct intersection_ray_info {
    Eigen::VectorXd t, u, v;
    Eigen::Array<bool, Eigen::Dynamic, 1> intersect;
};

class utility {

public:
    Eigen::Array<bool, Eigen::Dynamic, 1>
    rayIntersectsPolygon(const Eigen::Vector3d &point, const Eigen::Vector3d &direction, const Eigen::MatrixXd &V1,
                         const Eigen::MatrixXd &V2, const Eigen::MatrixXd &V3);

};


#endif //INC_3DRANDOMWALK_UTILITY_H
