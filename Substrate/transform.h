//
// Created by Alemany Juvanteny, Ignasi on 26/02/2023.
//

#ifndef INC_3DRANDOMWALK_TRANSFORM_H
#define INC_3DRANDOMWALK_TRANSFORM_H

#include <Eigen/Dense>
#include "utility_substrate.h"
#include <cmath>

struct transform_info{
    Eigen::Vector3d local_position;
    //TODO: Might delete angle 
    double angle;
    int iX,iY,iZ;
};

class transform {
public:
    transform()=default;
    transform(double rot_in_y, bool shift_block, bool isIdentity): deg_rot_per_m_in_Y(rot_in_y), shift_block(shift_block), isIdentity(isIdentity){};

    transform_info global2local(Eigen::Vector3d &global_position);
    Eigen::Vector3d local2global(Eigen::Vector3d &local_position, int iX, int iY, int iZ);

private:

    double find_yslice(double y_global);
    Eigen::Vector3d rotate_y(Eigen::Vector3d &position, double &theta);
    bool shift_block = true;
    double deg_rot_per_m_in_Y = 0.;
    double dx,dy,dz;
    bool isIdentity;
    Eigen::MatrixXd y_slice_minmax;
};


#endif //INC_3DRANDOMWALK_TRANSFORM_H
