//
// Created by Alemany Juvanteny, Ignasi on 26/02/2023.
//

#ifndef INC_3DRANDOMWALK_TRANSFORM_H
#define INC_3DRANDOMWALK_TRANSFORM_H

#include <Eigen/Dense>

struct transform_info{
    Eigen::Vector3d local_position;
    //TODO: Might delete angle 
    double angle;
    int iX,iY,iZ;
};

class transform {
public:
    transform()=default;
    transform(double rot_in_y, bool shift_block): deg_rot_per_m_in_Y(rot_in_y), shift_block(shift_block) {};

    transform_info global2local(Eigen::Vector3d &global_position);
    transform_info local2global(Eigen::Vector3d &local_position, int iX, int iY, int iZ);

private:

    bool shift_block = true;
    double deg_rot_per_m_in_Y = 0.;
    
};


#endif //INC_3DRANDOMWALK_TRANSFORM_H
