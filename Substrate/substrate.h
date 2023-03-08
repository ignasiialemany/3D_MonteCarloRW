//
// Created by Alemany Juvanteny, Ignasi on 26/02/2023.
//

#ifndef INC_3DRANDOMWALK_SUBSTRATE_H
#define INC_3DRANDOMWALK_SUBSTRATE_H

#include "transform.h"
#include "../Geometry/polygon.h"

class substrate {

public:
    substrate()=default;
    substrate(std::vector<Eigen::MatrixXd> myo_vertices, std::vector<Eigen::MatrixXd> myo_faces);
    void setTransform(transform &t);
    transform_info getGlobalToLocal(Eigen::Vector3d &global_position){return _transform.global2local(global_position);};
    Eigen::Vector3d getLocalToGlobal(Eigen::Vector3d &local_position, int iX, int iY, int iZ){return _transform.local2global(local_position,iX,iY,iZ);};
    
private:
    std::vector<polygon> _myocytes;
    transform _transform;
};

#endif //INC_3DRANDOMWALK_SUBSTRATE_H
