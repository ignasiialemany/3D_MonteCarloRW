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

private:
    std::vector<polygon> _myocytes;
    
};

#endif //INC_3DRANDOMWALK_SUBSTRATE_H
