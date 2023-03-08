//
// Created by Alemany Juvanteny, Ignasi on 26/02/2023.
//

#include "substrate.h"

substrate::substrate(std::vector<Eigen::MatrixXd> myo_vertices, std::vector<Eigen::MatrixXd> myo_faces)
{
    //Loop through myo_vertices and myo_faces and create a polygon for each
    for (int i = 0; i < myo_vertices.size(); ++i) {
        polygon myo(myo_vertices[i],myo_faces[i]);
        _myocytes.push_back(myo);
    }
}

void substrate::setTransform(transform &t)
{
}
