//
// Created by Alemany Juvanteny, Ignasi on 27/02/2023.
//

#ifndef INC_3DRANDOMWALK_POLYGON_H
#define INC_3DRANDOMWALK_POLYGON_H

#include <Eigen/Dense>
#include "boundingbox.h"
#include "utility.h"
#include <iostream>

class polygon {

public:
    polygon() = default;

    polygon(const Eigen::MatrixXd &vertices_input, const Eigen::MatrixXd &faces_input);
    //polygon(std::string box_type, Eigen::MatrixXd bounding_box);

    double computeVolume();

    double computeSurface();

    std::pair<int,double> intersection(const Eigen::Vector3d &orig, const Eigen::Vector3d &dir);


    double _volume=0;
    double _surface=0;
    int n_vertices;
    int n_faces;
    Eigen::MatrixXd vertices;
    Eigen::MatrixXd faces;
    boundingbox bounding_box;
    Eigen::MatrixXd V1;
    Eigen::MatrixXd V2;
    Eigen::MatrixXd V3;

};


#endif //INC_3DRANDOMWALK_POLYGON_H
