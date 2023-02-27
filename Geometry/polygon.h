//
// Created by Alemany Juvanteny, Ignasi on 27/02/2023.
//

#ifndef INC_3DRANDOMWALK_POLYGON_H
#define INC_3DRANDOMWALK_POLYGON_H

#include <Eigen/Dense>
#include "utility.h"

class polygon {

public:
    polygon() = default;

    polygon(Eigen::MatrixXd vertices_input, Eigen::MatrixXd faces_input);
    //polygon(std::string box_type, Eigen::MatrixXd bounding_box);

    Eigen::MatrixXd vertices;
    Eigen::MatrixXd faces;
    int n_vertices;
    int n_faces;

    static double computeVolume(Eigen::MatrixXd &vertices_input, Eigen::MatrixXd &faces_input);

    static double computeSurface(Eigen::MatrixXd &vertices_input, Eigen::MatrixXd &faces_input);

    void intersection(Eigen::Vector3d &orig, Eigen::Vector3d &dir);


private:

    Eigen::MatrixXd getOrderedVertices(int column);

    double _volume;
    double _surface;
};


#endif //INC_3DRANDOMWALK_POLYGON_H
