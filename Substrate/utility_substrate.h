//
// Created by Alemany Juvanteny, Ignasi on 08/03/2023.
//

#ifndef INC_3D_RANDOMWALK_UTILITY_SUBSTRATE_H
#define INC_3D_RANDOMWALK_UTILITY_SUBSTRATE_H

#include <matioCpp/matioCpp.h>
#include <Eigen/Dense>

class utility_substrate
{

public:
    static std::pair<std::vector<Eigen::MatrixXd>, std::vector<Eigen::MatrixXd>> read_mat_file(std::string &file);

    static Eigen::MatrixXd vertex_conversion(matvar_t *field);

    static Eigen::MatrixXd face_conversion(matvar_t *field);

    static double mod(double a, double b);

    static double find_yslice(double y_global, const Eigen::MatrixXd &y_slice_minmax);

    static Eigen::Vector3d rotate_y(const Eigen::Vector3d &position, double theta);
};

#endif // INC_3D_RANDOMWALK_UTILITY_SUBSTRATE_H
