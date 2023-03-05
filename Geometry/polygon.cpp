//
// Created by Alemany Juvanteny, Ignasi on 27/02/2023.
//

#include "polygon.h"

polygon::polygon(const Eigen::MatrixXd &vertices_input, const Eigen::MatrixXd &faces_input) : vertices(vertices_input),
                                                                                faces(faces_input),
                                                                                n_vertices(vertices_input.rows()),
                                                                                n_faces(faces_input.rows()),
                                                                                V1(utility::getOrderedVertices(
                                                                                        vertices_input, faces_input,
                                                                                        0)),
                                                                                V2(utility::getOrderedVertices(
                                                                                        vertices_input, faces_input,
                                                                                        1)),
                                                                                V3(utility::getOrderedVertices(
                                                                                        vertices_input, faces_input,
                                                                                        2)) {

    // Computing the bounding box range
    Eigen::VectorXd bb_range(6);

    Eigen::VectorXd minXYZ(3);
    minXYZ << vertices_input.col(0).minCoeff(), vertices_input.col(1).minCoeff(), vertices_input.col(2).minCoeff();

    Eigen::VectorXd maxXYZ(3);
    maxXYZ << vertices_input.col(0).maxCoeff(), vertices_input.col(1).maxCoeff(), vertices_input.col(2).maxCoeff();

    bb_range << minXYZ, maxXYZ;

    std::cout << bb_range << std::endl;
    //Init bounding box
    bounding_box.initialize(bb_range);
}

double polygon::computeVolume() {
    if (!_volume) {
        double volume = 0.0;
        const Eigen::Vector3d centroid = vertices.colwise().mean();
        for (int i = 0; i < faces.rows(); i++) {
            
            // Extract the vertices of the face
            Eigen::Vector3d v1 = vertices.row(faces(i, 0)-1).transpose();
            Eigen::Vector3d v2 = vertices.row(faces(i, 1)-1).transpose();
            Eigen::Vector3d v3 = vertices.row(faces(i, 2)-1).transpose();

            // Compute the volume of the pyramid formed by the face and the centroid
            double pyramid_volume = fabs((v2 - v1).cross(v3 - v1).dot(centroid - v1)) / 6.0;

            // Add the volume of the pyramid to the total volume
            volume += pyramid_volume;
        }
        _volume = volume;
        return volume;
    } else {
        return _volume;
    }

}

double polygon::computeSurface() {
    if (!_surface) {
        double surface = 0.0;
        for (int i = 0; i < faces.rows(); i++) {
            // Extract the vertices of the face
            Eigen::Vector3d v1 = vertices.row(faces(i, 0)-1).transpose();
            Eigen::Vector3d v2 = vertices.row(faces(i, 1)-1).transpose();
            Eigen::Vector3d v3 = vertices.row(faces(i, 2)-1).transpose();

            double face_area = 0.5 * (v2 - v1).cross(v3 - v1).norm();

            // Add the volume of the pyramid to the total surface
            surface += face_area;
        }
        _surface = surface;
        return surface;
    } else {
        return _surface;
    }
}

std::pair<int, double> polygon::intersection(const Eigen::Vector3d &orig, const Eigen::Vector3d &dir) {
    if (V1.isZero()) {
        std::cout << "Polygon not initialized properly" << std::endl;
        exit(EXIT_FAILURE);
    }
    std::pair<int, double> intersection_info = utility::intersection(orig, dir, V1, V2, V3);
    return intersection_info;
}
