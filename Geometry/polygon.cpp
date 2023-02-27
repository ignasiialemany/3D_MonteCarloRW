//
// Created by Alemany Juvanteny, Ignasi on 27/02/2023.
//

#include "polygon.h"

polygon::polygon(Eigen::MatrixXd vertices_input, Eigen::MatrixXd faces_input) {
    // Assigning vertices and faces
    vertices = vertices_input;
    faces = faces_input;

    // Number of faces and vertices
    n_vertices = vertices_input.rows();
    n_faces = faces_input.rows();

    _volume = polygon::computeVolume(vertices_input, faces_input);
    _surface = polygon::computeSurface(vertices_input, faces_input);

    // Computing the bounding box range
    bb_range.resize(6);
    minXYZ << vertices_input(Eigen::all, 0).minCoeff(), vertices_input(Eigen::all, 1).minCoeff(), vertices_input(
            Eigen::all, 2).minCoeff();
    maxXYZ << vertices_input(Eigen::all, 0).maxCoeff(), vertices_input(Eigen::all, 1).maxCoeff(), vertices_input(
            Eigen::all, 2).maxCoeff();
    bb_range << minXYZ(0), maxXYZ(0), minXYZ(1), maxXYZ(1), minXYZ(2), maxXYZ(2);

    boundingbox.initialize(bb_range);
}

double polygon::computeVolume(Eigen::MatrixXd &vertices_input, Eigen::MatrixXd &faces_input) {
    double volume = 0.0;
    Eigen::Vector3d centroid = vertices_input.colwise().mean();
    for (int i = 0; i < faces_input.rows(); i++) {
        // Extract the vertices of the face
        Eigen::Vector3d v1 = vertices_input.row(faces_input(i, 0)).transpose();
        Eigen::Vector3d v2 = vertices_input.row(faces_input(i, 1)).transpose();
        Eigen::Vector3d v3 = vertices_input.row(faces_input(i, 2)).transpose();

        // Compute the volume of the pyramid formed by the face and the centroid
        double pyramid_volume = fabs((v2 - v1).cross(v3 - v1).dot(centroid - v1)) / 6.0;

        // Add the volume of the pyramid to the total volume
        volume += pyramid_volume;
    }
    return volume;
}

double polygon::computeSurface(Eigen::MatrixXd &vertices_input, Eigen::MatrixXd &faces_input) {
    double surface = 0.0;
    for (int i = 0; i < faces_input.rows(); i++) {
        // Extract the vertices of the face
        Eigen::Vector3d v1 = vertices_input.row(faces_input(i, 0)).transpose();
        Eigen::Vector3d v2 = vertices_input.row(faces_input(i, 1)).transpose();
        Eigen::Vector3d v3 = vertices_input.row(faces_input(i, 2)).transpose();

        double face_area = 0.5 * (v2 - v1).cross(v3 - v1).norm();

        // Add the volume of the pyramid to the total volume
        surface += face_area;
    }
    return surface;
}

Eigen::MatrixXd polygon::getOrderedVertices(int column) {
    return vertices(faces(Eigen::all, column).array() - 1, Eigen::all);
}

void polygon::intersection(Eigen::Vector3d &orig, Eigen::Vector3d &dir) {

}
