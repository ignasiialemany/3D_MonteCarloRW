//
// Created by Alemany Juvanteny, Ignasi on 27/02/2023.
//

#include "utility.h"

std::vector<std::pair<int,double>> utility::rayIntersectsPolygon(const Eigen::Vector3d &point, const Eigen::Vector3d &direction, const Eigen::MatrixXd &V1, const Eigen::MatrixXd &V2,const Eigen::MatrixXd &V3) {
    
    const double EPSILON = 1e-10;
    const int num_faces = V1.rows();
    const double step_length = direction.norm();

    //Store the index and distance to face
    std::vector<std::pair<int, double>> distance_to_face;

    for (int i = 0; i < num_faces; i++) {

        Eigen::Vector3d edge1 = V2.row(i) - V1.row(i);
        Eigen::Vector3d edge2 = V3.row(i) - V1.row(i);
        Eigen::Vector3d tvec = point - V1.row(i).transpose();
        Eigen::Vector3d pvec = direction.cross(edge2);
        double det = edge1.dot(pvec);

        //Checks that ray is not in face
        if (abs(det) <= EPSILON){
            continue;
        }

        Eigen::Vector3d qvec = tvec.cross(edge1);
        double u = tvec.dot(pvec)/det;
        double v = direction.dot(qvec)/det;
        double w = 1 - u - v; 
        double t = edge2.dot(qvec)/det;
        double t2 = edge1.dot(qvec)/det;

        //If face is too close to an edge or out of the face
        if (u < EPSILON || u > 1 - EPSILON || v < EPSILON || v > 1 - EPSILON || w < EPSILON || w > 1 - EPSILON){
            continue;
        }
        
        //If the distance is too close or the remaining step is too short
        if (t < EPSILON || t > 1 - EPSILON){
            continue;
        }
        
        distance_to_face.emplace_back(i, t*step_length);
    }
    return distance_to_face;
}


std::pair<int,double> utility::intersection(const Eigen::Vector3d &point, const Eigen::Vector3d &direction, const Eigen::MatrixXd &V1, const Eigen::MatrixXd &V2, const Eigen::MatrixXd &V3) {
    //returns a vector of pairs containing the face index and distance to membrane
    std::vector<std::pair<int,double>> info = utility::rayIntersectsPolygon(point,direction,V1,V2,V3);
    if (info.size()==0){
        //TODO: see if we can handle this with error after
        return std::make_pair(-1,-1);
    }
    double min_distance_to_face = info[0].second;
    int face_index = info[0].first;
    if (info.size() > 1){
        //Do this check and for now just output the error
        std::cout << "Multiple faces are intersecting?" << std::endl;

        //For now compute the smallest distance_to_face but let's see if this error shows up a lot
        for (auto & it : info) {
            if (it.second < min_distance_to_face) {
                min_distance_to_face = it.second;
                face_index = it.first;
            }
        }
    }

    return std::make_pair(face_index+1,min_distance_to_face);
}



Eigen::MatrixXd utility::getOrderedVertices(const Eigen::MatrixXd &vertices, const Eigen::MatrixXd &faces, int index) {
    Eigen::MatrixXd output = vertices(faces(Eigen::all, index).array() - 1, Eigen::all);
    return output;
}

Eigen::MatrixXd utility::crossMat(Eigen::MatrixXd &a, Eigen::MatrixXd &b) {
    Eigen::MatrixXd output = Eigen::MatrixXd::Zero(a.rows(), 3);
    for (int i = 0; i < a.rows(); i++){
        output(i, {0, 1, 2}) = a(i, {0, 1, 2}).cross(b(i, {0, 1, 2}));
    }
    return output;
}
