//
// Created by Alemany Juvanteny, Ignasi on 27/02/2023.
//

#include "utility.h"

intersection_ray_info
utility::rayIntersectsPolygon(const Eigen::Vector3d &point, const Eigen::Vector3d &direction, const Eigen::MatrixXd &V1,
                              const Eigen::MatrixXd &V2, const Eigen::MatrixXd &V3) {

    //Intersection is computed using V1, V2 and V3
    intersection_ray_info output;
    int Nverts = V1.rows();
    Eigen::MatrixXd orig = Eigen::MatrixXd::Ones(Nverts, 3).array().rowwise() * point.transpose().array();
    Eigen::MatrixXd dir = Eigen::MatrixXd::Ones(Nverts, 3).array().rowwise() * direction.transpose().array();

    // tolerances
    double eps = 1e-20;
    double zero = 0.0;

    output.intersect.resize(Nverts, 1);
    output.intersect.fill(false);
    // output.intersect = Eigen::VectorXd::Zero(Nverts);
    output.t = Eigen::VectorXd::Ones(Nverts) * -1; // infinite is replaced by -1 in matlab
    output.u = Eigen::VectorXd::Ones(Nverts) * -1;
    output.v = Eigen::VectorXd::Ones(Nverts) * -1;

    // some pre-calculations
    Eigen::MatrixXd edge1 = V2 - V1; // find vectors for two edges sharing V1
    Eigen::MatrixXd edge2 = V3 - V1;
    Eigen::MatrixXd tvec = orig - V1; // vector from V1 to ray origin

    Eigen::MatrixXd pvec = crossMat(dir, edge2);
    Eigen::MatrixXd qvec = crossMat(tvec, edge1);
    Eigen::VectorXd det = (edge1.array() *
                           pvec.array()).rowwise().sum(); // determinant of the matrix M = dot(edge1, pvec)

    // find faces parallel to the ray
    // Eigen::Array<bool, Eigen::Dynamic, 1> angleOK = (det.array().abs2() > eps).array(); // if det ~ 0 then ray lies in the triangle plane
    Eigen::Array<bool, Eigen::Dynamic, 1> angleOK = (det.array().abs() >
                                                     eps).array(); // if det ~ 0 then ray lies in the triangle plane

    if ((angleOK == false).all()) {
        printf("polygon::TriangleRayIntersection::angle not within tolerance\n");
        return output;
    }
    // To do : change to avoid division by zero

    // calculate all variables for all line/triangle pairs
    output.u = (tvec.array() * pvec.array()).rowwise().sum().array() / det.array();
    output.v = (dir.array() * qvec.array()).rowwise().sum().array() / det.array();
    output.t = (edge2.array() * qvec.array()).rowwise().sum().array() / det.array();

    // test if line/plane intersection is within the triangle
    Eigen::Array<bool, Eigen::Dynamic, 1> ok = (angleOK and (output.u.array() >= -zero).array() and
                                                (output.v.array() >= -zero).array() and
                                                ((output.u + output.v).array() <= 1.0 + zero).array());
    // std::cout << ok.cast<double>().sum() << std::endl;
    output.intersect = (ok and (output.t.array() >= -zero).array() and (output.t.array() <= 1.0 + zero).array());

    return output;
}

