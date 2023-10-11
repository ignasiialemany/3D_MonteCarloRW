//
// Created by Alemany Juvanteny, Ignasi on 26/02/2023.
//

#ifndef INC_3DRANDOMWALK_TRANSFORM_H
#define INC_3DRANDOMWALK_TRANSFORM_H

#include <Eigen/Dense>
#include "utility_substrate.h"
#include <cmath>
#include <CGAL/Bbox_3.h>
#include "../Geometry/polyhedronSet.h"

struct transform_info{
    Eigen::Vector3d local_position;
    //TODO: Might delete angle 
    double angle;
    int iX,iY,iZ;
};


class transform {

public:
    transform()=default;
    transform(double rot_in_y, bool shift_block, bool isIdentity): deg_rot_per_m_in_Y(rot_in_y), shift_block(shift_block), isIdentity(isIdentity){};

    void precomputeStrainedBlock();
    transform_info global2local(const Eigen::Vector3d &global_position) const;
    Eigen::Vector3d local2global(const Eigen::Vector3d &local_position, int iX, int iY, int iZ) const;
    void set_block(double dx, double dy, double dz, double minY = -2000. , double maxY = 5000.);
    Eigen::Vector3d get_block_size()const {return Eigen::Vector3d(dx,dy,dz);};
    void precomputeTransform(Eigen::VectorXd &sequence_dt, std::function<double(double)> &strain);
    bool isTransformIdentity() const {return isIdentity;};
    Eigen::Vector3d get_block_centroid() const {return _centroid_block;};
    const polygon& getBlockCurrentTime(int sequence_index) const {return _block.getPolygon(sequence_index);};

private:

    void create_block(Kernel::Point_3 min_point, Kernel::Point_3 max_point);
    double deg_rot_per_m_in_Y = 0.01;
    void calculate_minmax_slices(Eigen::MatrixXd& slice_minmax, double minrange, double maxrange);
    Eigen::MatrixXd y_slice_minmax;
    Eigen::MatrixXd x_slice_minmax;
    Eigen::MatrixXd z_slice_minmax;

    bool shift_block = true;
    double dx,dy,dz;
    bool isIdentity;
    std::vector<Kernel::Triangle_3> triangles;
    polyhedronSet _block;
    Eigen::Vector3d _centroid_block;
};


#endif //INC_3DRANDOMWALK_TRANSFORM_H
