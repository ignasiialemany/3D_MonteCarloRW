//
// Created by Alemany Juvanteny, Ignasi on 26/02/2023.
//

#ifndef INC_3DRANDOMWALK_TRANSFORM_H
#define INC_3DRANDOMWALK_TRANSFORM_H

#include <Eigen/Dense>
#include "utility_substrate.h"
#include <cmath>
#include <CGAL/Bbox_3.h>
#include "../Geometry/polygon.h"

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

    transform_info global2local(const Eigen::Vector3d &global_position) const;
    Eigen::Vector3d local2global(const Eigen::Vector3d &local_position, int iX, int iY, int iZ) const;
    void set_block(double dx, double dy, double dz, double minY = -2000. , double maxY = 5000.);
    CGAL::Bbox_3 get_solid_block()const {return _solid_block;};
    Polyhedron get_block()const {return _block;};
    std::shared_ptr<Tree_AABB> get_AABBtree()const {return _AABBtree;};
    Eigen::Vector3d get_block_size()const {return Eigen::Vector3d(dx,dy,dz);};

private:

    void create_block(Kernel::Point_3 min_point, Kernel::Point_3 max_point);
    double deg_rot_per_m_in_Y = 0.01;
    bool shift_block = true;
    double dx,dy,dz;
    bool isIdentity;
    std::vector<Kernel::Triangle_3> triangles;
    Eigen::MatrixXd y_slice_minmax;
    CGAL::Bbox_3 _solid_block;
    Polyhedron _block;
    std::shared_ptr<Tree_AABB> _AABBtree;
};


#endif //INC_3DRANDOMWALK_TRANSFORM_H
