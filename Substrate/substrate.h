//
// Created by Alemany Juvanteny, Ignasi on 26/02/2023.
//

#ifndef INC_3DRANDOMWALK_SUBSTRATE_H
#define INC_3DRANDOMWALK_SUBSTRATE_H

#include "transform.h"
#include "../Geometry/polyhedronSet.h"
#include <list>
#include <vector>
#include <map>
#include <list>
#include <CGAL/point_generators_3.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <memory>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Bounded_normal_change_filter.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_cost.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/angle_and_area_smoothing.h>
#include <CGAL/Polygon_mesh_processing/detect_features.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
// #include <CGAL/Polygon_mesh_processing/smooth_mesh.h>
#include <CGAL/Polygon_mesh_processing/angle_and_area_smoothing.h>
#include <CGAL/Polygon_mesh_processing/smooth_shape.h>
#include <CGAL/Named_function_parameters.h>
#include <functional>
#include <chrono>
#include <limits>
#include <boost/variant.hpp>
#include <boost/filesystem.hpp>

typedef CGAL::Search_traits_3<Kernel> TreeTraits;
typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
typedef Neighbor_search::Tree Tree;

class substrate
{

public:
    // add deconstructor
    substrate() = default;
    substrate(std::vector<Eigen::MatrixXd> &myo_vertices, std::vector<Eigen::MatrixXd> &myo_faces, transform &t, std::function<double(double)> strain_f = [](double){return 0.;});
    substrate(const std::vector<std::string>& polygons, transform &t, std::function<double(double)> strain_f = [](double){return 0.;});

    transform_info getLocalFromGlobal(const Eigen::Vector3d &global_position) const { return _transform.global2local(global_position); };
    Eigen::Vector3d getGlobalFromLocal(const Eigen::Vector3d &local_position, int iX, int iY, int iZ) const { return _transform.local2global(local_position, iX, iY, iZ); };
    bool containsPoint(int index_polygon, const Eigen::Vector3d &point, int index_sequence);
    boost::optional<std::tuple<int, double, Eigen::Vector3d>> intersectPolygon(const Eigen::Vector3d &point, const Eigen::Vector3d &step, int index_sequence) const;
    boost::optional<std::tuple<int, double, Eigen::Vector3d>> intersectionBlock(const Eigen::Vector3d &point, const Eigen::Vector3d &step, int index_sequence) const;
    int searchPolygon(const Eigen::Vector3d &point,  int index_sequence, const std::string &frameOfReference = "local") const;
    Eigen::Vector3d get_block_size() const { return _transform.get_block_size(); };
    void setVoxel(const Eigen::VectorXd &voxel) { _voxel = voxel; };
    Eigen::VectorXd getVoxel() const { return _voxel; };
    void setBoundaryType(const std::string &type){_boundary_type = type;};
    double getStrain(double time) const {return _strain(time);};
    bool isTransformIdentity() const { return _transform.isTransformIdentity(); };
    std::string getBoundaryType() const { return _boundary_type; };
    Eigen::Vector3d get_block_centroid() const { return _transform.get_block_centroid(); };
    void preComputeSubstrate(Eigen::VectorXd sequence_dt);
    mutable std::vector<polyhedronSet> _myocytes;

private:
    std::string _boundary_type;
    std::vector<Kernel::Point_3> _points;
    transform _transform;
    std::unique_ptr<Tree> _tree;
    std::map<Kernel::Point_3, int> _map_centroid_to_polygon;
    Eigen::VectorXd _voxel;
    std::function<double(double)> _strain;

    void write_ply_polyhedron(polygon &poly, const std::string &filename);
    void write_ply_surface_mesh(const Mesh &mesh, const std::string &filename);
};

#endif // INC_3DRANDOMWALK_SUBSTRATE_H
