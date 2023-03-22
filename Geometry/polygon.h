//
// Created by Alemany Juvanteny, Ignasi on 27/02/2023.
//

#ifndef INC_3DRANDOMWALK_POLYGON_H
#define INC_3DRANDOMWALK_POLYGON_H

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/bounding_box.h>
#include <CGAL/Ray_3.h>
#include <CGAL/Triangle_3.h>
#include <CGAL/Vector_3.h>
#include <CGAL/intersection_3.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Point_3.h>
#include <Eigen/Dense>
#include <cassert>
#include <vector>
#include <boost/optional.hpp>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_profile.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_cost.h>
#include <CGAL/pca_estimate_normals.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/mst_orient_normals.h>
#include <CGAL/Poisson_reconstruction_function.h>
#include <CGAL/property_map.h>
#include <CGAL/compute_average_spacing.h>
#include <vector>
#include <fstream>
#include <Eigen/Dense>
#include <CGAL/Object.h>
#include <memory>
#include <CGAL/IO/PLY.h>
#include <fstream>
#include <CGAL/Polygon_mesh_processing/remesh.h>

//typedef CGAL::Simple_cartesian<long double> Kernel;
typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef Polyhedron::HalfedgeDS HalfedgeDS;
typedef CGAL::Surface_mesh<Kernel::Point_3> Mesh;

// Define the AABB traits
typedef std::vector<Kernel::Triangle_3>::iterator Iterator;
typedef CGAL::AABB_triangle_primitive<Kernel, Iterator> Primitive;
typedef CGAL::AABB_traits<Kernel, Primitive> AABB_triangle_traits;
typedef CGAL::AABB_tree<AABB_triangle_traits> Tree_AABB;

class polygon
{
public:
    polygon() = default;
    //Add deconstructor
    polygon(const Eigen::MatrixXd &vertices_input, const Eigen::MatrixXd &faces_input);
    polygon(Polyhedron &poly);

    double computeVolume();
    double computeSurface();
    Polyhedron getPolyhedron() { return _poly; };
    Polyhedron getBbox() { return _bbox; };
    CGAL::Bbox_3 getSolidBbox() { return _solid_bbox; };
    // boost::variant<bool, std::pair<int, double>> intersection(const Eigen::Vector3d &point, const Eigen::Vector3d &step);
    boost::optional<std::pair<int, double>> intersection(const Eigen::Vector3d &point, const Eigen::Vector3d &step) const;

    // bool containsPoint(const Eigen::Vector3d &point);
    bool containsPoint(const Eigen::Vector3d &point) const;
    Eigen::Vector3d getNormalVector(int index_face) const;
    std::vector<CGAL::Point_3<Kernel>> getVertices() const { return _vertices; };
    std::vector<std::vector<std::size_t>> getFaces() const { return _faces; };
    bool isPolygonClosed();
    Mesh getMesh() { return mesh; };

private:
    // TODO: Might delete some of these variables, probably _vertices and _faces
    Polyhedron _poly;
    Polyhedron _bbox;
    CGAL::Bbox_3 _solid_bbox;
    Mesh mesh;
    void createBbox(Kernel::Point_3 min_point, Kernel::Point_3 max_point);
    void createPolygon(const Eigen::MatrixXd &vertices, const Eigen::MatrixXd &faces);
    std::vector<CGAL::Point_3<Kernel>> _vertices;
    std::vector<std::vector<std::size_t>> _faces;
    std::unique_ptr<Tree_AABB> _AABBtree;
    std::vector<Kernel::Triangle_3> triangle_faces;

};

#endif // INC_3DRANDOMWALK_POLYGON_H
