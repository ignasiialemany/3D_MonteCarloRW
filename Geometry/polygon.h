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
#include <CGAL/Object.h>
#include <memory>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef Polyhedron::HalfedgeDS HalfedgeDS;

// Define the AABB traits
typedef std::vector<Kernel::Triangle_3>::iterator Iterator;
typedef CGAL::AABB_triangle_primitive<CGAL::Simple_cartesian<double>, Iterator> Primitive;
typedef CGAL::AABB_traits<CGAL::Simple_cartesian<double>, Primitive> AABB_triangle_traits;
typedef CGAL::AABB_tree<AABB_triangle_traits> Tree_AABB;

class polygon
{
public:
    polygon() = default;

    polygon(const Eigen::MatrixXd &vertices_input, const Eigen::MatrixXd &faces_input);

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


private:
    // TODO: Might delete some of these variables, probably _vertices and _faces
    Polyhedron _poly;
    Polyhedron _bbox;
    CGAL::Bbox_3 _solid_bbox;
    void createBbox(Kernel::Point_3 min_point, Kernel::Point_3 max_point);
    void createPolygon(const Eigen::MatrixXd &vertices, const Eigen::MatrixXd &faces);

    std::vector<CGAL::Point_3<Kernel>> _vertices;
    std::vector<std::vector<std::size_t>> _faces;
    std::unique_ptr<Tree_AABB> _AABBtree;
    std::vector<Kernel::Triangle_3> triangle_faces;

};

#endif // INC_3DRANDOMWALK_POLYGON_H
