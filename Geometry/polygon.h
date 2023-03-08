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
#include <cassert>
#include <vector>
#include <CGAL/Ray_3.h>
#include <CGAL/Triangle_3.h>
#include <CGAL/Vector_3.h>
#include <CGAL/intersection_3.h>
#include <CGAL/Bbox_3.h>
#include <boost/optional.hpp>

typedef CGAL::Simple_cartesian<double>     Kernel;
typedef CGAL::Polyhedron_3<Kernel>         Polyhedron;
typedef Polyhedron::HalfedgeDS             HalfedgeDS;
typedef Kernel::Point_3                    Point_3;
class polygon
{
public:

    polygon() = default;

    polygon(const Eigen::MatrixXd &vertices_input, const Eigen::MatrixXd &faces_input);

    //double computeVolume();
    double computeSurface();

    boost::variant<bool,std::pair<int,double>> intersection(const Eigen::Vector3d &point, const Eigen::Vector3d &step);

    bool containsPoint(const Eigen::Vector3d &point);

private:
    Polyhedron _poly;
    std::vector<Point_3> _vertices;
    std::vector<std::vector<std::size_t>> _faces;
    CGAL::Bbox_3 _bbox;
};

#endif // INC_3DRANDOMWALK_POLYGON_H
