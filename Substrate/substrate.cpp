//
// Created by Alemany Juvanteny, Ignasi on 26/02/2023.
//

#include "substrate.h"
#include <chrono>
#include <limits>
#include <boost/variant.hpp>
substrate::substrate(std::vector<Eigen::MatrixXd> &myo_vertices, std::vector<Eigen::MatrixXd> &myo_faces) : _tree()
{
    auto start = std::chrono::high_resolution_clock::now(); // start timer

    _myocytes.reserve(myo_vertices.size());
    _points.reserve(myo_vertices.size());

    int i = 0;
    for (auto &myo_vertices_i : myo_vertices)
    {
        polygon myo(myo_vertices_i, myo_faces[i]);
        Polyhedron myo_poly = myo.getPolyhedron();
        CGAL::Bbox_3 bbox = myo.getSolidBbox();
        Kernel::Point_3 centroid = CGAL::midpoint(
            Kernel::Point_3(bbox.xmin(), bbox.ymin(), bbox.zmin()),
            Kernel::Point_3(bbox.xmax(), bbox.ymax(), bbox.zmax()));
        _myocytes.emplace_back(std::move(myo));
        _points.emplace_back(centroid);
        _map_centroid_to_polygon[centroid] = i;
        ++i;
    }

    _tree = std::make_unique<Tree>(_points.begin(), _points.end());
    auto end = std::chrono::high_resolution_clock::now(); // end timer

    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Elapsed time set substrate: " << elapsed.count() << " seconds" << std::endl;
}

void substrate::setTransform(transform &t)
{
    _transform = t;
}

boost::optional<std::tuple<int, double, Eigen::Vector3d>> substrate::intersectPolygon(const Eigen::Vector3d &point, const Eigen::Vector3d &step) const
{
    // We will place the point in the center of the block as the tree is built with the centroid of the polygons

    // TODO: This will have to be modified when considering a full geometry. Maybe input into to the Ktree the 8 points of bounding box?
    Kernel::Point_3 query((double)point(0), (double)point(1), 126.0 / 2);

    Neighbor_search search(*_tree, query, 3);

    double min_distance_children = std::numeric_limits<double>::max();
    int index_polygon = -1;
    Eigen::Vector3d normal_vector = Eigen::Vector3d::Zero();

    for (Neighbor_search::iterator it = search.begin(); it != search.end(); ++it)
    {
        boost::optional<std::pair<int, double>> intersection = _myocytes[_map_centroid_to_polygon.at(it->first)].intersection(point, step);
        if (intersection)
        {
            int index_face = boost::get<std::pair<int, double>>(*intersection).first;
            double distance_to_face = boost::get<std::pair<int, double>>(*intersection).second;

            if (distance_to_face < min_distance_children)
            {
                min_distance_children = distance_to_face;
                index_polygon = _map_centroid_to_polygon.at(it->first);
                normal_vector = _myocytes[_map_centroid_to_polygon.at(it->first)].getNormalVector(index_face);
            }
        }
    }

    if (min_distance_children == std::numeric_limits<double>::max())
    {
        return boost::optional<std::tuple<int, double, Eigen::Vector3d>>();
    }

    return std::make_tuple(index_polygon, min_distance_children, normal_vector);
}

boost::optional<double> substrate::intersectionBlock(const Eigen::Vector3d &point, const Eigen::Vector3d &step) const
{
    // TODO: If we want reflecting boundaries in the block we need to output the normal vector of the face intersected
    Polyhedron box = _transform.get_block();
    CGAL::Bbox_3 solid_box = _transform.get_solid_block();
    auto box_tree = _transform.get_AABBtree();
    CGAL::Point_3<Kernel> p1(point(0), point(1), point(2));
    CGAL::Point_3<Kernel> p2(point(0) + step(0), point(1) + step(1), point(2) + step(2));
    CGAL::Segment_3<Kernel> segment(p1, p2);

    // Check if point p2 is inside the solid_box without the boundaries
    if (p2.x() >= solid_box.xmin() && p2.x() <= solid_box.xmax() &&
        p2.y() >= solid_box.ymin() && p2.y() <= solid_box.ymax() &&
        p2.z() >= solid_box.zmin() && p2.z() <= solid_box.zmax())
    {
        return boost::optional<double>();
    }

    std::vector<std::pair<boost::variant<Kernel::Point_3, Kernel::Segment_3>, Primitive>> intersections;
    box_tree->all_intersections(segment, std::back_inserter(intersections));
    double min_distance = std::numeric_limits<double>::max();

    // No intersection found
    if (intersections.size() == 0)
    {
        throw std::runtime_error("Substrate::intersectionBlock -> No intersection found but segment.target() is outside block");
        // std::cout << "No intersection found but point p2 is outside the solid_box" << std::endl;
        // return boost::optional<double>();
    }

    for (const auto &intersection : intersections)
    {
        // If intersection is point
        if (intersection.first.which() == 0)
        {
            Kernel::Point_3 point_intersected = boost::get<Kernel::Point_3>(intersection.first);
            double distance = CGAL::sqrt(CGAL::squared_distance(point_intersected, segment.source()));
            if (distance < min_distance)
            {
                min_distance = distance;
            }
        }
        else
        {
            throw std::runtime_error("Substrate::intersectionBlock -> Intersection found but segment is parallel and lies on face");
        }
    }
    // TODO: Unify epsilon values in the code
    double t = min_distance / CGAL::sqrt(segment.squared_length());
    if ((1 - t) < 1e-6)
    {
        throw std::runtime_error("Substrate::intersectionBlock -> Remaining step is too short, might be uncertain");
    }

    return min_distance;
}

int substrate::searchPolygon(const Eigen::Vector3d &point, const std::string &frameOfReference) const
{
    Eigen::Vector3d point_local;
    if (frameOfReference == "global")
    {
        // TODO: Check if point is out of block and if it is check that point_local is. otherwise throw error means that transform is not init properly
        transform_info transform_data = getLocalFromGlobal(point);
        point_local = transform_data.local_position;
    }
    else
    {
        point_local = point;
    }

    // TODO: Query the whole point when not using the block. Take into account when getting rid of block or implementing the transform type
    //  We will place the point in the center of the block as the tree is built with the centroid of the polygons
    Kernel::Point_3 query((double)point_local(0), (double)point_local(1), 126.0 / 2);

    Neighbor_search search(*_tree, query, 3);

    for (Neighbor_search::iterator it = search.begin(); it != search.end(); ++it)
    {
        int index_polygon = _map_centroid_to_polygon.at(it->first);
        if (_myocytes[index_polygon].containsPoint(point_local))
        {
            return index_polygon;
        }
    }
    return -1;
}
