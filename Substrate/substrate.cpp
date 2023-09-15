//
// Created by Alemany Juvanteny, Ignasi on 26/02/2023.
//

#include "substrate.h"
#include <chrono>
#include <limits>
#include <boost/variant.hpp>
#include <boost/filesystem.hpp>
substrate::substrate(std::vector<Eigen::MatrixXd> &myo_vertices, std::vector<Eigen::MatrixXd> &myo_faces) : _tree()
{
    auto start = std::chrono::high_resolution_clock::now(); // start timer

    _myocytes.reserve(myo_vertices.size());
    _points.reserve(myo_vertices.size());

    int i = 0;
    for (auto &myo_vertices_i : myo_vertices)
    {
        polygon myo(myo_vertices_i, myo_faces[i]);
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
substrate::substrate(const std::vector<std::string>& polygons){

    //iterate over polygons
    int i=0;
    auto start = std::chrono::high_resolution_clock::now(); // start timer

    for (const auto &poly: polygons){
        polygon myo(poly);
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

void substrate::setStrain(std::function<double(double)> strain_f){
    strain = strain_f;
}


void substrate::write_ply_surface_mesh(const Mesh &mesh, const std::string &filename)
{

    std::ofstream out(filename);
    if (!out.is_open())
    {
        throw std::runtime_error("Could not open file " + filename);
        return;
    }

    out << "ply\n";
    out << "format ascii 1.0\n";
    out << "element vertex " << mesh.number_of_vertices() << "\n";
    out << "property float x\n";
    out << "property float y\n";
    out << "property float z\n";
    out << "element face " << mesh.number_of_faces() << "\n";
    out << "property list uchar int vertex_indices\n";
    out << "end_header\n";

    // Write vertices
    for (const auto &v : mesh.vertices())
    {
        const auto &point = mesh.point(v);
        out << point.x() << " " << point.y() << " " << point.z() << "\n";
    }

    // Write faces
    for (const auto &f : mesh.faces())
    {
        out << "3";
        CGAL::Vertex_around_face_circulator<Mesh> vcirc(mesh.halfedge(f), mesh), done(vcirc);
        do
        {
            out << " " << (*vcirc).idx();
        } while (++vcirc != done);
        out << "\n";
    }

    out.close();
}

void substrate::write_ply_polyhedron(polygon &poly, const std::string &filename)
{
    std::ofstream out(filename);
    if (!out.is_open())
    {
        throw std::runtime_error("Could not open file " + filename);
        return;
    }

    // Write header
    out << "ply\n";
    out << "format ascii 1.0\n";
    out << "element vertex " << poly.getVertices().size() << "\n";
    out << "property float x\n";
    out << "property float y\n";
    out << "property float z\n";
    out << "element face " << poly.getFaces().size() << "\n";
    out << "property list uchar int vertex_indices\n";
    out << "end_header\n";

    // Write vertices
    for (auto &vertex : poly.getVertices())
    {
        out << vertex.x() << " " << vertex.y() << " " << vertex.z() << "\n";
    }

    for (auto &face : poly.getFaces())
    {
        out << "3 " << face[0] << " " << face[1] << " " << face[2] << "\n";
    }
    out.close();
}

void substrate::setTransform(transform &t)
{
    _transform = t;
}

boost::optional<std::tuple<int, double, Eigen::Vector3d>> substrate::intersectPolygon(const Eigen::Vector3d &point, const Eigen::Vector3d &step, const double time)
{
    // We will place the point in the center of the block as the tree is built with the centroid of the polygons

    // TODO: This will have to be modified when considering a full geometry. Maybe input into to the Ktree the 8 points of bounding box?
    Kernel::Point_3 query((double)point(0), (double)point(1), 126.0 / 2);

    Neighbor_search search(*_tree, query, 5);

    double min_distance_children = std::numeric_limits<double>::max();
    int index_polygon = -1;
    Eigen::Vector3d normal_vector = Eigen::Vector3d::Zero();

    for (Neighbor_search::iterator it = search.begin(); it != search.end(); ++it)
    {
        //Pass in total_time here 
        boost::optional<std::pair<int, double>> intersection = _myocytes[_map_centroid_to_polygon.at(it->first)].intersection(point, step, time, strain);
        if (intersection)
        {
            int index_face = boost::get<std::pair<int, double>>(*intersection).first;
            Kernel::FT distance_to_face = boost::get<std::pair<int, double>>(*intersection).second;

            if (distance_to_face < min_distance_children)
            {
                min_distance_children = CGAL::to_double(distance_to_face);
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
    }

    for (const auto &intersection : intersections)
    {
        // If intersection is point
        if (intersection.first.which() == 0)
        {
            Kernel::Point_3 point_intersected = boost::get<Kernel::Point_3>(intersection.first);
            auto primitive = boost::get<Primitive>(intersection.second);
            Kernel::Triangle_3 triangle = primitive.datum();
            // check point_intersected inside triangle
            // if (triangle.has_on(point_intersected))
            //{
            Kernel::Segment_3 edge1(triangle.vertex(1), triangle.vertex(0));
            Kernel::Segment_3 edge2(triangle.vertex(2), triangle.vertex(0));
            Kernel::Segment_3 edge3(triangle.vertex(2), triangle.vertex(1));

            double min_edge_distance = std::min(CGAL::sqrt(CGAL::to_double(CGAL::squared_distance(point_intersected, edge1))), std::min(
                                                                                                                                   CGAL::sqrt(CGAL::to_double(CGAL::squared_distance(point_intersected, edge2))),
                                                                                                                                   CGAL::sqrt(CGAL::to_double(CGAL::squared_distance(point_intersected, edge3)))));

            // Check if the point is too close to an edge
            if (min_edge_distance < 1e-8)
            {
                throw std::runtime_error("Polygon::intersection -> Intersection found but point is too close to an edge, uncertain");
            }

            Kernel::FT squared_distance = CGAL::squared_distance(point_intersected, segment.source());
            double distance = CGAL::sqrt(CGAL::to_double(squared_distance));
            if (distance < min_distance)
            {
                min_distance = distance;
            }
            //}
            // else
            //{
            //   throw std::runtime_error("Polygon::intersection -> Intersection found but point is not on face, uncertain");
            //}
        }
        else
        {
            throw std::runtime_error("Substrate::intersectionBlock -> Intersection found but segment is parallel and lies on face");
        }
    }

    if (min_distance == std::numeric_limits<double>::max())
    {
        throw std::runtime_error("Substrate::intersectionBlock -> Intersection error.");
    }
    // TODO: Unify epsilon values in the code
    double t = CGAL::to_double(min_distance) / step.norm();
    if ((1 - t) < 1e-8 || t < 1e-8)
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

    Neighbor_search search(*_tree, query, 10);

    std::vector<std::pair<int, double>> indices_and_squared_distances;

    for (Neighbor_search::iterator it = search.begin(); it != search.end(); ++it)
    {
        int index_polygon = _map_centroid_to_polygon.at(it->first);
        double squared_distance = CGAL::to_double(it->second); // Get the squared distance directly
        indices_and_squared_distances.push_back(std::make_pair(index_polygon, squared_distance));
    }

    // Sort the vector by squared distance
    std::sort(indices_and_squared_distances.begin(), indices_and_squared_distances.end(),
              [](const std::pair<int, double> &a, const std::pair<int, double> &b)
              {
                  return a.second < b.second;
              });

    // Check if the point is contained in any of the sorted polygons
    for (const auto &index_squared_distance_pair : indices_and_squared_distances)
    {
        int index_polygon = index_squared_distance_pair.first;
        if (_myocytes[index_polygon].containsPoint(point_local))
        {
            return index_polygon;
        }
    }

    return -1;
}
