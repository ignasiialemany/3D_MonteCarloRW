//
// Created by Alemany Juvanteny, Ignasi on 26/02/2023.
//

#include "substrate.h"
#include <chrono>
#include <limits>
#include <boost/variant.hpp>
substrate::substrate(std::vector<Eigen::MatrixXd> &myo_vertices, std::vector<Eigen::MatrixXd> &myo_faces) : _tree()
{
    auto start = std::chrono::high_resolution_clock::now();  // start timer

    _myocytes.reserve(myo_vertices.size());
    _points.reserve(myo_vertices.size());
    
    int i = 0;
    for (auto& myo_vertices_i : myo_vertices)
    {
        polygon myo(myo_vertices_i, myo_faces[i]);
        Polyhedron myo_poly = myo.getPolyhedron();
        CGAL::Bbox_3 bbox = myo.getBbox();
        Kernel::Point_3 centroid = CGAL::midpoint(
            Kernel::Point_3(bbox.xmin(), bbox.ymin(), bbox.zmin()),
            Kernel::Point_3(bbox.xmax(), bbox.ymax(), bbox.zmax()));
        _myocytes.emplace_back(std::move(myo));
        _points.emplace_back(centroid);
        _map_centroid_to_polygon[centroid] = i;
        ++i;
    }

    _tree = std::make_unique<Tree>(_points.begin(), _points.end());
                auto end = std::chrono::high_resolution_clock::now();  // end timer

    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Elapsed time set substrate: " << elapsed.count() << " seconds" << std::endl;
}

void substrate::setTransform(transform &t)
{
    _transform = t;
}

boost::variant<bool, std::pair<int, double>> substrate::intersection(const Eigen::Vector3d &point, const Eigen::Vector3d &step) const
{
    //We will place the point in the center of the block as the tree is built with the centroid of the polygons
    Kernel::Point_3 query((double)point(0), (double)point(1), 126.0/2);

    Neighbor_search search(*_tree, query, 3);

    double min_distance_children = std::numeric_limits<double>::max();
    int index_polygon = -1;

    for (Neighbor_search::iterator it = search.begin(); it != search.end(); ++it)
    {
        boost::variant<bool,std::pair<int,double>> intersection = _myocytes[_map_centroid_to_polygon.at(it->first)].intersection(point, step);
        if (intersection.which() == 1)
        {
            //auto end = std::chrono::high_resolution_clock::now();  // end timer
            //std::chrono::duration<double> elapsed = end - start;
            //std::cout << "Elapsed time to intersect myocyte: " << elapsed.count() << " seconds" << std::endl;
            if (boost::get<std::pair<int,double>>(intersection).second < min_distance_children)
            {
                min_distance_children = boost::get<std::pair<int,double>>(intersection).second;
                index_polygon = _map_centroid_to_polygon.at(it->first);
            }
        }
    }
    
    //auto end = std::chrono::high_resolution_clock::now();  // end timer
    //std::chrono::duration<double> elapsed = end - start;
    //std::cout << "Elapsed time to intersect myocyte: " << elapsed.count() << " seconds" << std::endl;
    
    if (min_distance_children == std::numeric_limits<double>::max())
    {
        return false;
    }
    else
    {
        return std::make_pair(index_polygon, min_distance_children);
    }
   
    return false;
}

int substrate::searchPolygon(Eigen::Vector3d &point)
{
    //auto start = std::chrono::high_resolution_clock::now();  // start timer

    //We will place the point in the center of the block as the tree is built with the centroid of the polygons
    Kernel::Point_3 query((double)point(0), (double)point(1), 126.0/2);

    Neighbor_search search(*_tree, query, 3);

    for (Neighbor_search::iterator it = search.begin(); it != search.end(); ++it)
    {
        int index_polygon = _map_centroid_to_polygon[it->first];
        if (_myocytes[index_polygon].containsPoint(point))
        {
            //auto end = std::chrono::high_resolution_clock::now();  // end timer
            //std::chrono::duration<double> elapsed = end - start;
            //std::cout << "Elapsed time: " << elapsed.count() << " seconds" << std::endl;
            return index_polygon;
        }
    }
    
    //auto end = std::chrono::high_resolution_clock::now();  // end timer
    //std::chrono::duration<double> elapsed = end - start;
    //std::cout << "Elapsed time to search myocyte: " << elapsed.count() << " seconds" << std::endl;
    
    return -1;
}
