//
// Created by Alemany Juvanteny, Ignasi on 26/02/2023.
//

#ifndef INC_3DRANDOMWALK_SUBSTRATE_H
#define INC_3DRANDOMWALK_SUBSTRATE_H

#include "transform.h"
#include "../Geometry/polygon.h"
#include <list>
#include <vector>
#include <map>
#include <list>
#include <CGAL/point_generators_3.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <memory>

typedef CGAL::Search_traits_3<Kernel> TreeTraits;
typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
typedef Neighbor_search::Tree Tree;

class substrate
{

public:
    substrate() = default;
    substrate(std::vector<Eigen::MatrixXd> &myo_vertices, std::vector<Eigen::MatrixXd> &myo_faces);
    void setTransform(transform &t);

    transform_info getLocalFromGlobal(const Eigen::Vector3d &global_position) const { return _transform.global2local(global_position); };
    Eigen::Vector3d getGlobalFromLocal(const Eigen::Vector3d &local_position, int iX, int iY, int iZ) const { return _transform.local2global(local_position, iX, iY, iZ); };

    bool containsPoint(int index_polygon, const Eigen::Vector3d &point) const { return _myocytes[index_polygon].containsPoint(point); };
    boost::optional<std::tuple<int, double, Eigen::Vector3d>> intersectPolygon(const Eigen::Vector3d &point, const Eigen::Vector3d &step) const;
    boost::optional<double> intersectionBlock(const Eigen::Vector3d &point, const Eigen::Vector3d &step) const;
    int searchPolygon(const Eigen::Vector3d &point, const std::string &frameOfReference = "local") const;

private:
    std::vector<polygon> _myocytes;
    std::vector<Kernel::Point_3> _points;
    transform _transform;
    std::unique_ptr<Tree> _tree;
    std::map<Kernel::Point_3, int> _map_centroid_to_polygon;
};

#endif // INC_3DRANDOMWALK_SUBSTRATE_H
