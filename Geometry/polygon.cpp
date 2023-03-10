//
// Created by Alemany Juvanteny, Ignasi on 27/02/2023.
//

#include "polygon.h"

polygon::polygon(const Eigen::MatrixXd &vertices, const Eigen::MatrixXd &faces)
{
    // Create a bounding box from min max vertices

    for (int i = 0; i < vertices.rows(); i++)
    {
        _vertices.push_back(CGAL::Point_3<Kernel>((double)vertices(i, 0), (double)vertices(i, 1), (double)vertices(i, 2)));
    }

    _bbox = CGAL::bbox_3(_vertices.begin(), _vertices.end());

    // Initialize the polyhedron builder
    CGAL::Polyhedron_incremental_builder_3<HalfedgeDS> B(_poly.hds(), true);
    B.begin_surface(_vertices.size(), faces.rows());

    // Add vertices
    for (const auto &point : _vertices)
    {
        B.add_vertex(point);
    }

    // Add faces and triangle faces
    for (int i = 0; i < faces.rows(); i++)
    {
        std::vector<std::size_t> face_input;
        for (int j = 0; j < 3; j++)
        {
            face_input.push_back(faces(i, j) - 1);
        }
        B.begin_facet();
        for (const auto &index : face_input)
        {
            B.add_vertex_to_facet(index);
        }
        B.end_facet();

        Kernel::Triangle_3 triangle_face(_vertices[face_input[0]], _vertices[face_input[1]], _vertices[face_input[2]]);
        triangle_faces.push_back(triangle_face);
    }
    B.end_surface();

    _AABBtree = std::make_unique<Tree_AABB>(triangle_faces.begin(), triangle_faces.end());
}

double polygon::computeVolume()
{
    double volume = 0;

    Polyhedron::Point centroid = CGAL::centroid(_vertices.begin(), _vertices.end());
    // Compute the volume of each tetrahedron and add it to the volume
    for (Polyhedron::Facet_iterator facet_it = _poly.facets_begin(); facet_it != _poly.facets_end(); ++facet_it)
    {
        Polyhedron::Halfedge_around_facet_circulator he_circ = facet_it->facet_begin();
        Polyhedron::Point_3 p1 = he_circ->vertex()->point();
        Polyhedron::Point_3 p2 = he_circ->next()->vertex()->point();
        Polyhedron::Point_3 p3 = he_circ->next()->next()->vertex()->point();

        // Compute the pyramid volume of p1,p2,p3,centroid
        volume += CGAL::abs(CGAL::volume(p1, p2, p3, centroid));
    }
    return volume;
}

double polygon::computeSurface()
{
    double surface_area = 0;

    // Compute the area of each triangle and add it to the surface area
    for (Polyhedron::Facet_iterator facet_it = _poly.facets_begin(); facet_it != _poly.facets_end(); ++facet_it)
    {
        Polyhedron::Halfedge_around_facet_circulator he_circ = facet_it->facet_begin();
        Polyhedron::Point_3 p1 = he_circ->vertex()->point();
        Polyhedron::Point_3 p2 = he_circ->next()->vertex()->point();
        Polyhedron::Point_3 p3 = he_circ->next()->next()->vertex()->point();

        surface_area += CGAL::sqrt(CGAL::squared_area(p1, p2, p3));
    }
    return surface_area;
}
/*
boost::variant<bool, std::pair<int, double>> polygon::intersection(const Eigen::Vector3d &point, const Eigen::Vector3d &step)
{
    // Define a segment that starts at point and finishes at step+point;
    CGAL::Segment_3<Kernel> segment(CGAL::Point_3<Kernel>(point(0), point(1), point(2)), CGAL::Point_3<Kernel>(step(0) + point(0), step(1) + point(1), step(2) + point(2)));

    double min_distance = std::numeric_limits<double>::max();
    int min_index = -1;
    int n_intersections = 0;

    // Check if the segment intersects any of the faces
    for (Polyhedron::Facet_iterator facet_it = _poly.facets_begin(); facet_it != _poly.facets_end(); ++facet_it)
    {
        Polyhedron::Halfedge_around_facet_circulator he_circ = facet_it->facet_begin();
        Polyhedron::Point_3 p1 = he_circ->vertex()->point();
        Polyhedron::Point_3 p2 = he_circ->next()->vertex()->point();
        Polyhedron::Point_3 p3 = he_circ->next()->next()->vertex()->point();
        CGAL::Triangle_3<Kernel> triangle(p1, p2, p3);
        boost::optional<boost::variant<CGAL::Point_3<Kernel>, CGAL::Segment_3<Kernel>>> intersection = CGAL::intersection(segment, triangle);

        // If init point is right at the face
        if (triangle.has_on(segment.source()))
        {
            continue;
        }

        // Check if intersection is type Segment_3 (is intersecting with edge) and return false as intersection is not certain
        if (intersection && intersection->which() == 1)
        {
            return false;
        }

        // Check if intersection is type Point_3 and update min_distance and min_index
        if (intersection && intersection->which() == 0)
        {
            n_intersections++;
            CGAL::Point_3<Kernel> intersection_point;
            if (CGAL::assign(intersection_point, *intersection))
            {
                double distance = CGAL::sqrt(CGAL::squared_distance(intersection_point, segment.source()));

                // If the distance is 0, the intersection is the starting point of the segment, so return false
                if (distance < 1e-8)
                {
                    return false;
                }

                // If the remaining step is also smaller than 1e-1, return false
                if (CGAL::sqrt(segment.squared_length()) - distance < 1e-8)
                {
                    return false;
                }

                if (distance < min_distance)
                {
                    min_distance = distance;
                    min_index = std::distance(_poly.facets_begin(), facet_it) + 1;
                }
            }
        }
    }

    if (n_intersections > 1)
    {
        // TODO: Remove this warning if it is not needed
        std::cout << "WARNING: More than one intersection found." << std::endl;
    }

    // If no intersection is found, return false
    if (min_distance == std::numeric_limits<double>::max())
    {
        return false;
    }

    // If an intersection is found, return the index of the face and the distance to the intersection
    return std::pair<int, double>(min_index, min_distance);
}
*/
boost::variant<bool, std::pair<int, double>> polygon::intersection(const Eigen::Vector3d &point, const Eigen::Vector3d &step)
{
    Kernel::Segment_3 segment(Kernel::Point_3(point(0), point(1), point(2)), Kernel::Point_3(step(0) + point(0), step(1) + point(1), step(2) + point(2)));
    std::vector<std::pair<boost::variant<Kernel::Point_3, Kernel::Segment_3>, Primitive>> intersections;
    _AABBtree->all_intersections(segment, std::back_inserter(intersections));
    double min_distance = std::numeric_limits<double>::max();
    int min_index = -1;

    // No intersection found
    if (intersections.size() == 0)
    {
        return false;
    }

    for (const auto &intersection : intersections)
    {
        if (intersection.first.which() == 0)
        {
            Kernel::Point_3 point_intersected = boost::get<Kernel::Point_3>(intersection.first);
            double distance = CGAL::sqrt(CGAL::squared_distance(point_intersected, segment.source()));
            if (distance == 0)
            {
                continue;
            }
            if (distance < min_distance)
            {
                auto iter = intersection.second.id();
                min_index = std::distance(triangle_faces.begin(), iter) + 1;
                min_distance = distance;
            }
        }
        else
        {
            // Handle intersection with segment, return false as it's not certain
            return false;
        }
    }

    // If there is only one intersection, but it is right at the face min_distance is not updated
    if (min_distance == std::numeric_limits<double>::max())
    {
        return false;
    }
    // If is too close to face or the remaining step is too short
    if (min_distance < 1e-8 || CGAL::sqrt(segment.squared_length()) - min_distance < 1e-8)
    {
        return false;
    }

    return std::pair<int, double>(min_index, min_distance);
}

/*
bool polygon::containsPoint(const Eigen::Vector3d &point)
{
    // Check if point is within min and max of the bounding box
    if (_bbox.xmin() > point(0) || _bbox.xmax() < point(0) || _bbox.ymin() > point(1) || _bbox.ymax() < point(1) || _bbox.zmin() > point(2) || _bbox.zmax() < point(2))
    {
        return false;
    }

    // Define a ray that starts at the point and goes in the positive x direction
    CGAL::Ray_3<Kernel> ray(CGAL::Point_3<Kernel>(point(0), point(1), point(2)), CGAL::Vector_3<Kernel>(1, 0, 0));
    int intersection_count = 0;

    // Check if the point lies on any of the faces
    CGAL::Point_3<Kernel> origin(point(0), point(1), point(2));
    for (Polyhedron::Facet_iterator facet_it = _poly.facets_begin(); facet_it != _poly.facets_end(); ++facet_it)
    {
        Polyhedron::Halfedge_around_facet_circulator he_circ = facet_it->facet_begin();
        Polyhedron::Point_3 p1 = he_circ->vertex()->point();
        Polyhedron::Point_3 p2 = he_circ->next()->vertex()->point();
        Polyhedron::Point_3 p3 = he_circ->next()->next()->vertex()->point();
        CGAL::Triangle_3<Kernel> triangle(p1, p2, p3);
        if (triangle.has_on(origin))
        {
            return true;
        }
    }

    for (Polyhedron::Facet_iterator facet_it = _poly.facets_begin(); facet_it != _poly.facets_end(); ++facet_it)
    {
        Polyhedron::Halfedge_around_facet_circulator he_circ = facet_it->facet_begin();
        Polyhedron::Point_3 p1 = he_circ->vertex()->point();
        Polyhedron::Point_3 p2 = he_circ->next()->vertex()->point();
        Polyhedron::Point_3 p3 = he_circ->next()->next()->vertex()->point();
        CGAL::Triangle_3<Kernel> triangle(p1, p2, p3);
        boost::optional<boost::variant<CGAL::Point_3<Kernel>, CGAL::Segment_3<Kernel>>> intersection = CGAL::intersection(ray, triangle);
        if (intersection)
        {
            intersection_count++;
        }
    }
    return (intersection_count % 2) == 1;
}
*/
bool polygon::containsPoint(const Eigen::Vector3d &point)
{
    if (_bbox.xmin() > point(0) || _bbox.xmax() < point(0) || _bbox.ymin() > point(1) ||
        _bbox.ymax() < point(1) || _bbox.zmin() > point(2) || _bbox.zmax() < point(2))
    {
        return false;
    }

    // Define a ray that starts at the point and goes in the positive x direction
    CGAL::Ray_3<Kernel> ray(CGAL::Point_3<Kernel>(point(0), point(1), point(2)), CGAL::Vector_3<Kernel>(1, 0, 0));

    // Find all intersections between the ray and the triangles in the AABB tree
    std::vector<std::pair<boost::variant<Kernel::Point_3, Kernel::Segment_3>, Primitive>> intersections;
    _AABBtree->all_intersections(ray, std::back_inserter(intersections));

    for (const auto &intersection : intersections)
    {
        if (intersection.first.which() == 0)
        {
            Kernel::Point_3 point_intersected = boost::get<Kernel::Point_3>(intersection.first);
            Primitive intersection_primitive = intersection.second;
            auto iter = intersection_primitive.id();
            int index = iter - triangle_faces.begin();
            if (triangle_faces[index].has_on(point_intersected))
            {
                return true;
            }
        }
    }

    // Check if the number of intersections is odd
    return (intersections.size() % 2) == 1;
}
