//
// Created by Alemany Juvanteny, Ignasi on 27/02/2023.
//

#include "polygon.h"

polygon::polygon(const Eigen::MatrixXd &vertices, const Eigen::MatrixXd &faces)
{
    // Create a bounding box from min max vertices

    for (int i = 0; i < vertices.rows(); i++)
    {
        _vertices.push_back(Point_3(vertices(i, 0), vertices(i, 1), vertices(i, 2)));
    }

    _bbox = CGAL::bbox_3(_vertices.begin(), _vertices.end());

    for (int i = 0; i < faces.rows(); i++)
    {
        std::vector<std::size_t> face_input;
        for (int j = 0; j < 3; j++)
        {
            face_input.push_back(faces(i, j) - 1);
        }
        _faces.push_back(face_input);
    }
    // Postcondition: hds is a valid polyhedral surface.
    CGAL::Polyhedron_incremental_builder_3<HalfedgeDS> B(_poly.hds(), true);
    B.begin_surface(_vertices.size(), _faces.size());

    // Add vertices
    for (const auto &point : _vertices)
    {
        B.add_vertex(point);
    }

    // Add faces
    for (const auto &face : _faces)
    {
        B.begin_facet();
        for (const auto &index : face)
        {
            B.add_vertex_to_facet(index);
        }
        B.end_facet();
    }

    B.end_surface();
}

double polygon::computeVolume(){
    double volume=0;

    Polyhedron::Point_3 centroid = CGAL::centroid(_vertices.begin(), _vertices.end());
    // Compute the volume of each tetrahedron and add it to the volume
    for (Polyhedron::Facet_iterator facet_it = _poly.facets_begin(); facet_it != _poly.facets_end(); ++facet_it)
    {
        Polyhedron::Halfedge_around_facet_circulator he_circ = facet_it->facet_begin();
        Polyhedron::Point_3 p1 = he_circ->vertex()->point();
        Polyhedron::Point_3 p2 = he_circ->next()->vertex()->point();
        Polyhedron::Point_3 p3 = he_circ->next()->next()->vertex()->point();

        //Compute the pyramid volume of p1,p2,p3,centroid
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

boost::variant<bool, std::pair<int, double>> polygon::intersection(const Eigen::Vector3d &point, const Eigen::Vector3d &step)
{
    // Define a segment that starts at point and finishes at step+point;
    CGAL::Segment_3<Kernel> segment(Kernel::Point_3(point(0), point(1), point(2)), Kernel::Point_3(step(0) + point(0), step(1) + point(1), step(2) + point(2)));

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

bool polygon::containsPoint(const Eigen::Vector3d &point)
{
    // Check if point is within min and max of the bounding box
    if (_bbox.xmin() > point(0) || _bbox.xmax() < point(0) || _bbox.ymin() > point(1) || _bbox.ymax() < point(1) || _bbox.zmin() > point(2) || _bbox.zmax() < point(2))
    {
        return false;
    }

    // Define a ray that starts at the point and goes in the positive x direction
    CGAL::Ray_3<Kernel> ray(Kernel::Point_3(point(0), point(1), point(2)), Kernel::Vector_3(1, 0, 0));
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
        if (CGAL::coplanar(p1, p2, p3, origin) && triangle.has_on(origin))
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
