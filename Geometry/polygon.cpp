//
// Created by Alemany Juvanteny, Ignasi on 27/02/2023.
//

#include "polygon.h"
#include <CGAL/enum.h>

polygon::polygon(const Eigen::MatrixXd &vertices, const Eigen::MatrixXd &faces)
{
    createPolygon(vertices, faces);
    if (!_poly.is_valid() || !_poly.is_closed())
    {
        std::cerr << "Error: The input mesh is not valid or not closed." << std::endl;
        return;
    }
    // Fix non-manifold issues
    CGAL::Polygon_mesh_processing::duplicate_non_manifold_vertices(_poly);
    CGAL::Polygon_mesh_processing::stitch_borders(_poly);
    _AABBtree = std::make_unique<Tree_AABB>(triangle_faces.begin(), triangle_faces.end());
    _solid_bbox = _AABBtree->bbox();
    Kernel::Point_3 centroid_poly = CGAL::centroid(_vertices.begin(), _vertices.end());
    centroid = Eigen::Vector3d(centroid_poly.x(), centroid_poly.y(), centroid_poly.z());
}


polygon::polygon(Polyhedron &poly)
{
    _poly = poly;
    // Iterate over the Polyhedron_3 and push back Triangle_3 to triangle_faces
    for (Polyhedron::Facet_iterator facet_it = _poly.facets_begin(); facet_it != _poly.facets_end(); ++facet_it)
    {
        Polyhedron::Halfedge_around_facet_circulator he_circ = facet_it->facet_begin();
        Polyhedron::Point_3 p1 = he_circ->vertex()->point();
        Polyhedron::Point_3 p2 = he_circ->next()->vertex()->point();
        Polyhedron::Point_3 p3 = he_circ->next()->next()->vertex()->point();
        // create std::vector<std::size_t> for each face
        triangle_faces.push_back(Kernel::Triangle_3(p1, p2, p3));
    }
    // Iterate over vertices and get min and max vertices
    _AABBtree = std::make_unique<Tree_AABB>(triangle_faces.begin(), triangle_faces.end());
    _solid_bbox = _AABBtree->bbox();
    
    // Fill in _vertices and _faces
    for (Polyhedron::Vertex_iterator vertex_it = _poly.vertices_begin(); vertex_it != _poly.vertices_end(); ++vertex_it)
    {
        _vertices.push_back(vertex_it->point());
    }

    for (Polyhedron::Facet_iterator facet_it = _poly.facets_begin(); facet_it != _poly.facets_end(); ++facet_it)
    {
        Polyhedron::Halfedge_around_facet_circulator he_circ = facet_it->facet_begin();
        Polyhedron::Point_3 p1 = he_circ->vertex()->point();
        Polyhedron::Point_3 p2 = he_circ->next()->vertex()->point();
        Polyhedron::Point_3 p3 = he_circ->next()->next()->vertex()->point();
        // create std::vector<std::size_t> for each face
        _faces.push_back(std::vector<std::size_t>{(std::size_t)std::distance(_vertices.begin(), std::find(_vertices.begin(), _vertices.end(), p1)),
                                                  (std::size_t)std::distance(_vertices.begin(), std::find(_vertices.begin(), _vertices.end(), p2)),
                                                  (std::size_t)std::distance(_vertices.begin(), std::find(_vertices.begin(), _vertices.end(), p3))});
    }
}

// Copy constructor
polygon::polygon(const polygon& other)
    : _poly(other._poly), _bbox(other._bbox), _solid_bbox(other._solid_bbox),
      mesh(other.mesh), _vertices(other._vertices), _faces(other._faces),
      centroid(other.centroid)
{
    double xafsdf=2343;
    // If you have dynamically allocated resources (e.g., pointers), you should copy them here.
    // For example, if you have a dynamically allocated array:
    // _someData = new SomeType[other.size];
    // Copy the data from 'other' to '_someData' here.
}

// Copy assignment operator
polygon& polygon::operator=(const polygon& other)
{
    if (this != &other) // Self-assignment check
    {
        // Copy all data members from 'other' to 'this' object.
        _poly = other._poly;
        _bbox = other._bbox;
        _solid_bbox = other._solid_bbox;
        mesh = other.mesh;
        _vertices = other._vertices;
        _faces = other._faces;
        centroid = other.centroid;

        // If you have dynamically allocated resources (e.g., pointers), you should copy them here.
        // For example, if you have a dynamically allocated array:
        // delete[] _someData; // Delete the existing data (if any).
        // _someData = new SomeType[other.size];
        // Copy the data from 'other' to '_someData' here.
    }
    return *this;
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
        // volume += CGAL::abs(CGAL::volume(p1, p2, p3, centroid));
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

        // surface_area += CGAL::sqrt(CGAL::squared_area(p1, p2, p3));
    }
    return surface_area;
}


boost::optional<std::pair<int, double>> polygon::intersection(const Eigen::Vector3d &point, const Eigen::Vector3d &step) const
{
    // Convert point and step to ExactKernel types
    Kernel::Point_3 source(point(0), point(1), point(2));
    Kernel::Point_3 target(step(0) + point(0), step(1) + point(1), step(2) + point(2));
    Kernel::Segment_3 segment(source, target);

    // Now use exact computations to find intersections
    std::vector<std::pair<boost::variant<Kernel::Point_3, Kernel::Segment_3>, Primitive>> intersections;
    _AABBtree->all_intersections(segment, std::back_inserter(intersections));

    if (intersections.size() == 0)
    {
        return boost::optional<std::pair<int, double>>();
    }

    std::set<Kernel::Point_3> unique_intersection_points;
    double min_distance = std::numeric_limits<double>::max();  // Start with a large value
    Kernel::Point_3 closest_point;
    int index_closest_face=-1;

    for(const auto& intersection : intersections)
    {
        if(auto *pt = boost::get<Kernel::Point_3>(&intersection.first))
        {
            if (unique_intersection_points.find(*pt) == unique_intersection_points.end())  // Check if the point is new
            {
                double distance = CGAL::sqrt(CGAL::squared_distance(source, *pt));
                unique_intersection_points.insert(*pt);

                if (distance < min_distance)  // Update the closest point if needed
                {
                    min_distance = distance;
                    closest_point = *pt;
                    index_closest_face = (intersection.second.id() - triangle_faces.begin());
                }
            }
        }
    }
    if (min_distance == std::numeric_limits<double>::max() || index_closest_face == -1)
    {
        throw std::runtime_error("Polygon::intersection -> Segment overlaps an edge");
    }
    return std::pair<int,double>(index_closest_face, min_distance);
}

/*
boost::optional<std::pair<int, double>> polygon::intersectionSecond(const Eigen::Vector3d &point, const Eigen::Vector3d &step) const
{
    Kernel::Segment_3 segment(Kernel::Point_3(point(0), point(1), point(2)), Kernel::Point_3(step(0) + point(0), step(1) + point(1), step(2) + point(2)));
    std::vector<std::pair<boost::variant<Kernel::Point_3, Kernel::Segment_3>, Primitive>> intersections;
    _AABBtree->all_intersections(segment, std::back_inserter(intersections));
    auto check_box_intersection_updated = _AABBtree->bbox();
    double min_distance = std::numeric_limits<double>::max();
    int min_index = -1;

    // No intersection found
    if (intersections.size() == 0)
    {
        return boost::optional<std::pair<int, double>>();
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
            //if (min_edge_distance < 1e-18)
            //{
              //throw std::runtime_error("Polygon::intersection -> Intersection found but point is too close to an edge, uncertain");
            //}

            Kernel::FT squared_distance = CGAL::squared_distance(point_intersected, segment.source());
            double distance = CGAL::sqrt(CGAL::to_double(squared_distance));
            if (distance < min_distance)
            {
                auto iter = intersection.second.id();
                min_index = (iter - triangle_faces.begin());
                min_distance = distance;
            } 
        }
        else{
            throw std::runtime_error("Polygon::intersection -> Intersection found but segment is parallel and lies on face");
        }
    }

    if (min_distance == std::numeric_limits<double>::max())
    {
        throw std::runtime_error("Polygon: Intersection error.");
    }
    // TODO: Unify epsilon values in the code
    double t = CGAL::to_double(min_distance) / step.norm();
    if ((1 - t) < 1e-8 || t < 1e-8)
    {
       throw std::runtime_error("Polygon::intersection -> Remaining step is too short, might be uncertain");
    }

    return std::pair<int, double>(min_index, min_distance);
}
*/
bool polygon::containsPoint(const Eigen::Vector3d &point) const
{
    if (_solid_bbox.xmin() > point(0) || _solid_bbox.xmax() < point(0) || _solid_bbox.ymin() > point(1) ||
        _solid_bbox.ymax() < point(1) || _solid_bbox.zmin() > point(2) || _solid_bbox.zmax() < point(2))
    {
        return false;
    }
    // Define a ray that starts at the point and goes in the positive x direction
    CGAL::Ray_3<Kernel> ray(CGAL::Point_3<Kernel>(point(0), point(1), point(2)), CGAL::Vector_3<Kernel>(1, 0, 0));

    // Find all intersections between the ray and the triangles in the AABB tree
    std::vector<std::pair<boost::variant<Kernel::Point_3, Kernel::Segment_3>, Primitive>> intersections;

   _AABBtree->all_intersections(ray, std::back_inserter(intersections));

    std::set<Kernel::Point_3> unique_intersection_points;

    for (const auto &intersection : intersections)
    {
        if (intersection.first.which() == 0)
        {
            Kernel::Point_3 point_intersected = boost::get<Kernel::Point_3>(intersection.first);

            // Rounding to some desired precision (for example, 6 decimal places)
            double x = std::round(point_intersected.x() * 10000) / 10000;
            double y = std::round(point_intersected.y() * 10000) / 10000;
            double z = std::round(point_intersected.z() * 10000) / 10000;

            Kernel::Point_3 rounded_point(x, y, z);

            if (unique_intersection_points.find(rounded_point) == unique_intersection_points.end())
            {
                unique_intersection_points.insert(rounded_point);

                Primitive intersection_primitive = intersection.second;
                auto iter = intersection_primitive.id();
                int index = iter - triangle_faces.begin();

                // Now, process the intersection...
            }
        }
    }
    // Check if the number of intersections is odd
    return (unique_intersection_points.size() % 2) == 1;
}

Eigen::Vector3d polygon::getNormalVector(int index_face) const
{
    // Helper functions
    auto is_counterclockwise_oriented = [](const Kernel::Triangle_3 &triangle)
    {
        if (triangle.is_degenerate())
        {
            throw std::runtime_error("Polygon::getNormalVector -> Triangle is degenerate");
        }

        const Kernel::Point_3 A = triangle.vertex(0);
        const Kernel::Point_3 B = triangle.vertex(1);
        const Kernel::Point_3 C = triangle.vertex(2);

        Kernel::Vector_3 normal = CGAL::cross_product(B - A, C - A);
        Kernel::Point_3 D = A + normal;

        // Check if the triangle_3 is oriented counterclockwise
        CGAL::Orientation orientation = CGAL::orientation(A, B, C, D);

        return orientation == CGAL::POSITIVE;
    };

    auto opposite = [](const Kernel::Triangle_3 &triangle)
    {
        const Kernel::Point_3 A = triangle.vertex(0);
        const Kernel::Point_3 B = triangle.vertex(1);
        const Kernel::Point_3 C = triangle.vertex(2);
        return Kernel::Triangle_3(A, C, B);
    };

    // Obtain triangle from index
    Kernel::Triangle_3 triangle = triangle_faces[index_face];

    auto isoriented = is_counterclockwise_oriented(triangle);
    Kernel::Triangle_3 new_triangle = triangle;
    if (!isoriented)
    {
        new_triangle = opposite(triangle);
        bool is_new_oriented = is_counterclockwise_oriented(new_triangle);
        if (!is_new_oriented)
        {
            throw std::runtime_error("Polygon::getNormalVector -> Triangle is not oriented");
        }
    }

    // Obtain the normal vector from Kernel::Triangle_3
    Kernel::Vector_3 normal_vector = new_triangle.supporting_plane().orthogonal_vector();

    double normal_vector_norm = CGAL::to_double(CGAL::sqrt(CGAL::to_double(normal_vector.squared_length())));
    // Obtain norm of normal_vector in double
    normal_vector = normal_vector / normal_vector_norm;

    // Transform to Eigen
    Eigen::Vector3d normal_vector_eigen(CGAL::to_double(normal_vector.x()), CGAL::to_double(normal_vector.y()), CGAL::to_double(normal_vector.z()));

    return normal_vector_eigen;
}

bool polygon::isPolygonClosed()
{
    // Iterate over halfedges in _poly
    for (auto he = _poly.halfedges_begin(); he != _poly.halfedges_end(); he++)
    {
        if (he->is_border())
        {
            // If the halfedge is a border, it can't have a twin
            return false;
        }
        if (he->opposite() == nullptr)
        {
            // If the halfedge doesn't have a twin, the polyhedron is not closed
            return false;
        }
    }
    // If all halfedges have twins, the polyhedron is closed
    return true;
}

void polygon::createPolygon(const Eigen::MatrixXd &vertices, const Eigen::MatrixXd &faces)
{
    for (int i = 0; i < vertices.rows(); i++)
    {
        _vertices.push_back(CGAL::Point_3<Kernel>((double)vertices(i, 0), (double)vertices(i, 1), (double)vertices(i, 2)));
    }
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
        _faces.push_back(face_input);
    }
    B.end_surface();
}