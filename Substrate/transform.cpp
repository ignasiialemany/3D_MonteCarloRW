//
// Created by Alemany Juvanteny, Ignasi on 26/02/2023.
//

#include "transform.h"
#include <unordered_set>

transform_info transform::global2local(const Eigen::Vector3d &global_position) const
{
    transform_info output;

    if (isIdentity)
    {
        output.angle = 0;
        output.local_position = global_position;
        return output;
    }
    else
    {
        double y_slice = utility_substrate::find_yslice(global_position(1), y_slice_minmax);
        output.angle = (deg_rot_per_m_in_Y * M_PI / 180.) * (y_slice); // This angle is in rad now

        Eigen::Vector3d position_rotated = utility_substrate::rotate_y(global_position, -output.angle);

        double position_y_sliced = utility_substrate::mod(position_rotated(1), dy);

        output.iY = 1 + std::floor(position_rotated(1) / dy);
        // TODO: Add sinusoidal displacement to position_rotated(2) and then calculate iZ
        output.iZ = 1 + std::floor(position_rotated(2) / dz);
        if (shift_block)
        {
            // If the block is in even Z rows, shift the X coordinate by half a dx
            position_rotated(0) = position_rotated(0) - utility_substrate::mod(output.iZ, 2) * dx / 2;
        }
        output.iX = 1 + std::floor(position_rotated(0) / dx);

        double ddxx = (output.iX - 1) * dx;
        double ddzz = (output.iZ - 1) * dz;

        output.local_position << position_rotated(0) - ddxx, position_y_sliced, position_rotated(2) - ddzz;

        // throw error if local_position is not within block
        if (output.local_position(0) < 0 || output.local_position(0) > dx || output.local_position(1) < 0 || output.local_position(1) > dy || output.local_position(2) < 0 || output.local_position(2) > dz)
        {
            throw std::runtime_error("Transform::global2local -> Local position is not within block");
        }

        return output;
    }
}

Eigen::Vector3d transform::local2global(const Eigen::Vector3d &local_position, int iX, int iY, int iZ) const
{
    int iX_new = iX - 1;
    int iY_new = iY - 1;
    int iZ_new = iZ - 1;

    // Translation offset
    double shift = 0;
    if (shift_block)
    {
        shift = utility_substrate::mod(iZ, 2) / 2;
    }

    Eigen::Vector3d offset = Eigen::Vector3d(iX_new * dx + shift * dx, iY_new * dy, iZ_new * dz);

    Eigen::Vector3d position_global = local_position + offset;

    double y_slice = dy * iY_new;
    double angle = (deg_rot_per_m_in_Y * M_PI / 180.) * (y_slice); // This angle is in rad now

    // TODO: Add sinusoidal displacement to position_rotated(2)

    // Invert rotation with +angle
    return utility_substrate::rotate_y(position_global, angle);
}

void transform::set_block(double dx, double dy, double dz, double minY, double maxY)
{
    this->dx = dx;
    this->dy = dy;
    this->dz = dz;
    Kernel::Point_3 min_point = Kernel::Point_3(0, 0, 0);
    Kernel::Point_3 max_point = Kernel::Point_3(dx, dy, dz);

    _solid_block = CGAL::Bbox_3(0, 0, 0, dx, dy, dz);
    create_block(min_point, max_point);
    // Create y_slice_minmax
    // create a set of unique y-values
    std::unordered_set<double> y_minvals_set;
    for (double y = 0.0; y <= maxY; y += dy)
    {
        y_minvals_set.insert(y);
    }
    for (double y = -dy; y >= minY; y -= dy)
    {
        y_minvals_set.insert(y);
    }

    Eigen::VectorXd y_minvals(y_minvals_set.size());
    std::copy(y_minvals_set.begin(), y_minvals_set.end(), y_minvals.data());
    std::sort(y_minvals.data(), y_minvals.data() + y_minvals.size());

    this->y_slice_minmax.resize(2, y_minvals.size() - 1);
    this->y_slice_minmax.row(0) = y_minvals.segment(0, y_minvals.size() - 1);
    this->y_slice_minmax.row(1) = y_minvals.tail(y_minvals.size() - 1);
}

void transform::create_block(Kernel::Point_3 min_point, Kernel::Point_3 max_point)
{
    CGAL::Polyhedron_incremental_builder_3<HalfedgeDS> builder(_block.hds(), true);
    builder.begin_surface(8, 12);

    // Add vertices
    std::vector<Kernel::Point_3> vertices = {
        min_point,
        Kernel::Point_3(max_point.x(), min_point.y(), min_point.z()),
        Kernel::Point_3(max_point.x(), max_point.y(), min_point.z()),
        Kernel::Point_3(min_point.x(), max_point.y(), min_point.z()),
        Kernel::Point_3(min_point.x(), min_point.y(), max_point.z()),
        Kernel::Point_3(max_point.x(), min_point.y(), max_point.z()),
        max_point,
        Kernel::Point_3(min_point.x(), max_point.y(), max_point.z())};

    for (const auto &vertex : vertices)
    {
        builder.add_vertex(vertex);
    }

    // Add faces (as triangles)
    std::vector<std::vector<int>> triangle_face_indices = {
        {0, 1, 2}, // Face 0
        {0, 2, 3}, // Face 1
        {4, 6, 5}, // Face 2
        {4, 7, 6}, // Face 3
        {0, 4, 5}, // Face 4
        {0, 5, 1}, // Face 5
        {1, 5, 6}, // Face 6
        {1, 6, 2}, // Face 7
        {2, 6, 7}, // Face 8
        {2, 7, 3}, // Face 9
        {3, 7, 4}, // Face 10
        {3, 4, 0}  // Face 11
    };
    for (const auto &face : triangle_face_indices)
    {
        builder.begin_facet();
        for (int vertex_index : face)
        {
            builder.add_vertex_to_facet(vertex_index);
        }
        Kernel::Triangle_3 face_triangle = Kernel::Triangle_3(vertices[face[0]], vertices[face[1]], vertices[face[2]]);
        triangles.push_back(face_triangle);
        builder.end_facet();
    }
    builder.end_surface();
    _AABBtree = std::make_shared<Tree_AABB>(triangles.begin(), triangles.end());
}
