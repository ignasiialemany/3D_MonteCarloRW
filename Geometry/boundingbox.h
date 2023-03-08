//
// Created by Alemany Juvanteny, Ignasi on 27/02/2023.
//

#ifndef INC_3DRANDOMWALK_BOUNDINGBOX_H
#define INC_3DRANDOMWALK_BOUNDINGBOX_H

#include <Eigen/Dense>

class boundingbox
{
public:
    boundingbox() = default;

    boundingbox(Eigen::VectorXd bb_range)
    {
        boundingbox::initialize(bb_range);
    };

    void initialize(Eigen::VectorXd bb_range_input)
    {
        bb_range = bb_range_input;
        DxDyDz << (bb_range(3) - bb_range(0)), (bb_range(4) - bb_range(1)), (bb_range(5) - bb_range(2));
        // Computing bounding box property
        if (DxDyDz(0) > 0 and DxDyDz(1) > 0 and DxDyDz(2) > 0)
        { // Check if any length is negative
            bb_volume = boundingbox::bounding_box_volume(DxDyDz);
            bb_surface_area = boundingbox::bounding_box_surface(DxDyDz);
        }
        else
        {
            printf("geometry::boundingbox::bounding box must be in increasing pairs");
        }
    }

    bool containsPoint(Eigen::Vector3d &point)
    {
        bool isInside = (point(0) > bb_range(0) && point(0) < bb_range(1))     // inside in x, and
                        && (point(1) > bb_range(2) && point(1) < bb_range(3))  // inside in y, and
                        && (point(2) > bb_range(4) && point(2) < bb_range(5)); // inside in
        return isInside;
    }

    // Bounding box range
    Eigen::VectorXd bb_range = Eigen::VectorXd::Zero(6);

    // Bounding box volume
    double bb_volume;
    // Bounding box surface area
    double bb_surface_area;

private:
    // Length of the bounding box
    Eigen::Vector3d DxDyDz;

    // Calculating the bounding box volume
    double bounding_box_volume(Eigen::Vector3d DxDyDz)
    {
        return DxDyDz.prod();
    }

    // Calculating the bounding box surface area
    double bounding_box_surface(Eigen::Vector3d DxDyDz)
    {
        return 2 * (DxDyDz(0) * DxDyDz(1) + DxDyDz(0) * DxDyDz(2) + DxDyDz(1) * DxDyDz(2));
    }
};

#endif // INC_3DRANDOMWALK_BOUNDINGBOX_H
