//
// Created by Alemany Juvanteny, Ignasi on 26/02/2023.
//

#include "transform.h"

transform_info transform::global2local(Eigen::Vector3d &global_position)
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
        output.angle = (deg_rot_per_m_in_Y * M_PI / 180.) * (y_slice); //This angle is in rad now

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

        double ddxx = (output.iX-1)*dx;
        double ddzz = (output.iZ-1)*dz;

        output.local_position << position_rotated(0) - ddxx, position_y_sliced, position_rotated(2) - ddzz;

        return output;
    }
}

Eigen::Vector3d transform::local2global(Eigen::Vector3d &local_position, int iX, int iY, int iZ)
{
    int iX_new = iX -1;
    int iY_new = iY -1;
    int iZ_new = iZ -1;

    //Translation offset
    double shift = 0;
    if(shift_block){
        shift = utility_substrate::mod(iZ,2)/2;
    }

    Eigen::Vector3d offset = Eigen::Vector3d(iX_new*dx+shift*dx,iY_new*dy,iZ_new*dz);

    Eigen::Vector3d position_global = local_position + offset;

    //TODO: Add sinusoidal displacement to position_rotated(2) 

    //Invert rotation 
    return utility_substrate::rotate_y(position_global, (deg_rot_per_m_in_Y * M_PI / 180.)*position_global(1));
}
