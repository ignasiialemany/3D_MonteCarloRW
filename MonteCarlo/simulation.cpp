//
// Created by Alemany Juvanteny, Ignasi on 23/02/2023.
//

#include "simulation.h"

simulation::simulation(const walkers& particles_input) {
    this->_particles = particles_input;
}

bool simulation::seedParticlesInBox(const Eigen::VectorXd &bounding_box) {

    // Number of boxes
    try {
        checkBoundingBox(bounding_box);
    }
    catch (const std::runtime_error &e) {
        return false;
    }

    //TODO Implement so that if max_point==min_point, this will ensure the code works for 2D,1D boxes
    Eigen::VectorXd max_points = bounding_box({1, 3, 5});
    Eigen::VectorXd min_points = bounding_box({0, 2, 4});
    Eigen::VectorXd side_lengths = max_points - min_points;

    //Seeding particles in bounding box
    seeding(bounding_box);

    return true;
}

std::variant<bool,std::runtime_error> simulation::checkBoundingBox(const Eigen::VectorXd &box) {
    int dimensions = (int) box.size() / 2;
    for (int j = 0; j < dimensions; j++) {
        bool isIncosistent = box(2 * j) > box(2 * j + 1);
        if (isIncosistent) {
            throw std::runtime_error("Bounding Box is inconsistent, min > max.");;
        }
    }
    return true;
}

void simulation::seeding(const Eigen::VectorXd &box) {

    Eigen::MatrixXd initial_positions(_particles.get_Np(), 3);
    //Specifically for the seeding we will use the general seed as well
    std::mt19937 seeding(_particles.get_seed());
    std::uniform_real_distribution<> x_interval(box(0), box(1));
    std::uniform_real_distribution<> y_interval(box(2), box(3));
    std::uniform_real_distribution<> z_interval(box(4), box(5));

    for (int i = 0; i < initial_positions.rows(); i++) {
        initial_positions(i, 0) = x_interval(seeding);
        initial_positions(i, 1) = y_interval(seeding);
        initial_positions(i, 2) = z_interval(seeding);
    }

    _particles.set_positions(initial_positions);
}



