//
// Created by Alemany Juvanteny, Ignasi on 23/02/2023.
//

#include "simulation.h"

simulation::simulation(walkers particles) {
    this->particles = particles;
}

bool simulation::seedParticlesInBox(Eigen::VectorXd bounding_box) {

    // Init RNG based on particles seed
    std::random_device rd;
    std::mt19937 gen(rd());

    // Number of boxes
    try {
        checkBoundingBox(bounding_box);
    }
    catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }

    //TODO Implement so that if max_point==min_point the code works for 2D,1D boxes
    Eigen::VectorXd max_points = bounding_box({1,3,5});
    Eigen::VectorXd min_points = bounding_box({0,2,4});
    Eigen::VectorXd side_lengths = max_points - min_points;

    Eigen::MatrixXd initial_positions(particles.get_Np(),3);

    for(int i=0;i<initial_positions.rows();i++){
        for(int j=0;j<initial_positions.cols();j++){
            //TODO implement random distribution here
            initial_positions(i,j) = 2;
        }
    }

    particles.set_positions(initial_positions);











    /*
    // Check if and dimension has zero length in specified dimensions (e.g. only x and y)
    hasZeroDim = (sidelengths_raw(Eigen::all, Eigen::seqN(0, dimension)).array() == 0.0).any();
    if (hasZeroDim) {
        printf("simulation::seedParticlesInBox::inconsistent, each dimension must be either zero or non-zero across all boxes");
        flag = false;
    }

    // Removing zero dimensions (?????)
    sidelengths = sidelengths_raw(Eigen::all, Eigen::seqN(0, dimension));

    // Calculating the box volume of all the boxes
    boxVolumes = sidelengths.rowwise().prod();
    total_boxVolume = boxVolumes.sum();

    // Calculating the particles per box
    particlesPerBox_theo = boxVolumes / total_boxVolume * obj.get_N_p();
    particlesPerBox_prel = Eigen::floor(particlesPerBox_theo.array());
    int missing_particles = obj.get_N_p() - particlesPerBox_prel.sum();

    // Refilling missing particles
    particlesPerBox = simulation::refill(particlesPerBox_prel, missing_particles);
}

else if(particlesPerBox_input > 0){
particlesPerBox.
resize(N_b);
particlesPerBox = Eigen::VectorXd::Ones(N_b, 1) * particlesPerBox;
}
*/



return false;
}

bool simulation::checkBoundingBox(Eigen::VectorXd box) {
    int dimensions = (int) box.size() / 2;
    for (int j = 0; j < dimensions; j++) {
        bool isIncosistent = box(2 * j) > box(2 * j + 1);
        if (isIncosistent) {
            throw std::runtime_error("Bounding Box is inconsistent, min > max.");;
        }
    }
}



