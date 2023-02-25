//
// Created by Alemany Juvanteny, Ignasi on 23/02/2023.
//

#include "simulation.h"
#include <omp.h>

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
    //Eigen::VectorXd max_points = bounding_box({1, 3, 5});
    //Eigen::VectorXd min_points = bounding_box({0, 2, 4});
    //Eigen::VectorXd side_lengths = max_points - min_points;

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
    std::mt19937 seeding(_particles.get_global_seed());
    std::uniform_real_distribution<> x_interval(box(0), box(1));
    std::uniform_real_distribution<> y_interval(box(2), box(3));
    std::uniform_real_distribution<> z_interval(box(4), box(5));

    for (int i = 0; i < initial_positions.rows(); i++) {
        initial_positions(i, 0) = x_interval(seeding);
        initial_positions(i, 1) = y_interval(seeding);
        initial_positions(i, 2) = z_interval(seeding);
    }

    _particles.init_positions(initial_positions);
}

void simulation::performScan(double diffusion_time, double time_step, int cores) {
    int num_threads = cores; // Set the number of threads to the number of cores available
    int number_of_particles= _particles.get_Np();
    #pragma omp parallel num_threads(num_threads) default(none) shared(_particles, number_of_particles,diffusion_time,time_step)
    {
        Eigen::Vector3d position_i;
        Eigen::Vector3d phase_i;
        int flag_i;
        int seed_i;
        #pragma omp for private(position_i)
        for (int i = 0; i < number_of_particles; i++) {
            position_i = _particles.get_position(i);
            phase_i = _particles.get_phase(i);
            flag_i = _particles.get_flag(i);
            seed_i = _particles.get_seed(i);
            //TODO: Walker should have sequence gradient and substrate geometry as well and not diffusion/time_step
            one_walker(i,position_i,phase_i,flag_i,diffusion_time,time_step, seed_i);
            this->_particles.set_position(position_i,i);
        }
    }
    #pragma omp barrier
}

void simulation::one_walker(int index_particle, Eigen::Vector3d &position, Eigen::Vector3d &phase, int &flag, double diffusion_time, double time_step, int seed) {

    //TODO: This will have to be replaced with sequence gradient/timing
    Eigen::VectorXd time = Eigen::VectorXd::LinSpaced((int)diffusion_time/time_step,0,diffusion_time);
    std::mt19937 rng_engine(seed);

    for(int i=0; i<time.size(); i++){

        int counter = 0;
        bool step_success = false;
        //While the particle has not stepped
        while(!step_success){
            try{

                counter++;
                if (counter > 50){
                    std::cout << "Stepping has stopped. Particle flagged" << std::endl;
                    flag=2;
                    break;
                }

                //Eigen::VectorXd position_raw = one_dt(position, dt, substrate_input, myoindex, i_particle, i);
                //position = position_raw({0,1,2});
                //myoindex = position_raw(3);






            }
        }
    }
}

void
simulation::one_dt(int index_particle, Eigen::Vector3d &position, Eigen::Vector3d &phase, int &flag, double time_step,
                   std::mt19937 &rng_engine) {

    Eigen::Vector3d step;
    double D_coeff_old;
    double D_coeff_new;

    //TODO Implement function to get step based on stepType

    //TODO Retrieve myoindex and choose Diffusivity coefficient to D_coeff_old

    //TODO init loop till the step has finished. Do as many substeps as needed








}






