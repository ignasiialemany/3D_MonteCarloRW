//
// Created by Alemany Juvanteny, Ignasi on 23/02/2023.
//

#include "simulation.h"
#include <omp.h>

simulation::simulation(const walkers &particles_input) {
    _particles = particles_input;
}

bool simulation::seedParticlesInBox(const Eigen::VectorXd &bounding_box) {
    // Number of boxes
    try {
        checkBoundingBox(bounding_box);
    }
    catch (const std::runtime_error &e) {
        return false;
    }

    //TODO: Implement 2D logic in seeding max_point==min_point the code works for 2D,1D boxes
    //Eigen::VectorXd max_points = bounding_box({1, 3, 5});
    //Eigen::VectorXd min_points = bounding_box({0, 2, 4});
    //Eigen::VectorXd side_lengths = max_points - min_points;

    //Seeding particles in bounding box
    seeding(bounding_box);
    return true;
}

std::variant<bool, std::runtime_error> simulation::checkBoundingBox(const Eigen::VectorXd &box) {
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
    Eigen::MatrixXd initial_positions(_particles.get_number_of_particles(), 3);
    //Specifically for the seeding we will use the general seed as well
    std::mt19937 seeding(_particles.get_global_seed());
    std::uniform_real_distribution<> x_interval(box(0), box(1));
    std::uniform_real_distribution<> y_interval(box(2), box(3));
    std::uniform_real_distribution<> z_interval(box(4), box(5));
    for (int i = 0; i < _particles.get_number_of_particles(); i++) {
        Particle& particle = _particles.get_particle(i);
        particle.position(0) = x_interval(seeding);
        particle.position(1) = y_interval(seeding);
        particle.position(2) = z_interval(seeding);
    }
}

void simulation::performScan(const substrate &substrate, const sequence &sequence) {
#pragma omp parallel num_threads(_params.cores) default(none) shared(_particles, _params, substrate, sequence)
    {
#pragma omp for
        for (int index_particle = 0; index_particle < _particles.get_number_of_particles(); index_particle++) {
            Particle& particle = _particles.get_particle(index_particle);
            //TODO: Walker should have sequence gradient and substrate geometry as well and not diffusion/time_step
            one_walker(particle,substrate, sequence);
        }
    }
#pragma omp barrier
}

void simulation::one_walker(Particle &particle, const substrate &substrate, const sequence &sequence) {

    //TODO: This will have to be replaced with sequence gradient/timing
    std::mt19937 rng_engine(particle.seed);

    for (int i = 0; i < sequence.parameters.dt.size(); i++) {
        int counter = 0;
        bool step_success = false;
        //While the particle has not stepped
        while (!step_success) {
            try {
                counter++;
                if (counter > 50) {
                    std::cout << "Stepping has stopped. Particle flagged" << std::endl;
                    particle.flag = 2;
                    break;
                }
                //Update position
                double dt_magnitude = sequence.parameters.dt(i);
                double dG_magnitude = sequence.parameters.gG(i);
                std::pair<double,double> step_magnitude = std::make_pair(dt_magnitude, dG_magnitude);
                one_dt(particle, substrate, sequence, rng_engine, step_magnitude);
                //If any throw error is called by one_dt will be caught here
                step_success = true;
            }
            catch (const std::exception &ex) {

            }


        }
    }
}

void simulation::one_dt( Particle &particle, const substrate &substrate, const sequence &sequence, std::mt19937 &rng_engine, std::pair<double,double> &step_magnitude) {
    Eigen::Vector3d step = getStep(rng_engine,_params.dimension,_params.step_type);
    double dt_magnitude = step_magnitude.first;
    double dG_magnitude = step_magnitude.second;
    double D_coeff_old;

    if (particle.myocyte_index!= 1){
        D_coeff_old = _params.D_ecs;
    }
    else{
        D_coeff_old = _params.D_ics;
    }
    double D_coeff_new = D_coeff_old;
    step = step * std::sqrt(2 * dt_magnitude * D_coeff_old);

    double probability_of_transit, ds, term; // Distrance normal to boundary
    double D_low, l_low, p_fieremans, p_maruyama;
    int counter = 0;
    try {
        double convergence_eps = 1e-12;
        while (step.norm() > convergence_eps) {
            D_coeff_old = D_coeff_new;

            if (counter > 50) { throw std::logic_error("Stepping error, stopped after trying 50 times"); };

            //TODO: Transform position from global 2 local
            //TODO: rotate step
            //TODO: check intersection

            if (!intersectInfo.empty) {
                next_position = position + step;


            }


        }


    }
    catch (const std::exception &ex) {


    }

    //TODO Retrieve myoindex and choose Diffusivity coefficient to D_coeff_old

    //TODO init loop till the step has finished. Do as many substeps as needed

}

Eigen::Vector3d simulation::getStep(std::mt19937 &rng_engine, int dimension, std::string step_type) {
    Eigen::Vector3d step = Eigen::Vector3d::Zero(3);
    std::unordered_map<std::string, int> string_map = {
            {"constant", 0},
            {"normal",   1}
    };
    bool isValid = false;
    std::uniform_int_distribution<int> uniform_dist(0, 1);
    std::normal_distribution<double> normal(0, 1);
    while (!isValid) {
        switch (string_map[step_type]) {
            case 0:
                for (int i = 0; i < dimension; i++) {
                    //Maps -1 and 1 from 0 , 1
                    step(i) = 2 * uniform_dist(rng_engine) - 1;
                }
                break;
            case 1:
                for (int i = 0; i < dimension; i++) {
                    step(i) = normal(rng_engine);
                }
                break;
        }
        isValid = (step.array() * step.array()).sum() > max_step;
    }
    return step;
}

void simulation::set_parameters(int cores, int dimension, std::string &step_type) {
    _params.cores = cores;
    _params.dimension = dimension;
    _params.step_type = step_type;
}








