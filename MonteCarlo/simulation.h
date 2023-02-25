#ifndef INC_3DRANDOMWALK_SIMULATION_H
#define INC_3DRANDOMWALK_SIMULATION_H

#include "walkers.h"
#include <iostream>
#include <variant>
#include <stdexcept>
#include <cstdlib>

class simulation {

public:

    simulation() = default;
    explicit simulation(const walkers& particles_input);
    // Copy constructor
    simulation(const simulation& other) {
        this->_particles = other._particles;
    }

    //returning *this allows for chaining assignment operations
    simulation& operator=(const simulation& other) {
        if (this != &other) {
            // Copy all data members from the other object
            this->_particles = other._particles;
        }
        return *this;
    }


    bool seedParticlesInBox(const Eigen::VectorXd &boundingbox);
    Eigen::MatrixXd getPositions(){return _particles.get_positions();};
    void performScan(double diffusion_time, double time_step, int cores);


private:
    walkers _particles;
    static std::variant<bool,std::runtime_error> checkBoundingBox(const Eigen::VectorXd &box);
    static void one_walker(int index_particle, Eigen::Vector3d &position, Eigen::Vector3d &phase, int &flag, double diffusion_time, double time_step, int seed);
    static void one_dt(int index_particle, Eigen::Vector3d &position, Eigen::Vector3d &phase, int &flag, double time_step, std::mt19937 &rng_engine);
    void seeding(const Eigen::VectorXd &box);

};


#endif //INC_3DRANDOMWALK_SIMULATION_H
