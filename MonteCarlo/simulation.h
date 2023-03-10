#ifndef INC_3DRANDOMWALK_SIMULATION_H
#define INC_3DRANDOMWALK_SIMULATION_H

#include <iostream>
#include <variant>
#include <stdexcept>
#include <cstdlib>
#include "../MRI/sequence.h"
#include "walkers.h"
#include "../Substrate/substrate.h"
#include <boost/variant.hpp>

struct parameters
{
    int cores;
    int dimension;
    std::string step_type;
    std::string transit_model;
    double D_ecs;
    double D_ics;
};

class simulation
{

public:
    simulation() = default;
    explicit simulation(const walkers &particles_input);
    // Copy constructor
    simulation(const simulation &other)
    {
        this->_particles = other._particles;
    }

    // returning *this allows for chaining assignment operations
    simulation &operator=(const simulation &other)
    {
        if (this != &other)
        {
            // Copy all data members from the other object
            this->_particles = other._particles;
            this->_params = other._params;
        }
        return *this;
    }
    bool seedParticlesInBox(const Eigen::VectorXd &boundingbox);
    walkers getParticles() { return _particles; };
    void performScan(const substrate &substrate, const sequence &sequence);
    void set_parameters(int cores, int dimension, std::string &step_type);

private:
    walkers _particles;
    parameters _params;
    static boost::variant<bool, std::runtime_error> checkBoundingBox(const Eigen::VectorXd &box);
    void one_walker(Particle &particle, const substrate &substrate, const sequence &sequence);
    void one_dt(Particle &particle, const substrate &substrate, const sequence &sequence, std::mt19937 &rng_engine, std::pair<double, double> &dTdG);
    void seeding(const Eigen::VectorXd &box);
    static Eigen::Vector3d getStep(std::mt19937 &rng_engine, int dimension, std::string step_type);
    double max_step = 5;
};

#endif // INC_3DRANDOMWALK_SIMULATION_H
