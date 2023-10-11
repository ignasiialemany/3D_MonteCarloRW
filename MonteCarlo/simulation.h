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
#include <omp.h>
#include "../testing/oneGenerator.h"
#include <memory>

struct parameters
{
    int cores=1;
    int dimension = 3;
    std::string step_type = "constant";
    std::string transit_model;
    double D_ecs = 2.5;
    double D_ics = 1.;
    double kappa = 0.05;
    bool isOutput = false;
};

class simulation
{

public:

    simulation():_particles(nullptr){};

    simulation(int N_p, int seed) : _particles(std::make_shared<walkers>(N_p,seed)){};

    simulation(const simulation &other) : _particles(other._particles){};
    // returning *this allows for chaining assignment operations
    //simulation &operator=(const simulation &other)
    //{
        //if (this != &other)
        //{
            // Copy all data members from the other object
            //this->_particles = other._particles;
            //this->params = other.params;
        //}
        //return *this;
    //}

    bool seedParticlesInBox(const substrate &substrate);
    void performScan(substrate &substrate, const sequence &sequence);
    void precomputeSubstrate(substrate &substrate, Eigen::VectorXd sequence_dt);
    void writeParticlesState(const std::string &file_path, substrate &sub){_particles->writeParameters(file_path,sub);};
    Particle& get_particle(int i){return _particles->get_particle(i);};
    template <typename URNG>
    void one_dt(Particle &particle, const substrate &substrate, URNG &rng_engine, double dt, int index_sequence, double total_time);

    void one_walker(Particle &particle, const substrate &substrate, const sequence &sequence);
    parameters params;

private:
    std::shared_ptr<walkers> _particles;
    static boost::variant<bool, std::runtime_error> checkBoundingBox(const Eigen::VectorXd &box);
    void seeding(const Eigen::VectorXd &box);
    void writeToFile(Particle &particle);
    template <typename URNG>
    static Eigen::Vector3d getStep(URNG &rng_engine, int dimension, std::string step_type);
    double max_step = 5;
};

#endif // INC_3DRANDOMWALK_SIMULATION_H
