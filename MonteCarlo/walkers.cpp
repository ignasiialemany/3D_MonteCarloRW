//
// Created by Alemany Juvanteny, Ignasi on 21/02/2023.
//

#include "walkers.h"

walkers::walkers(const int N_p, const int seed, const std::string &stepType) {
    this->N_p = N_p;
    this->rng_seed = seed;
    this->stepType = stepType;

    //Generate unique seeds for each particle and init positions/phase/flags
    generate_unique_seeds();
    positions = Eigen::MatrixXd::Zero(N_p, 3);
    phase = Eigen::MatrixXd::Zero(N_p, 3);
    flag = Eigen::VectorXd::Zero(N_p);
}

void walkers::generate_unique_seeds() {
    //Generate unique seeds for each particle using a RNG generator with the global seed
    std::mt19937 gen(rng_seed);
    int counter = 0;
    seeds = Eigen::VectorXd::Zero(N_p);
    while(seed_set.size() < N_p){
        unsigned int seed = gen();
        if (seed_set.find(seed) == seed_set.end()){
            seeds(counter) = seed;
            seed_set.insert(seed);
            counter = counter + 1;
        }
    }
}

void walkers::init_positions(Eigen::MatrixXd &pos) {
    this->positions = pos;
}

void walkers::set_position(Eigen::Vector3d& pos, int index) {
    this->positions(index,Eigen::all) = pos.transpose();
}



