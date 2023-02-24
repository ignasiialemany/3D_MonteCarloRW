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
    flag = Eigen::MatrixXd::Zero(N_p, 3);
}

void walkers::generate_unique_seeds() {
    //Generate unique seeds for each particle using a RNG generator with the global seed
    std::mt19937 gen(rng_seed);
    while(seed_set.size() < N_p){
        unsigned int seed = gen();
        if (seed_set.find(seed) == seed_set.end()){
            seed_set.insert(seed);
        }
    }
}

std::set<unsigned int> walkers::get_seed_set() {
    return seed_set;
}

void walkers::set_positions(Eigen::MatrixXd &pos) {
    this->positions = pos;
}



