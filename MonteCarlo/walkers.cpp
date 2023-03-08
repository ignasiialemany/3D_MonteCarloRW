//
// Created by Alemany Juvanteny, Ignasi on 21/02/2023.
//

#include "walkers.h"

walkers::walkers(const int Np, const int seed)
{
    _number_of_particles = Np;
    _rng_seed = seed;
    // Generate unique seeds for each particle and init positions/phase/flags
    generate_unique_seeds();
}

void walkers::generate_unique_seeds()
{
    // Generate unique seeds for each particle using a RNG generator with the global seed
    std::mt19937 gen(_rng_seed);
    int counter = 0;
    while (_seed_set.size() < _number_of_particles)
    {
        unsigned int seed = gen();
        if (_seed_set.find(seed) == _seed_set.end())
        {
            Particle &particle = get_particle(counter);
            particle.seed = seed;
            counter++;
        }
    }
}
