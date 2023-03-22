#ifndef UNTITLED_WALKERS_H
#define UNTITLED_WALKERS_H

#include <set>
#include <random>
#include <Eigen/Dense>
#include "../Utility/fileWriter.h"

struct Particle
{
    Eigen::Vector3d position = Eigen::Vector3d::Zero(3);
    Eigen::Vector3d phase = Eigen::Vector3d::Zero(3);
    unsigned int flag = 0;
    unsigned int seed = 0;
    int myocyte_index = -1;
    int index = -1;
    fileWriter file;
};

class walkers
{

public:
    walkers() = default;
    ~walkers() = default;
    walkers(int Np, int seed); // Construtor definition
    Particle &get_particle(int index) { return _particles.at(index); };
    int get_number_of_particles() const { return _number_of_particles; };
    int get_global_seed() const { return _rng_seed; };
    void openFile(const std::string &filename, int particle_index);
    void writeParameters(const std::string &filepath);
    void initializeFiles();

private:
    int _number_of_particles = 0; // Number of particles
    int _rng_seed = 0;            // Seed
    bool _initialized = false;    // Flag for initialized (might delete)
    std::set<unsigned int> _seed_set;
    std::vector<Particle> _particles; // Positions
    void generate_unique_seeds();
    std::string output_path_var;
};

#endif // UNTITLED_WALKERS_H
