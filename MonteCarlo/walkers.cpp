//
// Created by Alemany Juvanteny, Ignasi on 21/02/2023.
//

#include "walkers.h"
#include <string>
#include <cstdlib>
#include <chrono>
#include <ctime>
#include <sstream>
#include <iomanip>
#include <algorithm>

walkers::walkers(const int Np, const int seed)
{
    _number_of_particles = Np;
    _rng_seed = seed;
    _particles.resize(_number_of_particles);
    // Generate unique seeds for each particle and init positions/phase/flags
    boost::filesystem::path curr_path = boost::filesystem::current_path();
    boost::filesystem::path parent_path = curr_path.parent_path();

    // Obtain the grandparent path using substring of parent_path
    std::string grand_parent_path = parent_path.string().substr(0, parent_path.string().size() - 5);

    // Create a directory with current time, if exists, throw an error
    auto now = std::chrono::system_clock::now();
    std::time_t now_c = std::chrono::system_clock::to_time_t(now);

    std::stringstream ss;
    ss << std::put_time(std::localtime(&now_c), "%H.%M");
    std::string directoryName = ss.str();

    output_path_var = grand_parent_path + "output_seed_" + std::to_string(_rng_seed) + "_" + directoryName + "/";
    std::string command = "mkdir " + output_path_var;
    auto command_char = command.c_str();
    auto status = std::system(command_char);
    if (status != 0)
    {
        throw std::runtime_error("Directory exists: " + output_path_var + ". You should delete it first");
    }

    //Generate seeds for each particle
    generate_unique_seeds();
}

void walkers::openFile(const std::string &filename, int particle_index)
{
    //Open file for each particle and write header
    Particle &particle = get_particle(particle_index);
    std::string filepath = output_path_var + filename;
    particle.file = fileWriter(filepath);
    std::string header = "positionX, positionY, positionZ, phaseX, phaseY, phaseZ, myocyte_index, flag \n";
    particle.file.addHeader(header);
}

void walkers::writeParameters(const std::string &filename, substrate &sub)
{
    auto current_path = boost::filesystem::current_path();
    auto parent_path = current_path.parent_path();
    std::string parent_path_str = parent_path.string();
    // Get substring of parent_path_str subtracting the last 5 characters
    std::string grandparent_path_str = parent_path_str.substr(0, parent_path_str.size() - 5);
    std::string filepath = grandparent_path_str + "seed_" + std::to_string(_rng_seed) + "_" + filename;
    // iterate through particles and write their positions
    fileWriter file(filepath);
    std::string header = "positionX, positionY, positionZ, phaseX, phaseY, phaseZ, myocyte_index, flag \n";
    file.addHeader(header);
    for (int i = 0; i < _number_of_particles; i++)
    {
        Particle &particle = get_particle(i);
        auto transform = sub.getLocalFromGlobal(particle.position);
        auto local_pos = transform.local_position;
        std::string text = std::to_string(local_pos(0)) + "," + std::to_string(local_pos(1)) + "," + std::to_string(local_pos(2)) + "," + std::to_string(particle.phase(0)) + "," + std::to_string(particle.phase(1)) + "," + std::to_string(particle.phase(2)) + "," + std::to_string(particle.myocyte_index) + "," + std::to_string(particle.flag) + "\n";
        file.write(text);
    }
}

void walkers::initializeFiles()
{
    for(int i=0; i<_number_of_particles;i++){
        std::string filename = "/particle_" + std::to_string(i) + ".csv";
        openFile(filename, i);
    }
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
            _seed_set.insert(seed);
            counter++;
        }
    }
}
