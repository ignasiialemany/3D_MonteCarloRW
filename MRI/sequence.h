//
// Created by Alemany Juvanteny, Ignasi on 26/02/2023.
//

#ifndef INC_3DRANDOMWALK_SEQUENCE_H
#define INC_3DRANDOMWALK_SEQUENCE_H

#include <iostream>
#include <vector>
#include <yaml-cpp/yaml.h>
#include <Eigen/Dense>

struct sequence_parameters
{
    std::string type;
    double G_max;
    double Delta;
    double epsilon;
    double delta;
    double delta2;
    double alpha90;
    double alphaR0;
    double gamma;
    int number_of_timesteps;
    double dt_max_free;
    double dt_max_grad;
};

class sequence
{
public:
    sequence()=default;
    sequence(const std::string &filename)
    {
        YAML::Node config = YAML::LoadFile(filename);
        std::string sequence_type = config["sequence"]["type"].as<std::string>();
        parameters.type = sequence_type;
        YAML::Node params = config["sequence"][sequence_type];
        parameters.G_max = params["Gmax"].as<double>();
        parameters.epsilon = params["epsilon"].as<double>();
        if (parameters.type == "MCSE")
        {
            parameters.delta = params["delta1"].as<double>();
            parameters.delta2 = params["delta2"].as<double>();
        }
        else
        {
            parameters.Delta =params["Delta"].as<double>();
            parameters.delta = params["delta"].as<double>();
        }
        parameters.alpha90 = params["alpha90"].as<double>();
        parameters.alphaR0 = params["alphaRO"].as<double>();
        parameters.gamma = config["sequence"]["gamma"].as<double>();
        parameters.number_of_timesteps = config["sequence"]["N_t"].as<int>();
        parameters.dt_max_free = config["sequence"]["dt_max"][0].as<double>();
        parameters.dt_max_grad = config["sequence"]["dt_max"][1].as<double>();
    };
    sequence_parameters parameters;
    void create();
    void discretize(Eigen::VectorXd durations, Eigen::VectorXd ids);
    Eigen::VectorXd dt;
    Eigen::VectorXd gG;
};

#endif // INC_3DRANDOMWALK_SEQUENCE_H
