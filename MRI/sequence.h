//
// Created by Alemany Juvanteny, Ignasi on 26/02/2023.
//

#ifndef INC_3DRANDOMWALK_SEQUENCE_H
#define INC_3DRANDOMWALK_SEQUENCE_H

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <yaml-cpp/yaml.h>

struct sequence_parameters{
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
    Eigen::VectorXd dt;
    Eigen::MatrixXd gG;
};

class sequence {
public:
    sequence(const std::string& filename) {
        YAML::Node config = YAML::LoadFile(filename);
        parameters.type = config["sequence"]["type"].as<std::string>();
        parameters.G_max = config["sequence"][parameters.type]["G"].as<double>();
        parameters.Delta = config["sequence"][parameters.type]["Delta"].as<double>();
        parameters.epsilon = config["sequence"][parameters.type]["epsilon"].as<double>();
        if (parameters.type=="MCSE"){
            parameters.delta = config["sequence"][parameters.type]["delta1"].as<double>();
            parameters.delta2 = config["sequence"][parameters.type]["delta2"].as<double>();
        }
        else{
            parameters.delta = config["sequence"][parameters.type]["delta"].as<double>();
        }
        parameters.alpha90 = config["sequence"][parameters.type]["alpha90"].as<double>();
        parameters.alphaR0 = config["sequence"][parameters.type]["alphaR0"].as<double>();
        parameters.gamma = config["sequence"][parameters.type]["gamma"].as<double>();
        parameters.number_of_timesteps = config["sequence"][parameters.type]["N_t"].as<int>();
        parameters.dt_max_free = config["sequence"][parameters.type]["dt_max"][0].as<double>();
        parameters.dt_max_grad = config["sequence"][parameters.type]["dt_max"][1].as<double>();
    };
    sequence_parameters parameters;
    void create();
    void discretize(Eigen::VectorXd durations, Eigen::VectorXd ids);

};

#endif //INC_3DRANDOMWALK_SEQUENCE_H
