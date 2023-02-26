//
// Created by Alemany Juvanteny, Ignasi on 26/02/2023.
//

#include "sequence.h"

void sequence::create() {
    if (parameters.type == "PGSE" or parameters.type == "STEAM"){
        Eigen::VectorXd durations(9), ids(9);
        durations << parameters.alpha90, parameters.epsilon, parameters.delta, parameters.epsilon,
                parameters.Delta-(2*parameters.epsilon + parameters.delta),
                parameters.epsilon, parameters.delta, parameters.epsilon, parameters.alphaR0;
        ids << 0, 1, 2, 3, 0, -1, -2, -3, 0;
        discretize(durations, ids);
    }
    else if (parameters.type == "MCSE" or parameters.type == "M2SE"){
        Eigen::VectorXd durations(15), ids(15);
        double del1 = parameters.delta + 2*parameters.epsilon;
        double del2 = parameters.delta2 + 2*parameters.epsilon;
        parameters.Delta = (del2 * (-2*del1+parameters.epsilon) + del1*parameters.epsilon)/(del1-del2);
        durations << parameters.alpha90, parameters.epsilon, parameters.delta, parameters.epsilon, parameters.epsilon,
                parameters.delta2, parameters.epsilon, parameters.Delta-(del1+del2), parameters.epsilon,
                parameters.delta2, parameters.epsilon, parameters.epsilon, parameters.delta, parameters.epsilon, parameters.alphaR0;
        ids << 0, 1, 2, 3, -1, -2, -3, 0, 1, 2, 3, -1, -2, -3, 0;
        discretize(durations, ids);
    }
    else{
        parameters.gamma = 1;
        parameters.dt = Eigen::VectorXd::Ones(parameters.number_of_timesteps)*parameters.dt_max_free;
        parameters.gG = Eigen::MatrixXd::Zero(parameters.number_of_timesteps, 1);
    }
    parameters.gG = parameters.gG*parameters.gamma;
}

void sequence::discretize(Eigen::VectorXd durations, Eigen::VectorXd ids) {
// %   ids := designation of interval
// %       0 - flat (free, gradient OFF)
// %       % gradients: +/- sign of id represents sign of Gmax at the end of the gradient
// %       1 - (+/-) gradient ramp-up
// %       2 - (+/-) gradient flat
// %       3 - (+/-) gradient ramp-down
    Eigen::VectorXd Nt_intervals(ids.rows());

    // calculate the target time step dt
    double dt_aim, dt_free, dt_grad;
    dt_aim = durations.sum()/parameters.number_of_timesteps;
    dt_free = std::fmin(dt_aim, parameters.dt_max_free);
    dt_grad = std::fmin(dt_aim, parameters.dt_max_grad);

    for (int i = 0; i < ids.rows(); i++){
        if (ids(i) == 0){
            Nt_intervals(i) = std::ceil(durations(i)/dt_free);
        }
        else{
            Nt_intervals(i) = std::ceil(durations(i)/dt_grad);
        }
    }

    parameters.dt = Eigen::VectorXd::Zero(0);
    parameters.gG = Eigen::MatrixXd::Zero(0,0);

    // Temporarily store data
    double gA, gB;
    for (int i = 0; i < durations.rows(); i++){
        double Nt_i = Nt_intervals(i);
        double dt_i = durations(i)/Nt_i;
        // For repetition
        Eigen::VectorXd dt_temp = parameters.dt;
        Eigen::VectorXd dt_rep = Eigen::VectorXd::Ones(Nt_i)*dt_i;
        parameters.dt.resize(dt_temp.rows() + Nt_i);
        parameters.dt << dt_temp, dt_rep;

        // gradient
        double id_i = ids(i);
        switch (std::abs((int)id_i)){ // use absolute value for switch and use sign(id_i) inside cases
            case 0: // flat, gradient off
                gA = 0;
                gB = gA;
                break;
            case 1: // ramp-up gradient
                gA = 0;
                gB = std::signbit(id_i)*parameters.G_max;
                break;
            case 2: // flat, gradient on
                gA = std::signbit(id_i)*parameters.G_max;
                gB = gA;
                break;
            case 3: // ramp-down gradient
                gA = std::signbit(id_i)*parameters.G_max;
                gB = 0;
                break;
            default:
                printf("sequence::Ids input out of bound, please check the values of Ids");
        }
        Eigen::MatrixXd gG_temp = parameters.gG;
        Eigen::MatrixXd gvals = Eigen::VectorXd::LinSpaced(Nt_i, gA, gB - (gB-gA)/(Nt_i));
        // gvals.resize(Nt_i, 1); // Removing the last element
        parameters.gG.resize(gG_temp.rows() + (Nt_i),1);
        if (i ==0){
            parameters.gG << gvals;
        }
        else{
            parameters.gG << gG_temp, gvals;
        }
    }
}



