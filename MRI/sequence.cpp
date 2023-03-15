//
// Created by Alemany Juvanteny, Ignasi on 26/02/2023.
//

#include "sequence.h"

void sequence::create()
{
    if (parameters.type == "PGSE" or parameters.type == "STEAM")
    {
        Eigen::VectorXd durations(9), ids(9);
        durations << parameters.alpha90, parameters.epsilon, parameters.delta, parameters.epsilon,
            parameters.Delta - (2 * parameters.epsilon + parameters.delta),
            parameters.epsilon, parameters.delta, parameters.epsilon, parameters.alphaR0;
        ids << 0, 1, 2, 3, 0, -1, -2, -3, 0;
        discretize(durations, ids);
    }
    else if (parameters.type == "MCSE" or parameters.type == "M2SE")
    {
        Eigen::VectorXd durations(15), ids(15);
        double del1 = parameters.delta + 2 * parameters.epsilon;
        double del2 = parameters.delta2 + 2 * parameters.epsilon;
        parameters.Delta = (del2 * (-2 * del1 + parameters.epsilon) + del1 * parameters.epsilon) / (del1 - del2);
        durations << parameters.alpha90, parameters.epsilon, parameters.delta, parameters.epsilon, parameters.epsilon,
            parameters.delta2, parameters.epsilon, parameters.Delta - (del1 + del2), parameters.epsilon,
            parameters.delta2, parameters.epsilon, parameters.epsilon, parameters.delta, parameters.epsilon, parameters.alphaR0;
        ids << 0, 1, 2, 3, -1, -2, -3, 0, 1, 2, 3, -1, -2, -3, 0;
        discretize(durations, ids);
    }
    else
    {
        //TODO: Check this is correct
        //parameters.gamma = 1;
        dt = Eigen::VectorXd::Ones(parameters.number_of_timesteps) * parameters.dt_max_free;
        gG = Eigen::VectorXd::Zero(parameters.number_of_timesteps);
    }
    gG = gG * parameters.gamma;
}

void sequence::discretize(Eigen::VectorXd durations, Eigen::VectorXd ids)
{
    // %   ids := designation of interval
    // %       0 - flat (free, gradient OFF)
    // %       % gradients: +/- sign of id represents sign of Gmax at the end of the gradient
    // %       1 - (+/-) gradient ramp-up
    // %       2 - (+/-) gradient flat
    // %       3 - (+/-) gradient ramp-down
    Eigen::VectorXd Nt_intervals(ids.rows());

    // calculate the target time step dt
    double dt_aim, dt_free, dt_grad;
    dt_aim = durations.sum() / parameters.number_of_timesteps;
    dt_free = std::fmin(dt_aim, parameters.dt_max_free);
    dt_grad = std::fmin(dt_aim, parameters.dt_max_grad);

    for (int i = 0; i < ids.rows(); i++)
    {
        if (ids(i) == 0)
        {
            Nt_intervals(i) = std::ceil(durations(i) / dt_free);
        }
        else
        {
            Nt_intervals(i) = std::ceil(durations(i) / dt_grad);
        }
    }

    dt = Eigen::VectorXd::Zero(0);
    gG = Eigen::VectorXd::Zero(0);

    // Temporarily store data
    double gA, gB;
    for (int i = 0; i < durations.rows(); i++)
    {
        double Nt_i = Nt_intervals(i);
        double dt_i = durations(i) / Nt_i;
        // For repetition
        Eigen::VectorXd dt_temp = dt;
        Eigen::VectorXd dt_rep = Eigen::VectorXd::Ones(Nt_i) * dt_i;
        dt.resize(dt_temp.rows() + Nt_i);
        dt << dt_temp, dt_rep;

        // gradient
        double id_i = ids(i);
        double sign = std::signbit((int)id_i) ? -1 : 1;
        switch (std::abs((int)id_i))
        {       // use absolute value for switch and use sign(id_i) inside cases
        case 0: // flat, gradient off
            gA = 0;
            gB = gA;
            break;
        case 1: // ramp-up gradient
            gA = 0;
            gB = sign * parameters.G_max;
            break;
        case 2: // flat, gradient on
            gA = sign * parameters.G_max;
            gB = gA;
            break;
        case 3: // ramp-down gradient
            gA = sign * parameters.G_max;
            gB = 0;
            break;
        default:
            printf("sequence::Ids input out of bound, please check the values of Ids");
        }
        Eigen::VectorXd gvals = Eigen::VectorXd::LinSpaced(Nt_i+1, gA, gB);
        // gvals.resize(Nt_i, 1); // Removing the last element
        Eigen::VectorXd gvals_segmented = gvals.segment(0,gvals.rows()-1);
        gG.conservativeResize(gG.rows() + gvals_segmented.rows());
        gG.tail(gvals_segmented.rows()) = gvals_segmented;
    }
}
