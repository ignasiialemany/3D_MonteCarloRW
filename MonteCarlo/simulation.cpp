//
// Created by Alemany Juvanteny, Ignasi on 23/02/2023.
//

#include "simulation.h"
#include <omp.h>

simulation::simulation(const walkers &particles_input)
{
    _particles = particles_input;
}

bool simulation::seedParticlesInBox(const Eigen::VectorXd &bounding_box)
{
    // Number of boxes
    try
    {
        checkBoundingBox(bounding_box);
    }
    catch (const std::runtime_error &e)
    {
        return false;
    }

    // TODO: Implement 2D logic in seeding max_point==min_point the code works for 2D,1D boxes
    // Seeding particles in bounding box
    seeding(bounding_box);
    return true;
}

boost::variant<bool, std::runtime_error> simulation::checkBoundingBox(const Eigen::VectorXd &box)
{
    int dimensions = (int)box.size() / 2;
    for (int j = 0; j < dimensions; j++)
    {
        bool isIncosistent = box(2 * j) > box(2 * j + 1);
        if (isIncosistent)
        {
            throw std::runtime_error("Bounding Box is inconsistent, min > max.");
            ;
        }
    }
    return true;
}

void simulation::seeding(const Eigen::VectorXd &box)
{
    Eigen::MatrixXd initial_positions(_particles.get_number_of_particles(), 3);
    // Specifically for the seeding we will use the general seed as well
    std::mt19937 seeding(_particles.get_global_seed());
    std::uniform_real_distribution<> x_interval(box(0), box(1));
    std::uniform_real_distribution<> y_interval(box(2), box(3));
    std::uniform_real_distribution<> z_interval(box(4), box(5));
    for (int i = 0; i < _particles.get_number_of_particles(); i++)
    {
        Particle &particle = _particles.get_particle(i);
        particle.position(0) = x_interval(seeding);
        particle.position(1) = y_interval(seeding);
        particle.position(2) = z_interval(seeding);
    }
}

void simulation::performScan(const substrate &substrate, const sequence &sequence)
{
#pragma omp parallel num_threads(_params.cores) default(none) shared(_particles, _params, substrate, sequence)
    {
#pragma omp for
        for (int index_particle = 0; index_particle < _particles.get_number_of_particles(); index_particle++)
        {
            Particle &particle = _particles.get_particle(index_particle);
            one_walker(particle, substrate, sequence);
        }
    }
#pragma omp barrier
}

void simulation::one_walker(Particle &particle, const substrate &substrate, const sequence &sequence)
{
    // TODO: Fill myo index for each particle, with try/catch for errors.

    std::mt19937 rng_engine(particle.seed);

    for (int i = 0; i < sequence.parameters.dt.size(); i++)
    {
        // Get dT and dG values
        double dt_magnitude = sequence.parameters.dt(i);
        double dG_magnitude = sequence.parameters.gG(i);
        std::pair<double, double> dTdG = std::make_pair(dt_magnitude, dG_magnitude);

        // TODO: Calculate phase from dG and dT and position

        // Initialize counter and flags
        int counter = 0;
        bool step_success = false;
        while (!step_success)
        {
            counter++;
            if (counter > 50)
            {
                particle.flag = 2;
                throw std::runtime_error("Runtime error: Stepping has stopped. Particle flagged");
            }

            try
            {
                // TODO: Probably delete sequence input from one_dt
                one_dt(particle, substrate, sequence, rng_engine, dTdG);
                step_success = true;
            }
            catch (const std::exception &ex)
            {
                // TODO: Identify exceptions type from one_dt
            }
        }
    }
}

void simulation::one_dt(Particle &particle, const substrate &substrate, const sequence &sequence, std::mt19937 &rng_engine, std::pair<double, double> &dTdG)
{
    // Get step
    Eigen::Vector3d step = simulation::getStep(rng_engine, _params.dimension, _params.step_type);
    double dt_magnitude = dTdG.first;
    double dG_magnitude = dTdG.second;
    double D_coeff_old, D_coeff_new;

    // Get D coefficient
    D_coeff_old = (particle.myocyte_index != 1) ? _params.D_ecs : _params.D_ics;
    D_coeff_new = D_coeff_old;
    step = step * std::sqrt(2 * dt_magnitude * D_coeff_old);

    // Get variables for Transit Model (TODO: Maybe transit model as a _params enumerator or class?)
    double probability_of_transit, ds, term, D_low, l_low, p_fieremans, p_maruyama;

    int counter = 0;
    // TODO: Unify convergence epsilon in class, maybe simulation parameter?
    double convergence_eps = 1e-12;
    while (step.norm() > convergence_eps)
    {
        D_coeff_old = D_coeff_new;

        // Throw logic_error if we we try to intersect more than 50 times
        if (counter > 50) throw std::logic_error("Stepping error, stopped after trying 50 times");

        // Get transform and local position/step
        transform_info transform_data = substrate.getGlobalToLocal(particle.position);
        Eigen::Vector3d local_position = transform_data.local_position;
        Eigen::Vector3d local_step = utility_substrate::rotate_y(step, -transform_data.angle);

        // Get intersection data
        boost::variant<bool, std::pair<int, double>> intersection_data = substrate.intersection(local_position, local_step);

        // If intersection is type bool
        if (intersection_data.which() == 0)
        {
            Eigen::Vector3d new_position = local_position + local_step;

            //TODO: Check if particle wants to leave the block, if it does update position till intersection and continue while loop

            //TODO: If not just full update position and it will break the loop
        }
        else
        {
            // If there is intersection then get the distance to membrane (normal to face, so need normal vector)


            
        }

    //TODO: Update position for the substep, make sure it is done last as any throw error would be caught in one_walker
    }

}

Eigen::Vector3d simulation::getStep(std::mt19937 &rng_engine, int dimension, std::string step_type)
{
    Eigen::Vector3d step = Eigen::Vector3d::Zero(3);
    std::unordered_map<std::string, int> string_map = {
        {"constant", 0},
        {"normal", 1}};
    bool isValid = false;
    double max_step = 5;
    std::uniform_int_distribution<int> uniform_dist(0, 1);
    std::normal_distribution<double> normal(0, 1);
    while (!isValid)
    {
        switch (string_map[step_type])
        {
        case 0:
            for (int i = 0; i < dimension; i++)
            {
                // Maps -1 and 1 from 0 , 1
                step(i) = 2 * uniform_dist(rng_engine) - 1;
            }
            break;
        case 1:
            for (int i = 0; i < dimension; i++)
            {
                step(i) = normal(rng_engine);
            }
            break;
        }
        isValid = (step.array() * step.array()).sum() > max_step;
    }
    return step;
}

void simulation::set_parameters(int cores, int dimension, std::string &step_type)
{
    _params.cores = cores;
    _params.dimension = dimension;
    _params.step_type = step_type;
}
