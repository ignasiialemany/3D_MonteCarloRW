//
// Created by Alemany Juvanteny, Ignasi on 23/02/2023.
//

#include "simulation.h"
#include <typeinfo>

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
    particle.myocyte_index = substrate.searchPolygon(particle.position);
    std::mt19937 rng_engine(particle.seed);

    for (int i = 0; i < sequence.dt.size(); i++)
    {
        // Get dT and dG values
        double dt_magnitude = sequence.dt(i);
        double dG_magnitude = sequence.gG(i);
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
                one_dt(particle, substrate, rng_engine, dt_magnitude);
                step_success = true;
            }
            catch (const std::exception &ex)
            {
                // If it is a std::logic_error break the loop it means it has tried one_dt for 50 times
                if (typeid(ex) == typeid(std::logic_error))
                {
                    particle.flag = 2;
                    std::cout << ex.what() << std::endl;
                    break;
                }
                else
                {
                    std::cout << ex.what() << std::endl;
                    continue;
                }
            }
        }
    }
}

template <typename URNG>
Eigen::Vector3d simulation::getStep(URNG &rng_engine, int dimension, std::string step_type)
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
            { // Maps -1 and 1 from 0 , 1
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
        isValid = (step.array() * step.array()).sum() <= max_step * max_step;
    }
    return step;
}

void simulation::set_parameters(int cores, int dimension, std::string &step_type, std::string &transit_model, double D_ecs, double D_ics, double kappa)
{
    _params.cores = cores;
    _params.dimension = dimension;
    _params.step_type = step_type;
    _params.transit_model = transit_model;
    _params.D_ecs = D_ecs;
    _params.D_ics = D_ics;
    _params.kappa = kappa;
}

template <typename URNG>
void simulation::one_dt(Particle &particle, const substrate &substrate, URNG &rng_engine, double dt)
{
    // Get step
    Eigen::Vector3d step = simulation::getStep(rng_engine, _params.dimension, _params.step_type);
    double dt_magnitude = dt;
    double D_coeff_old, D_coeff_new;

    // Get D coefficient
    D_coeff_old = (particle.myocyte_index != -1) ? _params.D_ics : _params.D_ecs;
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
        if (counter > 50)
            throw std::logic_error("Stepping error, stopped while loop in one_dt after trying 50 times");

        // Get transform and local position/step
        transform_info transform_data = substrate.getLocalFromGlobal(particle.position);
        Eigen::Vector3d local_position = transform_data.local_position;
        Eigen::Vector3d local_step = utility_substrate::rotate_y(step, -transform_data.angle);
        Eigen::Vector3d local_normalized_step = local_step / local_step.norm();
        Eigen::Vector3d remaining_step;

        // Get intersection data
        boost::optional<std::tuple<int, double, Eigen::Vector3d>> intersection_data = substrate.intersectPolygon(local_position, local_step);

        // If intersection is type bool and false
        if (intersection_data)
        {
            // Get intersection data
            int polygon_index = std::get<0>(*intersection_data);
            double distance_to_intersection = std::get<1>(*intersection_data);
            Eigen::Vector3d normal = std::get<2>(*intersection_data);

            // Obtain the remaining step and the step to the intersection
            double remaining_value = local_step.norm() - distance_to_intersection;
            remaining_step = local_normalized_step * remaining_value;
            Eigen::Vector3d step_to_intersection = local_normalized_step * distance_to_intersection;

            // dot product between normal and normalized step
            double dot_product = normal.dot(local_normalized_step);
            double cos_angle = dot_product > 0 ? dot_product : -dot_product;
            double distance_to_intersection_projected = distance_to_intersection * cos_angle;

            // TODO: Implement transit models, for now we consider the hybrid model

            if (D_coeff_old == _params.D_ecs)
            {
                D_low = _params.D_ics;
                l_low = distance_to_intersection_projected * std::sqrt(D_low / D_coeff_old);
                term = 2 * l_low * _params.kappa / D_low;
                p_fieremans = term / (term + 1);
                p_maruyama = std::sqrt(_params.D_ics / D_coeff_old) > 1 ? 1 : std::sqrt(_params.D_ics / D_coeff_old); // equivalent to matlab min(1,sqrt(substrate.D_i/D_old))
                probability_of_transit = p_fieremans * p_maruyama;
            }
            else
            {
                term = 2 * distance_to_intersection_projected * _params.kappa / _params.D_ics;
                p_fieremans = term / (term + 1);
                // p_maruyama = min(1,sqrt(D_new/D_old)); (in matlab)
                probability_of_transit = p_fieremans;
            }

            std::uniform_real_distribution<> transit(0, 1);
            double U;
            if (typeid(rng_engine) == typeid(oneGenerator))
            {
                U = 0.5; // Hardcoded value for U when rng_engine is of type oneGenerator
            }
            else
            {
                U = transit(rng_engine);
            }

            // Particle is crossing
            if (U < probability_of_transit)
            {
                // If the particle is in ECS update to myocyte index otherwise update to -1
                if (particle.myocyte_index == -1)
                {
                    particle.myocyte_index = polygon_index;
                    D_coeff_new = _params.D_ics;
                }
                else
                {
                    particle.myocyte_index = -1;
                    D_coeff_new = _params.D_ecs;
                }
                // Update remaining step to account for the change in D
                remaining_step = remaining_step * std::sqrt(D_coeff_new / D_coeff_old);
            }
            else
            {

                // Compute the reflection based on the remaining normalized step
                Eigen::Vector3d reflected_step = remaining_step / remaining_step.norm() - 2 * normal * normal.dot(remaining_step / remaining_step.norm());
                // Reflect the remaining_step
                remaining_step = reflected_step * remaining_step.norm();
            }

            local_position = local_position + step_to_intersection;
            local_position = local_position + remaining_step * 1e-8;
            remaining_step = remaining_step * (1 - 1e-8);
        }
        else
        {
            boost::optional<double> intersection_block = substrate.intersectionBlock(local_position, local_step);

            // If it intersects with block
            if (intersection_block)
            {
                double distance = (*intersection_block);
                remaining_step = local_normalized_step * (local_step.norm() - distance);
                Eigen::Vector3d step_to_intersection = local_normalized_step * distance;
                // Update position to intersection + epsilon
                local_position = local_position + step_to_intersection;
                local_position = local_position + remaining_step * 1e-8;
                remaining_step = remaining_step * (1 - 1e-8);
            }
            else
            {
                local_position = local_position + local_step;
                remaining_step = Eigen::Vector3d::Zero(3);
            }
        }

        Eigen::Vector3d next_global_position = substrate.getGlobalFromLocal(local_position, transform_data.iX, transform_data.iY, transform_data.iZ);
        Eigen::Vector3d global_remaining_step = utility_substrate::rotate_y(remaining_step, transform_data.angle);
        particle.position = next_global_position;
        step = global_remaining_step;
    }
}

template void simulation::one_dt<oneGenerator>(Particle &particle, const substrate &substrate, oneGenerator &rng_engine, double dt);
template void simulation::one_dt<std::mt19937>(Particle &particle, const substrate &substrate, std::mt19937 &rng_engine, double dt);
template Eigen::Vector3d simulation::getStep<oneGenerator>(oneGenerator &rng_engine, int dimension, std::string step_type);
template Eigen::Vector3d simulation::getStep<std::mt19937>(std::mt19937 &rng_engine, int dimension, std::string step_type);
