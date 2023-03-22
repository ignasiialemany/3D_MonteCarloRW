//
// Created by Alemany Juvanteny, Ignasi on 23/02/2023.
//

#include "simulation.h"
#include <typeinfo>

bool simulation::seedParticlesInBox(const substrate &substrate)
{
    Eigen::VectorXd voxel_boundingbox = substrate.getVoxel();
    // Number of boxes
    try
    {
        checkBoundingBox(voxel_boundingbox);
    }
    catch (const std::runtime_error &e)
    {
        return false;
    }

    // TODO: Implement 2D logic in seeding max_point==min_point the code works for 2D,1D boxes
    // Seeding particles in bounding box
    seeding(voxel_boundingbox);
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
    Eigen::MatrixXd initial_positions(_particles->get_number_of_particles(), 3);
    // Specifically for the seeding we will use the general seed as well
    std::mt19937 seeding(_particles->get_global_seed());
    std::uniform_real_distribution<> x_interval(box(0), box(1));
    std::uniform_real_distribution<> y_interval(box(2), box(3));
    std::uniform_real_distribution<> z_interval(box(4), box(5));
    for (int i = 0; i < _particles->get_number_of_particles(); i++)
    {
        Particle &particle = _particles->get_particle(i);
        particle.position(0) = x_interval(seeding);
        particle.position(1) = y_interval(seeding);
        particle.position(2) = z_interval(seeding);
    }
}

void simulation::writeToFile(Particle &particle)
{
    // Add particle.position, particle.myocyte_index, particle.seed, particle.index
    std::string text = std::to_string(particle.position(0)) + "," + std::to_string(particle.position(1)) + "," + std::to_string(particle.position(2)) + "," + std::to_string(particle.phase(0)) + "," + std::to_string(particle.phase(1)) + "," + std::to_string(particle.phase(2)) + "," + std::to_string(particle.myocyte_index) + "," + std::to_string(particle.flag) + "\n";
    particle.file.write(text);
}

void simulation::performScan(const substrate &substrate, const sequence &sequence)
{
    if (params.isOutput)
    {
        _particles->initializeFiles();
    }

    int completed_particles = 0;
#pragma omp parallel num_threads(params.cores) default(none) shared(_particles, params, substrate, sequence, std::cout, completed_particles)
    {
#pragma omp for
        for (int index_particle = 0; index_particle < _particles->get_number_of_particles(); index_particle++)
        {
            Particle &particle = _particles->get_particle(index_particle);
            particle.index = index_particle;
            try
            {
                one_walker(particle, substrate, sequence);
            }
            catch (const std::exception &ex)
            {
#pragma omp critical
                std::cout << "Fatal error : " << index_particle << ": " << ex.what() << std::endl;
            }
#pragma omp critical
            completed_particles++;
            if (index_particle % 300 == 0)
            {
                double percentage = 100 * (double)completed_particles / (double)_particles->get_number_of_particles();
                std::cout << "Completed " << percentage << "% of particles" << std::endl;
            }
        }
    }
#pragma omp barrier
    std::cout << "Simulation finished" << std::endl;
}

void simulation::one_walker(Particle &particle, const substrate &substrate, const sequence &sequence)
{
    particle.myocyte_index = substrate.searchPolygon(particle.position, "global");
    std::mt19937 rng_engine(particle.seed);

    for (int i = 0; i < sequence.dt.size(); i++)
    {
        // Get dT and dG values
        double dt_magnitude = sequence.dt(i);
        double dG_magnitude = sequence.gG(i);
        // TODO: Calculate phase from dG and dT and position
        particle.phase = particle.phase + (dG_magnitude * dt_magnitude) * particle.position;
        // Initialize counter and flags
        int counter = 0;
        bool step_success = false;

        while (!step_success)
        {
            counter++;
            if (counter > 50)
            {
                particle.flag = 2;
                // std::cout << "Particle index COUNTER REACHED" << particle.index << " " << std::endl;
                throw std::runtime_error("Runtime error: Stepping has stopped. Particle flagged");
            }

            try
            {

                // TODO: Probably delete sequence input from one_dt
                one_dt(particle, substrate, rng_engine, dt_magnitude);
                if (params.isOutput)
                {
                    writeToFile(particle);
                }
                step_success = true;
            }
            catch (const std::exception &ex)
            {
                if (typeid(ex) == typeid(std::logic_error))
                {
                    // TODO: Implement more logic errors in functions
                    particle.flag = 2;
                    std::cout << "Particle Index" << particle.index << " " << ex.what() << std::endl;
                    std::cout << "Time step :" << i << std::endl;
                    return;
                }
                else
                {
                    std::cout << "Particle Index" << particle.index << " " << ex.what() << std::endl;
                    std::string errorMessage = ex.what();
                    if (errorMessage.find("Transform") != std::string::npos)
                    {
                        // The particle local position lays very close to the block (epsilon value), just move a bit
                        particle.position = particle.position + Eigen::Vector3d::Constant(1e-6);
                        // Runtime error is in intersection, just repeat one_dt with another step
                        continue;
                    }
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

template <typename URNG>
void simulation::one_dt(Particle &particle, const substrate &substrate, URNG &rng_engine, double dt)
{
    // Get step
    Eigen::Vector3d step = simulation::getStep(rng_engine, params.dimension, params.step_type);
    // Eigen::Vector3d step = Eigen::Vector3d::Zero(3);
    // step(0) = 1;
    // step(1) = 0;
    // step(2) = 1;

    double dt_magnitude = dt;
    double D_coeff_old, D_coeff_new;

    // Get D coefficient
    D_coeff_old = (particle.myocyte_index != -1) ? params.D_ics : params.D_ecs;
    D_coeff_new = D_coeff_old;
    step = step * std::sqrt(2 * dt_magnitude * D_coeff_old);

    // Get variables for Transit Model (TODO: Maybe transit model as a params enumerator or class?)
    double probability_of_transit, term, D_low, l_low, p_fieremans, p_maruyama;

    int counter = 0;
    // TODO: Unify convergence epsilon in class, maybe simulation parameter?
    double convergence_eps = 1e-12;

    // We get local data
    transform_info transform_data = substrate.getLocalFromGlobal(particle.position);
    Eigen::Vector3d local_position = transform_data.local_position;
    Eigen::Vector3d local_step = utility_substrate::rotate_y(step, -transform_data.angle);

    while (local_step.norm() > convergence_eps)
    {
        Eigen::Vector3d local_normalized_step = local_step / local_step.norm();
        Eigen::Vector3d remaining_step;
        // std::cout << step.norm() << std::endl;

        D_coeff_old = D_coeff_new;

        //Check that the particle is in cardiomycocyte
        if (particle.myocyte_index!=-1){
            if (substrate.searchPolygon(local_position)==-1){
               D_coeff_old = params.D_ecs;
            }
            D_coeff_old = params.D_ics;
        }
        else{
            if (substrate.searchPolygon(local_position)!=-1){
                D_coeff_old = params.D_ics;
            }
            D_coeff_old = params.D_ecs;
        }

        counter++;
        // std::cout << "Counter: " << counter << std::endl;
        //  Throw logic_error if we we try to intersect more than 50 times
        if (counter > 50)
        {
            throw std::logic_error("Stepping error, stopped while loop in one_dt after trying 50 times");
        }
        // Get transform and local position/step

        // auto check_index = substrate.searchPolygon(local_position);
        // if (particle.myocyte_index != check_index)
        //{
        //   std::cout << "Particle index :" << particle.index << " - Particle myocyte index does not match the index of the polygon it is in" << std::endl;
        // throw std::logic_error("Particle myocyte index does not match the index of the polygon it is in");
        //}
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

            if (D_coeff_old == params.D_ecs)
            {
                D_low = params.D_ics;
                l_low = distance_to_intersection_projected * std::sqrt(D_low / D_coeff_old);
                term = 2 * l_low * params.kappa / D_low;
                p_fieremans = term / (term + 1);
                p_maruyama = std::sqrt(params.D_ics / D_coeff_old) > 1 ? 1 : std::sqrt(params.D_ics / D_coeff_old); // equivalent to matlab min(1,sqrt(substrate.D_i/D_old))
                probability_of_transit = p_fieremans * p_maruyama;
            }
            else
            {
                term = 2 * distance_to_intersection_projected * params.kappa / params.D_ics;
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
                    D_coeff_new = params.D_ics;
                }
                else
                {
                    particle.myocyte_index = -1;
                    D_coeff_new = params.D_ecs;
                }
                // Update remaining step to account for the change in D
                remaining_step = remaining_step * std::sqrt(D_coeff_new / D_coeff_old);
            }
            else
            {

                // Compute the reflection based on the remaining normalized step
                Eigen::Vector3d remaining_step_normalized = remaining_step / remaining_step.norm();
                Eigen::Vector3d reflected_step = remaining_step_normalized - 2 * normal * normal.dot(remaining_step_normalized);
                // Reflect the remaining_step
                remaining_step = reflected_step * remaining_step.norm();
            }

            local_position = local_position + step_to_intersection;
            local_position = local_position + remaining_step * 1e-6;
            remaining_step = remaining_step * (1 - 1e-6);
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
                local_position = local_position + remaining_step * 1e-6;
                remaining_step = remaining_step * (1 - 1e-6);

                // Get global position position now and update local_position
                Eigen::Vector3d global_position_after_block = substrate.getGlobalFromLocal(local_position, transform_data.iX, transform_data.iY, transform_data.iZ);
                Eigen::Vector3d global_step = utility_substrate::rotate_y(remaining_step, transform_data.angle);

                // Recalculate transform data and local data (note that local_step is also rotated back might be that block falls into new angle)
                transform_data = substrate.getLocalFromGlobal(global_position_after_block);

                //Can we do any type of check? Maybe check the highest diff between the two local positions?
                

                local_position = transform_data.local_position;
                remaining_step = utility_substrate::rotate_y(global_step, -transform_data.angle);

                int possible_myo = substrate.searchPolygon(local_position);

                if (possible_myo != -1)
                {
                    throw std::runtime_error("Particle has crossed block and directly to cardiomyocyte");
                }
            }
            else
            {
                local_position = local_position + local_step;
                remaining_step = Eigen::Vector3d::Zero(3);
            }
        }

        /*
        if (particle.myocyte_index != -1)
        {
            // Check if particle is in myocyte
            bool check_myo = substrate.containsPoint(particle.myocyte_index, local_position);
            if (!check_myo && substrate.searchPolygon(local_position) == -1)
            {
                //std::cout << "Local position: " << local_position.transpose() << std::endl;
                throw std::logic_error("FINAL ONE_DT: Particle is not in the myocyte");
            }
        }
        else
        {
            int index_myo = substrate.searchPolygon(local_position);
            if (index_myo != -1)
            {
                //std::cout << "Local position: " << local_position.transpose() << " index myo: " << index_myo << std::endl;
                throw std::logic_error("FINAL ONE_DT: Particle should be in ECS but is in myocyte");
            }
        }*/
        local_step = remaining_step;
    }

    // Update particle position
    Eigen::Vector3d next_global_position = substrate.getGlobalFromLocal(local_position, transform_data.iX, transform_data.iY, transform_data.iZ);
    particle.position = next_global_position;
}

template void simulation::one_dt<oneGenerator>(Particle &particle, const substrate &substrate, oneGenerator &rng_engine, double dt);
template void simulation::one_dt<std::mt19937>(Particle &particle, const substrate &substrate, std::mt19937 &rng_engine, double dt);
template Eigen::Vector3d simulation::getStep<oneGenerator>(oneGenerator &rng_engine, int dimension, std::string step_type);
template Eigen::Vector3d simulation::getStep<std::mt19937>(std::mt19937 &rng_engine, int dimension, std::string step_type);
