#ifndef INC_3DRANDOMWALK_SIMULATION_H
#define INC_3DRANDOMWALK_SIMULATION_H

#include "walkers.h"
#include <iostream>
#include <variant>
#include <stdexcept>
#include <cstdlib>

class simulation {

public:
    simulation() = default;
    explicit simulation(const walkers& particles_input);

    bool seedParticlesInBox(const Eigen::VectorXd &boundingbox);
    Eigen::MatrixXd getPositions(){return _particles.get_positions();};
    //void performScan(sequence sequence_input, substrate substrate_input);


private:
    walkers _particles;
    std::variant<bool,std::runtime_error> checkBoundingBox(const Eigen::VectorXd &box);
    void seeding(const Eigen::VectorXd &box);

};


#endif //INC_3DRANDOMWALK_SIMULATION_H
