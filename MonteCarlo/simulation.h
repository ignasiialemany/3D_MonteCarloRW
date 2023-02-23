#ifndef INC_3DRANDOMWALK_SIMULATION_H
#define INC_3DRANDOMWALK_SIMULATION_H

#include "walkers.h"
#include <iostream>
#include <cstdlib>

class simulation {

    simulation() = default;
    explicit simulation(walkers particles);

public:
    bool seedParticlesInBox(Eigen::VectorXd boundingbox);
    //void performScan(sequence sequence_input, substrate substrate_input);


private:
    walkers particles;

    static bool checkBoundingBox(Eigen::VectorXd box)




};


#endif //INC_3DRANDOMWALK_SIMULATION_H
