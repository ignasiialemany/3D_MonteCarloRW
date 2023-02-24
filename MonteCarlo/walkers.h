#ifndef UNTITLED_WALKERS_H
#define UNTITLED_WALKERS_H

#include <set>
#include <random>
#include <Eigen/Dense>

class walkers {

public:
    walkers()=default;
    ~walkers() = default;

    walkers(int N_p, int seed, const std::string &stepType); //Construtor definition

    std::set<unsigned int> get_seed_set();
    void init_positions(Eigen::MatrixXd &pos);
    void set_position(Eigen::VectorXd &pos , int index);
    Eigen::VectorXd get_position(int index) {return positions(index,Eigen::all).transpose();};
    Eigen::MatrixXd get_positions(){return positions;};
    [[nodiscard]] int get_Np() const {return N_p;};
    [[nodiscard]] int get_seed() const {return rng_seed;};

private:
    int N_p; //Number of particles
    int rng_seed; //Seed
    std::string stepType = "constant"; //Step type
    bool initialized = false; //Flag for initialized (might delete)
    std::set<unsigned int> seed_set;
    Eigen::MatrixXd positions; //Positions
    Eigen::MatrixXd phase; //Phase
    Eigen::MatrixXd flag; //Flags


    void generate_unique_seeds();
};


#endif //UNTITLED_WALKERS_H
