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

    Eigen::VectorXd get_seeds(){return seeds;};
    void init_positions(Eigen::MatrixXd &pos);
    void set_position(Eigen::Vector3d &pos , int index);
    void set_flag(int flag, int index){flag(index)=flag;};
    //TODO might be delete flags
    double get_flag(int index) {return flag(index);};
    int get_seed(int index){return seeds(index);};
    Eigen::Vector3d get_position(int index) {return positions(index,Eigen::all).transpose();};
    Eigen::Vector3d get_phase(int index) {return phase(index,Eigen::all).transpose();};
    Eigen::MatrixXd get_positions(){return positions;};
    Eigen::MatrixXd get_phases(){return phase;};
    std::string get_step_type(){return stepType;};
    [[nodiscard]] int get_Np() const {return N_p;};
    [[nodiscard]] int get_global_seed() const {return rng_seed;};

private:
    int N_p; //Number of particles
    int rng_seed; //Seed
    std::string stepType = "constant"; //Step type
    bool initialized = false; //Flag for initialized (might delete)
    std::set<unsigned int> seed_set;
    Eigen::MatrixXd positions; //Positions
    Eigen::MatrixXd phase; //Phase
    //TODO might be delete flags
    Eigen::VectorXd flag; //Flags
    //TODO might be useful to use a Vector of std RNG engines instead of just int seeds?
    Eigen::VectorXd seeds;


    void generate_unique_seeds();
};


#endif //UNTITLED_WALKERS_H
