#include <iostream>
#include "yaml-cpp/yaml.h"
#include "MonteCarlo/simulation.h"
#include <boost/filesystem.hpp>

int main()
{

   auto current_path = boost::filesystem::current_path();
   auto parent_path = current_path.parent_path();
   std::string parent_path_str = parent_path.string();
   //Get substring of parent_path_str subtracting the last 5 characters
   std::string grandparent_path_str = parent_path_str.substr(0, parent_path_str.size() - 5);
   std::string file_path = "geometry_1.mat";
   std::string seq_path = "sequence.yaml";
   std::string full_path = grandparent_path_str + file_path;
   std::string seq_full_path = grandparent_path_str + seq_path;

   std::cout << full_path << std::endl;

   //Set transform and substrate
   std::cout << "Setting substrate" << std::endl;
   auto data = utility_substrate::read_mat_file(full_path);
   substrate substrate(data.first, data.second);
   Eigen::VectorXd voxel(6);
   voxel << 0.0, 2800., 0.0, 2800., 0.0, 8000.;
   substrate.setVoxel(voxel);
   transform t(0.01, true, false);
   t.set_block(495.4, 392.3, 126.6);
   substrate.setTransform(t);

   //Create sequence
   sequence seq(seq_full_path);
   seq.create();

   std::cout << "Initializing simulation" << std::endl;
   //Create simulation
   simulation sim(5000,1234);
   sim.seedParticlesInBox(substrate);

   std::string file_init_positions = "init.csv";
   sim.writeParticlesState(file_init_positions);
   std::string step_type = "constant";
   std::string transit_model = "whatever";
   sim.set_parameters(8, 3, step_type, transit_model, 2.5, 1., 0.05);
   sim.performScan(substrate, seq);
   std::string file_final_positions = "final.csv";
   sim.writeParticlesState(file_final_positions);
   return 0;
}