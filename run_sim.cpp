#include <iostream>
#include "yaml-cpp/yaml.h"
#include "MonteCarlo/simulation.h"
#include <boost/filesystem.hpp>

int main()
{

   auto current_path = boost::filesystem::current_path();
   auto parent_path = current_path.parent_path();
   std::string parent_path_str = parent_path.string();
   // Get substring of parent_path_str subtracting the last 5 characters
   std::string grandparent_path_str = parent_path_str.substr(0, parent_path_str.size() - 5);
   std::string file_path = "2d_geo_rw.mat";
   std::string seq_path = "sequence.yaml";
   std::string full_path = grandparent_path_str + file_path;
   std::string seq_full_path = grandparent_path_str + seq_path;

   std::cout << full_path << std::endl;

   // Set transform and substrate
   std::cout << "Setting substrate" << std::endl;
   auto data = utility_substrate::read_mat_file(full_path);
   substrate substrate(data.first, data.second);
   // Create sequence
   sequence seq;

   Eigen::VectorXd voxel(6);
   //double buffer_zone = std::sqrt(6 * 2.5 * seq.dt.sum());
   voxel << 47. , 52., 42., 47., 50., 52.;
   substrate.setVoxel(voxel);
   transform t(0., false, false);
   t.set_block(109., 76., 126.6);
   substrate.setTransform(t);

   // Initialize simulation
   std::cout << "Initializing simulation" << std::endl;
   std::vector<double> kappas{0.05};

   Eigen::VectorXd durations(10000);
   durations.setConstant(501./10000.);

   //Create a VetorXd same length as durations with all elements equal to 1
   Eigen::VectorXd ids = Eigen::VectorXd::Ones(durations.size());
   seq.dt = durations;
   seq.gG = ids;


   for (auto &kappa : kappas)
   {
      // Create simulation
      int seed = rand() % 100000 + 1;  // Generate a random number between 1 and 100000
      simulation sim(500000, seed);
      // Seed particles
      sim.seedParticlesInBox(substrate);

      // Write initial positions
      std::string file_init_positions = "RW_seed_" + std::to_string(seed) + "_kappa_" + std::to_string(kappa) + "_fixed_step_init.csv";
      sim.writeParticlesState(file_init_positions,substrate);

      // Set parameters
      sim.params.cores = 14;
      sim.params.isOutput = false;
      sim.params.kappa = kappa;
      // Perform scan

      // start chrono
      auto start = std::chrono::high_resolution_clock::now();
      sim.performScan(substrate, seq);
      // end chrono
      auto finish = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> elapsed = finish - start;
      std::cout << "Elapsed simulation time: " << elapsed.count() << " seconds" << std::endl;

      // Write final positions
      std::string file_final_positions = "RW_seed_" + std::to_string(seed) + "_kappa_"  + std::to_string(kappa) + "_fixed_step_final.csv";
      sim.writeParticlesState(file_final_positions,substrate);
   }

   return 0;
}