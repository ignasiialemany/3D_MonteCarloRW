#include <iostream>
#include "yaml-cpp/yaml.h"
#include "MonteCarlo/simulation.h"
#include <boost/filesystem.hpp>
#include <cmath>

int main()
{
   //We first set the substrate using a std::vector<std::string> of polygon paths

   auto current_path = boost::filesystem::current_path();
   auto parent_path = current_path.parent_path();
   std::string parent_path_str = parent_path.string();

   std::string grandparent_path_str = parent_path_str.substr(0, parent_path_str.size() - 5);
   std::string path = "/Users/ia4118/CLionProjects/3DRandomWalk/poly_geos/rw_mesh.off";
   std::vector<std::string> polygons;
   polygons.push_back(path);

   //We need the transform and the strain to create the substrate.

   transform t(0., false, true);
   t.set_block(81., 81., 101.);
   std::function<double(double)> strain = [](double t) -> double {
      double sin_value = 2*M_PI*t*(1.0/800.0);
      double val = 0.*std::sin(sin_value);
      return val;
   };

   substrate sub(polygons, t, strain);
   
   Eigen::VectorXd voxel(6);
   //double buffer_zone = std::sqrt(6 * 2.5 * seq.dt.sum());
   voxel << 40.5, 40.5, 40.5 , 40.5, 50.5, 50.5;
   sub.setVoxel(voxel);
   sub.setBoundaryType("reflective");


   sequence seq;
   Eigen::VectorXd durations(20000);
   durations.setConstant(0.01);
   //Create a VetorXd same length as durations with all elements equal to 1
   Eigen::VectorXd ids = Eigen::VectorXd::Ones(durations.size());
   seq.dt = durations;
   seq.gG = ids;
   int seed = 12;

   simulation sim(300000, seed);
   sim.seedParticlesInBox(sub);
   auto start = std::chrono::high_resolution_clock::now();
   sim.precomputeSubstrate(sub, seq.dt);
   auto end = std::chrono::high_resolution_clock::now();
   std::chrono::duration<double> elapsed = end - start;
   std::cout << "Precomputed substrate took: " << elapsed.count() << " s\n";

   double kappa = 0.1;
   sim.params.isOutput = false;
   sim.params.kappa = kappa;
   sim.params.D_ecs = 2.5;
   sim.params.D_ics = 1.0;
   sim.params.cores = 1.;

   //std::string file_init_positions = "RW_seed_" + std::to_string(seed) + "_kappa_" + std::to_string(kappa) + "_fixed_step_init.csv";
   //sim.writeParticlesState(file_init_positions,sub);
   
   std::cout << "Starting simulation" << std::endl;
   sim.performScan(sub, seq);   
   std::string file_final_positions = "3d_largeperm_nostrain_RW_seed_" + std::to_string(seed) + "_kappa_"  + std::to_string(kappa) + "_fixed_step_final.csv";
   sim.writeParticlesState(file_final_positions,sub);
   /* 
   printf("path: %s\n", path.c_str());
   std::vector<std::string> polygons;
   //create new polygon in std::vector polygons
   polygons.push_back(path);

   substrate sub(polygons);

   //create lambda function with input double and output Eigen::Vector3d
   std::function<double(double)> strain = [](double t) -> double {
      double sin_value = 2*M_PI*t*(1.0/20.0);
      double val = 0.5*std::sin(sin_value);
      return val;
   };

   Eigen::VectorXd voxel(6);
   //double buffer_zone = std::sqrt(6 * 2.5 * seq.dt.sum());
   voxel << 40., 40., 40., 40., 59., 59.;
   sub.setVoxel(voxel);
   transform t(0., false, true);
   t.set_block(80., 80., 118.);

   sequence seq;

   Eigen::VectorXd durations(1000);
   durations.setConstant(0.01);

   //Create a VetorXd same length as durations with all elements equal to 1
   Eigen::VectorXd ids = Eigen::VectorXd::Ones(durations.size());
   seq.dt = durations;
   seq.gG = ids;
   int seed = 7;

   simulation sim(100000, seed);

   //Set transform that it could be either// shiftedBlock and rotatedBlock
   sub.setTransform(t);
   //set strain
   sub.setStrain(strain);
   // "none", "reflective", "periodic"
   sub.setBoundaryType("none");

   sim.seedParticlesInBox(sub);
   double kappa = 0.0;
   sim.params.isOutput = false;
   sim.params.kappa = kappa;
   sim.params.D_ecs = 0;
   sim.params.D_ics = 1.0;
   sim.params.cores = 8;

   std::string file_init_positions = "RW_seed_" + std::to_string(seed) + "_kappa_" + std::to_string(kappa) + "_fixed_step_init.csv";
   sim.writeParticlesState(file_init_positions,sub);

   sim.performScan(sub, seq);

   std::string file_final_positions = "RW_seed_" + std::to_string(seed) + "_kappa_"  + std::to_string(kappa) + "_fixed_step_final.csv";
   sim.writeParticlesState(file_final_positions,sub);
   */
   return 0;
}