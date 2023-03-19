#include <yaml-cpp/yaml.h>
#include <catch2/catch_all.hpp>
#include "../MRI/sequence.h"
#include <boost/filesystem.hpp>

TEST_CASE("Load sequence parameters from YAML file", "[sequence]")
{
    auto cwd = boost::filesystem::current_path();
    auto parent_path = cwd.parent_path();
    std::string parent_path_str = parent_path.string();
    std::string grandparent_path_str = parent_path_str.substr(0, parent_path_str.size() - 5);

    std::string file_path = "testing/data/sequence_testing.yaml";
    std::string full_path = grandparent_path_str + file_path;
    std::cout << full_path << std::endl;
    sequence seq(full_path);
    seq.create();
    REQUIRE(seq.parameters.type == "MCSE");
    //Checked comparing with matlab, just throw one case here 
    REQUIRE(std::abs(seq.gG(2)-0.00528018) < 1e-4);
}
