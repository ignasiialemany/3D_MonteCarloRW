#ifndef INC_3D_RANDOMWALK_ONE_GENERATOR_H
#define INC_3D_RANDOMWALK_ONE_GENERATOR_H

#include <random>

class oneGenerator : public std::mt19937
    {
    public:
        using result_type = std::mt19937::result_type;

        oneGenerator() : std::mt19937(default_seed) {}

        result_type operator()()
        {
            return max(); // You can change this to 0 if you always want it to return 0
        }

        static constexpr result_type min()
        {
            return std::mt19937::min();
        }

        static constexpr result_type max()
        {
            return std::mt19937::max();
        }

        const result_type default_seed = 12345;
    };

#endif