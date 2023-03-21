#ifndef RTWEEKEND_H
#define RTWEEKEND_H

#include <cmath>
#include <limits>
#include <memory>
#include <cstdlib>

const double infinity = std::numeric_limits<double>::infinity();

inline double random_double()
{
    // Returns a random real in [0,1).
    return rand() / (RAND_MAX + 1.0);
}

inline double random_double(double min, double max)
{
    // Returns a random real in [min,max).
    return min + (max - min) * random_double();
}

// inline glm::vec3 random_in_unit_sphere()
// {
//     glm::vec3 p;
//     do
//     {
//         p = 2.0f * glm::vec3(random_double(), random_double(), random_double()) - glm::vec3(1, 1, 1);
//     } while (glm::dot(p, p) >= 1.0);
//     return p;
// }

inline static glm::vec3 random_vector(double min, double max)
{
    return glm::vec3(random_double(min, max), random_double(min, max), random_double(min, max));
}

inline glm::vec3 random_in_unit_sphere()
{
    while (true)
    {
        auto p = random_vector(-1, 1);
        if ((p[0] * p[0] + p[1] * p[1] + p[2] * p[2]) >= 1)
            continue;
        return p;
    }
}

inline double clamp(double x, double min, double max)
{
    if (x < min)
        return min;
    if (x > max)
        return max;
    return x;
}

#endif RTWEEKEND_H