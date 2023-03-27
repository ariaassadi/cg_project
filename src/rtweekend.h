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

inline glm::vec3 reflect(const glm::vec3 &v, const glm::vec3 &n)
{
    return v - 2 * glm::dot(v, n) * n;
}

inline glm::vec3 mult(const glm::vec3 &e, const double t)
{
    return glm::vec3(e[0] * t, e[1] * t, e[2] * t);
}

inline glm::vec3 refract(const glm::vec3 &uv, const glm::vec3 &n, double etai_over_etat)
{
    auto cos_theta = fmin(glm::dot(-uv, n), 1.0);
    glm::vec3 r_out_perp = mult(uv + mult(n, cos_theta), etai_over_etat);
    glm::vec3 r_out_parallel = mult(n, -sqrt(fabs(1.0 - glm::dot(r_out_perp, r_out_perp))));
    return r_out_perp + r_out_parallel;
}

#endif RTWEEKEND_H