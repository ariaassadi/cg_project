#pragma once

#include "rtweekend.h"

namespace rt
{

    class material
    {
    public:
        virtual bool scatter(const Ray &r_in, const HitRecord &rec, glm::vec3 &attenuation, Ray &scattered) const = 0;
    };

    class Lambertian : public material
    {
    public:
        Lambertian(const glm::vec3 &a) : albedo(a) {}
        virtual bool scatter(const Ray &r_in, const HitRecord &rec, glm::vec3 &attenuation, Ray &scattered) const override
        {
            glm::vec3 target = rec.p + rec.normal + random_in_unit_sphere();
            scattered = Ray(rec.p, target - rec.p);
            attenuation = albedo;
            return true;
        }
        glm::vec3 albedo;
    };

    class Metal : public material
    {
    public:
        Metal(const glm::vec3 &a, float f) : albedo(a) { f < 1 ? fuzz = f : fuzz = 1; }
        virtual bool scatter(const Ray &r_in, const HitRecord &rec, glm::vec3 &attenuation, Ray &scattered) const override
        {
            glm::vec3 reflected = reflect(glm::normalize(r_in.direction()), rec.normal);
            scattered = Ray(rec.p, reflected);
            attenuation = albedo;
            return (glm::dot(scattered.direction(), rec.normal) > 0);
        }
        glm::vec3 albedo;
        float fuzz;
    };

} // namespace rt