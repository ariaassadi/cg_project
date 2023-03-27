#include "rt_raytracing.h"
#include "rt_ray.h"
#include "rt_hitable.h"
#include "rt_sphere.h"
#include "rt_triangle.h"
#include "rt_box.h"
#include "rtweekend.h"
#include "material.h"

#include "cg_utils2.h" // Used for OBJ-mesh loading
#include <stdlib.h>    // Needed for drand48()

namespace rt
{

    // Store scene (world) in a global variable for convenience
    struct Scene
    {
        Sphere ground;
        std::vector<Sphere> spheres;
        std::vector<Box> boxes;
        std::vector<Triangle> mesh;
        Box mesh_bbox;
        Sphere mesh_bsphere;
    } g_scene;

    // struct BoundingSphere
    // {
    //     glm::vec3 center;
    //     float radius;
    // } g_bounding_sphere;

    bool hit_world(RTContext &rtx, const Ray &r, float t_min, float t_max, HitRecord &rec)
    {
        HitRecord temp_rec;
        bool hit_anything = false;
        float closest_so_far = t_max;

        if (g_scene.ground.hit(r, t_min, closest_so_far, temp_rec))
        {
            hit_anything = true;
            closest_so_far = temp_rec.t;
            rec = temp_rec;
        }
        for (int i = 0; i < g_scene.spheres.size(); ++i)
        {
            if (g_scene.spheres[i].hit(r, t_min, closest_so_far, temp_rec))
            {
                hit_anything = true;
                closest_so_far = temp_rec.t;
                rec = temp_rec;
            }
        }
        for (int i = 0; i < g_scene.boxes.size(); ++i)
        {
            if (g_scene.boxes[i].hit(r, t_min, closest_so_far, temp_rec))
            {
                hit_anything = true;
                closest_so_far = temp_rec.t;
                rec = temp_rec;
            }
        }
        if (rtx.bsphere)
        {
            if (g_scene.mesh_bsphere.hit(r, t_min, closest_so_far, temp_rec))
            {
                for (int i = 0; i < g_scene.mesh.size(); ++i)
                {
                    if (g_scene.mesh[i].hit(r, t_min, closest_so_far, temp_rec))
                    {
                        hit_anything = true;
                        closest_so_far = temp_rec.t;
                        rec = temp_rec;
                    }
                }
            }
        }
        else if (rtx.bbox)
        {
            if (g_scene.mesh_bbox.hit(r, t_min, closest_so_far, temp_rec))
            {
                for (int i = 0; i < g_scene.mesh.size(); ++i)
                {
                    if (g_scene.mesh[i].hit(r, t_min, closest_so_far, temp_rec))
                    {
                        hit_anything = true;
                        closest_so_far = temp_rec.t;
                        rec = temp_rec;
                    }
                }
            }
        }
        else
        {
            for (int i = 0; i < g_scene.mesh.size(); ++i)
            {
                if (g_scene.mesh[i].hit(r, t_min, closest_so_far, temp_rec))
                {
                    hit_anything = true;
                    closest_so_far = temp_rec.t;
                    rec = temp_rec;
                }
            }
        }

        // bool spehere_hit = g_scene.spheres[0].hit(r, t_min, closest_so_far, temp_rec);
        // if (spehere_hit)
        // {
        //     // Check if bunny is hit
        //     for (int i = 0; i < g_scene.mesh.size(); ++i)
        //     {
        //         if (g_scene.mesh[i].hit(r, t_min, closest_so_far, temp_rec))
        //         {
        //             hit_anything = true;
        //             closest_so_far = temp_rec.t;
        //             rec = temp_rec;
        //         }
        //     }
        // }
        return hit_anything;
    }

    // This function should be called recursively (inside the function) for
    // bouncing rays when you compute the lighting for materials, like this
    //
    // if (hit_world(...)) {
    //     ...
    //     return color(rtx, r_bounce, max_bounces - 1);
    // }
    //
    // See Chapter 7 in the "Ray Tracing in a Weekend" book
    glm::vec3 color(RTContext &rtx, const Ray &r, int max_bounces)
    {
        if (max_bounces < 0)
            return glm::vec3(0.0f);

        HitRecord rec;
        if (hit_world(rtx, r, 0.001f, infinity, rec))
        {
            rec.normal = glm::normalize(rec.normal); // Always normalise before use!
            // rec.p = glm::normalize(rec.p);           // Always normalise before use!
            if (rtx.show_normals)
            {
                return rec.normal * 0.5f + 0.5f;
            }

            // Implement lighting for materials here
            Ray scattered;
            glm::vec3 attenuation;
            if (rec.mat_ptr->scatter(r, rec, attenuation, scattered))
            {
                return attenuation * color(rtx, scattered, max_bounces - 1);
            }
            return glm::vec3(0.0f);

            // glm::vec3 target = rec.p + rec.normal + glm::normalize(random_in_unit_sphere());
            // return 0.5f * color(rtx, Ray(rec.p, target - rec.p), max_bounces - 1);
        }

        // If no hit, return sky color
        glm::vec3 unit_direction = glm::normalize(r.direction());
        float t = 0.5f * (unit_direction.y + 1.0f);
        return (1.0f - t) * rtx.ground_color + t * rtx.sky_color;
    }

    // The function average_of_vertices() is used to compute the center of a mesh
    // (for the bounding box)
    glm::vec3 average_of_vertices(const std::vector<glm::vec3> &vertices)
    {
        glm::vec3 sum(0.0f);
        for (int i = 0; i < vertices.size(); ++i)
            sum += vertices[i];
        return sum / (float)vertices.size();
    }

    // The point furthers away from the center of the mesh
    glm::vec3 furthest_point(const std::vector<glm::vec3> &vertices, const glm::vec3 &center)
    {
        glm::vec3 furthest = vertices[0];
        float max_dist = glm::distance(furthest, center);
        for (int i = 1; i < vertices.size(); ++i)
        {
            float dist = glm::distance(vertices[i], center);
            if (dist > max_dist)
            {
                max_dist = dist;
                furthest = vertices[i];
            }
        }
        return furthest;
    }

    // MODIFY THIS FUNCTION!
    void setupScene(RTContext &rtx, const char *filename)
    {
        g_scene.ground = Sphere(glm::vec3(0.0f, -1000.5f, 0.0f), 1000.0f, new Lambertian(glm::vec3(0.1f, 0.5f, 0.1f)));
        // g_scene.spheres = {
        //     // Sphere(glm::vec3(0.0f, 0.2f, 0.0f), 1.1f, new HitBox(0.0f)),
        //     // Sphere(glm::vec3(-0.3f, -0.4f, 0.9f), 0.1f, new Metal(glm::vec3(0.8f, 0.8f, 0.8f), 0.6f)),
        //     Sphere(glm::vec3(0.0f, 0.0f, 0.0f), 0.5f, new Metal(glm::vec3(0.8f, 0.8f, 0.8f), 0.6f)),
        //     // Sphere(glm::vec3(-1.0f, -0.4f, 0.6f), 0.1f, new Lambertian(glm::vec3(0.8f, 0.8f, 0.8f))),
        //     // Sphere(glm::vec3(0.5f, -0.3f, 0.7f), 0.2f, new Metal(glm::vec3(0.8f, 0.8f, 0.8f), 0.0f)),
        // };
        cg::OBJMesh mesh;
        cg::objMeshLoad(mesh, filename);
        g_scene.mesh.clear();
        for (int i = 0; i < mesh.indices.size(); i += 3)
        {
            int i0 = mesh.indices[i + 0];
            int i1 = mesh.indices[i + 1];
            int i2 = mesh.indices[i + 2];
            glm::vec3 v0 = mesh.vertices[i0] + glm::vec3(0.0f, 0.135f, 0.0f);
            glm::vec3 v1 = mesh.vertices[i1] + glm::vec3(0.0f, 0.135f, 0.0f);
            glm::vec3 v2 = mesh.vertices[i2] + glm::vec3(0.0f, 0.135f, 0.0f);
            g_scene.mesh.push_back(Triangle(v0, v1, v2, rtx.mat_ptr));
        }
        glm::vec3 center_of_mesh = average_of_vertices(mesh.vertices);
        glm::vec3 furthest = furthest_point(mesh.vertices, center_of_mesh);
        float radius = glm::distance(furthest, center_of_mesh);

        g_scene.mesh_bbox = Box(center_of_mesh, glm::vec3(radius), new Metal(glm::vec3(0.8f, 0.8f, 0.8f), 0.0f));
        g_scene.mesh_bsphere = Sphere(center_of_mesh, radius, new Lambertian(glm::vec3(0.8f, 0.8f, 0.8f)));

        // g_scene.boxes = {g_scene.mesh_bbox};
        // g_scene.spheres = {g_scene.mesh_bsphere};
    }

    // MODIFY THIS FUNCTION!
    void updateLine(RTContext &rtx, int y)
    {
        int nx = rtx.width;
        int ny = rtx.height;
        float aspect = float(nx) / float(ny);
        glm::vec3 lower_left_corner(-1.0f * aspect, -1.0f, -1.0f);
        glm::vec3 horizontal(2.0f * aspect, 0.0f, 0.0f);
        glm::vec3 vertical(0.0f, 2.0f, 0.0f);
        glm::vec3 origin(0.0f, 0.0f, 0.0f);
        glm::mat4 world_from_view = glm::inverse(rtx.view);

        if (rtx.reflective)
        {
            for (int i = 0; i < g_scene.mesh.size(); i += 1)
            {
                g_scene.mesh.at(i).change_material(new Metal(glm::vec3(0.8f, 0.8f, 0.8f), 0.0f));
            }
        }
        else
        {
            for (int i = 0; i < g_scene.mesh.size(); i += 1)
            {
                g_scene.mesh.at(i).change_material(new Lambertian(glm::vec3(0.8f, 0.6f, 0.1f)));
            }
        }

        // You can try parallelising this loop by uncommenting this line:
        // #pragma omp parallel for schedule(dynamic)
        for (int x = 0; x < nx; ++x)
        {
            float u = 0;
            float v = 0;
            if (rtx.antialiasing_jitter)
            {
                v = (float(y) + random_double()) / float(ny - 1);
                u = (float(x) + random_double()) / float(nx - 1);
            }
            else
            {
                u = (float(x) + 0.5f) / float(nx - 1);
                v = (float(y) + 0.5f) / float(ny - 1);
            }

            Ray r(origin, lower_left_corner + u * horizontal + v * vertical);
            r.A = glm::vec3(world_from_view * glm::vec4(r.A, 1.0f));
            r.B = glm::vec3(world_from_view * glm::vec4(r.B, 0.0f));

            // Note: in the RTOW book, they have an inner loop for the number of
            // samples per pixel. Here, you do not need this loop, because we want
            // some interactivity and accumulate samples over multiple frames
            // instead (until the camera moves or the rendering is reset).

            if (rtx.current_frame <= 0)
            {
                // Here we make the first frame blend with the old image,
                // to smoothen the transition when resetting the accumulation
                glm::vec4 old = rtx.image[y * nx + x];
                rtx.image[y * nx + x] = glm::clamp(old / glm::max(1.0f, old.a), 0.0f, 1.0f);
            }
            glm::vec3 c = color(rtx, r, rtx.max_bounces);
            rtx.image[y * nx + x] += glm::vec4(c, 1.0f);
        }
    }

    void updateImage(RTContext &rtx)
    {
        if (rtx.freeze)
            return;                               // Skip update
        rtx.image.resize(rtx.width * rtx.height); // Just in case...

        // if (rtx.cube_reflective)
        // {
        //     g_scene.boxes = {
        //         Box(glm::vec3(-1.0f, -0.25f, 0.0f), glm::vec3(0.25f), new Metal(glm::vec3(0.8f, 0.8f, 0.8f), 0.0f))};
        // }
        // else
        // {
        //     g_scene.boxes = {
        //         Box(glm::vec3(-1.0f, -0.25f, 0.0f), glm::vec3(0.25f), new Lambertian(glm::vec3(0.8f, 0.8f, 0.8f)))};
        // }
        updateLine(rtx, rtx.current_line % rtx.height);

        if (rtx.current_frame < rtx.max_frames)
        {
            rtx.current_line += 1;
            if (rtx.current_line >= rtx.height)
            {
                rtx.current_frame += 1;
                rtx.current_line = rtx.current_line % rtx.height;
            }
        }
    }

    void resetImage(RTContext &rtx)
    {
        rtx.image.clear();
        rtx.image.resize(rtx.width * rtx.height);
        rtx.current_frame = 0;
        rtx.current_line = 0;
        rtx.freeze = false;
    }

    void resetAccumulation(RTContext &rtx)
    {
        rtx.current_frame = -1;
    }

} // namespace rt
