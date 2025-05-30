#include "../includes/particles_solver.h"

#include "../includes/structs.h"
#include "../includes/particles_arrays.h"
#include "../includes/utils.h"
#include <math.h>
#include <stdio.h>



// Solver constants
#define SOLVER_TUNER 0.3f // larger values can cause this phenomenon where extra energy/momentum is introduced to the system, values greater than 0.4 result is massive impulse spikes
#define EPSILON 1e-8f // small value for division by zero correction
#define DEPTH_CAP 0.1f
#define SLOP_FACTOR 0.005f
#define BAUMGARTE 0.2f
#define GRAVITY 9.81 // Negated so it can be used correctly (positive makes particles go down)
#define AIR_RESISTENCE 0.98
#define RESTITUTION 0.90
#define DAMPING_WALL 0.90
#define TIMESTEP_FACTOR 5 // time_step /= TIMESTEP_FACTOR
#define ENERGY_MIN 1e-4

// Array constants
#define NUM_INIT_SIZE 20
#define ARRAY_GROWTH_FACTOR 2

float energyCalc(uint32_t mass, Vector2 velocity) {
    return 0.5 * mass * (velocity.x * velocity.x + velocity.y * velocity.y);
}


// https://stackoverflow.com/a/35212639
// https://en.wikipedia.org/wiki/Elastic_collision
struct Tuple_Vector2_t particlesResolveVelocity(uint32_t m1, uint32_t m2, Vector2 v1, Vector2 v2, Vector2 p1, Vector2 p2, uint32_t r1, uint32_t r2, float time_step) {
    
    // Euclidean distance between the two particles
    float euclid_distance_squared = (p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y);
    euclid_distance_squared = fmaxf(euclid_distance_squared, EPSILON); // in case it is too close to zero
    float euclid_distance = sqrtf(euclid_distance_squared);
    // inverse of each mass
    float m1_inverse = 1.0f / m1;
    float m2_inverse = 1.0f / m2;

    // the contact normal between the two particles
    Vector2 contact_normal;
    if (euclid_distance < EPSILON) {
        euclid_distance = 1;
        // default values, change to be random values on the unit circle (random but still of length 1)
        
        contact_normal.x = 1.0f;
        contact_normal.y = 0.0f;
        // int random_angle = (float) rand() / (float) (RAND_MAX / (2 * PI));
        // contact_normal.x = 1.0f / cos(random_angle);
        // contact_normal.y = 1.0f / sin(random_angle);
    }
    else {
        contact_normal.x = (p1.x - p2.x) / euclid_distance;
        contact_normal.y = (p1.y - p2.y) / euclid_distance;
    }
    
    // positional correction if there is some penetration
    float velocity_correction_term = 0; // set to non zero if there is penetration
    float penetration_depth_corrected;
    float penetration_depth = r1 + r2 - euclid_distance;

    if (penetration_depth > 0.0f) {

        float penetration_depth_max = DEPTH_CAP * fminf(r1, r2);
        penetration_depth_corrected = penetration_depth * SOLVER_TUNER;
        penetration_depth_corrected = fmaxf(-penetration_depth_max, fminf(penetration_depth_corrected, penetration_depth_max));

        float p1_factor = (m1_inverse / (m1_inverse + m2_inverse)) * penetration_depth_corrected;
        float p2_factor = (m2_inverse / (m1_inverse + m2_inverse)) * penetration_depth_corrected;

        // correct positions
        p1.x += p1_factor * contact_normal.x;
        p1.y += p1_factor * contact_normal.y;
        p2.x -= p2_factor * contact_normal.x;
        p2.y -= p2_factor * contact_normal.y;

        // recompute
        euclid_distance = sqrtf(fmaxf((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y), EPSILON));
        contact_normal.x = (p1.x - p2.x) / euclid_distance;
        contact_normal.y = (p1.y - p2.y) / euclid_distance;

        float residual_depth = penetration_depth - penetration_depth_corrected;

        // set velocity_correction_term to a non zero value if the remaining depth fulfills the condition
        float slop = SLOP_FACTOR * fmaxf(r1, r2);
        if (residual_depth > slop) {
            velocity_correction_term = BAUMGARTE * (residual_depth - slop) / time_step;
        }
    }
    
    // Component of relative velocity between the particles will also be the same
    Vector2 relative_velocity = { v1.x - v2.x, v1.y - v2.y };
    float velocity_normal = relative_velocity.x * contact_normal.x + relative_velocity.y * contact_normal.y;

    if (velocity_normal > 0.0f && velocity_correction_term == 0.0f) {
        struct Tuple_Vector2_t old_velocities = { v1, v2 };
        return old_velocities;
    }

    float velocity_normal_desired = -RESTITUTION * velocity_normal + velocity_correction_term;

    float scalar_impulse = (velocity_normal_desired - velocity_normal) / (m1_inverse + m2_inverse);

    Vector2 vector_impulse = { scalar_impulse * contact_normal.x, scalar_impulse * contact_normal.y };

    Vector2 v1_prime = {v1.x + vector_impulse.x * m1_inverse, v1.y + vector_impulse.y * m1_inverse };
    Vector2 v2_prime = { v2.x - vector_impulse.x * m2_inverse, v2.y - vector_impulse.y * m2_inverse };
    
    struct Tuple_Vector2_t new_velocities = { v1_prime, v2_prime };

    return new_velocities;
}

// Maybe integrate the raylib function for partilce-particle collisions
// bool CheckCollisionCircles(Vector2 center1, float radius1, Vector2 center2, float radius2);  
void particlesCollisionsCheck(struct Particles_t *particles, const uint16_t window_width, const uint16_t window_height, float time_step) {
    
    for (size_t idx = 0; idx < particles->number_of_stored_particles; idx++) {

        if (particles->position[idx].x <= 0) {
            particles->position[idx].x = 1;
            particles->velocity[idx].x *= -1 * DAMPING_WALL;
        }

        if (particles->position[idx].x >= window_width) {
            particles->position[idx].x = window_width - 1;
            particles->velocity[idx].x *= -1 * DAMPING_WALL;
        }

        if (particles->position[idx].y < 0) {
            particles->position[idx].y = 1;
            particles->velocity[idx].y *= -1 * DAMPING_WALL;
            // particles->velocity[idx].y += particles->mass[idx] * GRAVITY;
        }

        if (particles->position[idx].y >= window_height) {
            particles->position[idx].y = window_height - 1;
            particles->velocity[idx].y *= -1 * DAMPING_WALL;
        }
    }

    for (size_t idx = 0; idx < particles->number_of_stored_particles - 1; idx++) {
        Vector2 position_idx = particles->position[idx];
        uint32_t radius_idx = particles->radius[idx];
        uint32_t mass_idx = particles->mass[idx];
        Vector2 velocity_idx = particles->velocity[idx];
        for (size_t jdx = idx + 1; jdx < particles->number_of_stored_particles; jdx++) {
            if (CheckCollisionCircles(position_idx, radius_idx, particles->position[jdx], particles->radius[jdx])) {
                Vector2 position_jdx = particles->position[jdx];
                uint32_t mass_jdx = particles->mass[jdx];
                uint32_t radius_jdx = particles->radius[jdx];
                Vector2 velocity_jdx = particles->velocity[jdx];
                // solve for final velocities
                struct Tuple_Vector2_t v1_v2_new = particlesResolveVelocity(mass_idx, mass_jdx, velocity_idx, velocity_jdx, position_idx, position_jdx, radius_idx, radius_jdx, time_step);
                particles->velocity[idx] = v1_v2_new.vec1;
                particles->velocity[jdx] = v1_v2_new.vec2;
            }
        }
    }
}    

void particlesUpdate(struct Particles_t *particles, const Vector2 acceleration, const float time_step) {
    float highest_energy_iteration = 0;
    for (size_t idx = 0; idx < particles->number_of_stored_particles; idx++) {
        // update the positions
        particles->position[idx].x += 0.5 * acceleration.x * time_step * time_step + particles->velocity[idx].x * time_step;
        particles->position[idx].y += 0.5 * acceleration.y * time_step * time_step + particles->velocity[idx].y * time_step;
    
        // update the velocities
        particles->velocity[idx].x += acceleration.x * time_step * AIR_RESISTENCE; 
        particles->velocity[idx].y += acceleration.y * time_step * AIR_RESISTENCE;

        // update the energy
        particles->energy[idx] = energyCalc(particles->mass[idx], particles->velocity[idx]);

        float energy_idx = particles->energy[idx];
        if (energy_idx > highest_energy_iteration) {
            highest_energy_iteration = energy_idx;
        }
        
    }

    particles->highest_energy = fmaxf(particles->highest_energy * 0.95, highest_energy_iteration);
    particles->highest_energy = fmaxf(particles->highest_energy, 50);

    // update the colors
    for (size_t idx = 0; idx < particles->number_of_stored_particles; idx++) {
        // float temp = (particles->energy[idx] + ENERGY_MIN) / (particles->highest_energy - ENERGY_MIN); 
        float temp = (particles->energy[idx] + ENERGY_MIN) / (100000 - ENERGY_MIN); 

        temp = clamp_f(temp, 0, 1);
        float nth_root = 3.0f;
        temp = pow(temp, 1.0f / nth_root);

        // const float TAU = 0.05f;
        // float temp_shoulder = (temp - TAU) / (1 - TAU);
        
        // float temp_biased = clamp_f(temp_shoulder, 0.0f, 1.0f); 
        
        particles->colors[idx] = TransformH(RED, 360 * (1 - temp));
        particles->colors[idx].a = 255;
    }
}

