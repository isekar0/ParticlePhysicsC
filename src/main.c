/*
 * =====================================================================================
 *
 *       Filename:  main.c
 *
 *    Description: Basic fluid/particle simulation in C using Raylib for rendering 
 *
 *        Version:  1.1
 *        Created:  05/11/2025 03:59:13 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Martin Hanna
 *   Organization:  N/A
 *
 * =====================================================================================
 */

// Includes
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <sys/types.h>
#include "raylib.h"

// Metadata
// Keep the ratio as 16:9 and scale it up and down
#define WINDOW_SIZE_FACTOR 50 
#define WINDOW_TITLE "Particle simulation :3"
// Testing
#define NUM_INIT_SIZE 20
#define ARRAY_GROWTH_FACTOR 2

// Physical constants
// Negated so it can be used correctly (positive makes particles go down)
#define GRAVITY 9.81
// #define AIR_RESISTENCE
#define DAMPING_PARTICLE 0.98
#define DAMPING_WALL 0.98




// Consider placing structures and functions relating to the physics of particles in a separate header
// Data structures
struct Particles_t {
    // Metadata
    size_t number_of_stored_particles;
    size_t arrays_size;
    // constants
    uint32_t *mass;
    uint32_t *radius;
    Color *colors;
    // variables
    Vector2 *position; // raylib 2D vector = { .x, .y }
    Vector2 *velocity;
};

struct Tuple_Vector2_t {
    Vector2 vec1, vec2;
};

void particlesFree(struct Particles_t *particles) {
    free(particles->mass);
    free(particles->radius);
    free(particles->colors);
    free(particles->position);
    free(particles->velocity);
    // free(particles);
}

struct Particles_t particlesCreate(size_t num_particles_init) {
    struct Particles_t particles;
    particles.mass = calloc(num_particles_init, sizeof(*particles.mass));
    if (!particles.mass) {
        fprintf(stderr, "Failed to initialize particles");
        particlesFree(&particles);
        exit(EXIT_FAILURE);
    }
    particles.radius = calloc(num_particles_init, sizeof(*particles.radius));
    if (!particles.radius) {
        fprintf(stderr, "Failed to initialize particles");
        particlesFree(&particles);
        exit(EXIT_FAILURE);
    }
    particles.colors = calloc(num_particles_init, sizeof(*particles.colors));
    if (!particles.colors) {
        fprintf(stderr, "Failed to initialize particles");
        particlesFree(&particles);
        exit(EXIT_FAILURE);
    }
    particles.position = calloc(num_particles_init, sizeof(*particles.position));
    if (!particles.position) {
        fprintf(stderr, "Failed to initialize particles");
        particlesFree(&particles);
        exit(EXIT_FAILURE);
    }
    particles.velocity = calloc(num_particles_init, sizeof(*particles.velocity));
    if (!particles.velocity) {
        fprintf(stderr, "Failed to initialize particles");
        particlesFree(&particles);
        exit(EXIT_FAILURE);
    }

    // Fill in metadata
    particles.number_of_stored_particles = 0;
    particles.arrays_size = num_particles_init;
    return particles;
} 



void particlesExpandArrays(struct Particles_t *particles) {
    size_t new_size = ((int) particles->arrays_size ) * ARRAY_GROWTH_FACTOR;
    printf("Expand array has been called, current size is %lu and new size will be %lu", particles->arrays_size, new_size);
    struct Particles_t temp_particles = *particles; 

    temp_particles.mass = realloc(temp_particles.mass, new_size * sizeof(*particles->mass));
    if (!temp_particles.mass) {
        fprintf(stderr, "Failed to initialize particles");
        particlesFree(&temp_particles);
        exit(EXIT_FAILURE);
    }
    temp_particles.radius = realloc(temp_particles.radius, new_size * sizeof(*particles->radius));
    if (!temp_particles.radius) {
        fprintf(stderr, "Failed to initialize particles");
        particlesFree(&temp_particles);
        exit(EXIT_FAILURE);
    }
    temp_particles.colors = realloc(temp_particles.colors, new_size * sizeof(*particles->colors));
    if (!temp_particles.colors) {
        fprintf(stderr, "Failed to initialize particles");
        particlesFree(&temp_particles);
        exit(EXIT_FAILURE);
    }
    temp_particles.position = realloc(temp_particles.position, new_size * sizeof(*particles->position));
    if (!temp_particles.position) {
        fprintf(stderr, "Failed to initialize particles");
        particlesFree(&temp_particles);
        exit(EXIT_FAILURE);
    }
    temp_particles.velocity = realloc(temp_particles.velocity, new_size * sizeof(*particles->velocity));
    if (!temp_particles.velocity) {
        fprintf(stderr, "Failed to initialize particles");
        particlesFree(&temp_particles);
        particlesFree(particles); // its joever
        exit(EXIT_FAILURE);
    }

    // once all tests pass, set particles to be the temp
    particles->arrays_size = new_size;
    particles->mass = temp_particles.mass;
    particles->radius = temp_particles.radius;
    particles->colors = temp_particles.colors;
    particles->position = temp_particles.position;
    particles->velocity = temp_particles.velocity;
}

void particlesAddParticle(struct Particles_t *particles, const Vector2 particle_i_position, const uint32_t particle_i_mass, const uint32_t particle_i_radius) {

    // Expand size of arrays if they are full 
    if (particles->number_of_stored_particles >= particles->arrays_size) {
        particlesExpandArrays(particles);
    }
    
    // update number of particles stored
    particles->number_of_stored_particles += 1;

    // add new particle
    size_t index = particles->number_of_stored_particles;
    particles->mass[index] = particle_i_mass;
    particles->radius[index] = particle_i_radius;;
    particles->position[index].x = (float) particle_i_position.x;
    particles->position[index].y = (float) particle_i_position.y;
    particles->velocity[index].x = 0.0f;
    particles->velocity[index].y = 0.0f;
}

void particlesUpdate(struct Particles_t *particles, const Vector2 acceleration) {
    
    for (size_t idx = 0; idx < particles->number_of_stored_particles; idx++) {
        // update the velocities
        particles->velocity[idx].x += acceleration.x;
        particles->velocity[idx].y += acceleration.y; 
        // update the positions
        particles->position[idx].x += particles->velocity[idx].x;
        particles->position[idx].y += particles->velocity[idx].y;
    }
}

// https://stackoverflow.com/a/35212639
// https://en.wikipedia.org/wiki/Elastic_collision
struct Tuple_Vector2_t particlesResolveVelocity(uint32_t m1, uint32_t m2, Vector2 v1, Vector2 v2, Vector2 p1, Vector2 p2) {
    // resolve first velocity
    Vector2 v1_diff_v = { v1.x - v2.x, v1.y - v2.y};
    Vector2 v1_diff_p = { p1.x - p2.x, p1.y - p2.y};
    float v1_mag_p = v1_diff_p.x*v1_diff_p.x + v1_diff_p.y*v1_diff_p.y;
    v1_mag_p = powf(v1_mag_p, 0.5);
    float v1_dot_prod = v1_diff_v.x * v1_diff_p.x + v1_diff_v.y * v1_diff_p.y;
    uint32_t v1_frac_m1m2 = 2 * m1 / (m1 + m2);
    float v1_frac_vx = v1_dot_prod / v1_mag_p;
    Vector2 v1_diff_v1 = { v1_diff_p.x * v1_frac_vx * v1_frac_m1m2, v1_diff_p.y * v1_frac_vx * v1_frac_m1m2 };
    Vector2 v1_prime = { v1.x - v1_diff_v1.x, v1.y - v1_diff_v1.y };
    v1_prime.x *= DAMPING_PARTICLE;
    v1_prime.y *= DAMPING_PARTICLE;

    // now second, just swap 1 and 2 around
    // you can probably optimise this, a - b = - (b - a)
    Vector2 v2_diff_v = { v2.x - v1.x, v2.y - v1.y };
    Vector2 v2_diff_p = { p2.x - p1.x, p2.y - p1.y};
    float v2_mag_p = v2_diff_p.x*v2_diff_p.x + v2_diff_p.y*v2_diff_p.y;
    v2_mag_p = powf(v2_mag_p, 0.5);
    float v2_dot_prod = v2_diff_v.x * v2_diff_p.x + v2_diff_v.y * v2_diff_p.y;
    uint32_t v2_frac_m2m1 = 2 * m2 / (m2 + m1);
    float v2_frac_vx = v2_dot_prod / v2_mag_p;
    Vector2 v2_diff_v2 = { v2_diff_p.x * v2_frac_vx * v2_frac_m2m1, v2_diff_p.y * v2_frac_vx * v2_frac_m2m1 };
    Vector2 v2_prime = { v1.x - v2_diff_v2.x, v1.y - v2_diff_v2.y };
    v2_prime.x *= DAMPING_PARTICLE;
    v2_prime.y *= DAMPING_PARTICLE;

    // return the values
    struct Tuple_Vector2_t resolved_velocities = { v1_prime, v2_prime };
    return resolved_velocities;
}

// Maybe integrate the raylib function for partilce-particle collisions
// bool CheckCollisionCircles(Vector2 center1, float radius1, Vector2 center2, float radius2);  
void particlesCollisionsCheck(struct Particles_t *particles, const uint16_t window_width, const uint16_t window_height) {
    
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
                uint32_t radius_jdx = particles->radius[jdx];
                uint32_t mass_jdx = particles->mass[jdx];
                Vector2 velocity_jdx = particles->velocity[jdx];
                // solve for final velocities
                struct Tuple_Vector2_t v1_v2_new = particlesResolveVelocity(mass_idx, mass_jdx, velocity_idx, velocity_jdx, position_idx, position_jdx);
                particles->velocity[idx] = v1_v2_new.vec1;
                particles->velocity[jdx] = v1_v2_new.vec2;
            }
        }
    }
}    

void particlesDraw(struct Particles_t *particles) {
    for (size_t idx = 0; idx < particles->number_of_stored_particles; idx++) {
        DrawCircleV(particles->position[idx], (float) particles->radius[idx], MAROON);
    }
}

int main(void) {
    FILE *p_file = fopen("log.txt", "a");
    // Initialize the raylib window 
    const uint16_t window_width = 16 * WINDOW_SIZE_FACTOR;
    const uint16_t window_height = 9 * WINDOW_SIZE_FACTOR;

    SetConfigFlags(FLAG_MSAA_4X_HINT);
    InitWindow(window_width, window_height, WINDOW_TITLE);

    unsigned int frames_counter = 0;

    SetTargetFPS(10);

    // Initialize 
    size_t num_particles_init = NUM_INIT_SIZE; // WILL CRASH IF NOT MULTIPLE OF 20
    struct Particles_t particles = particlesCreate(num_particles_init);
    struct Particles_t *p_particles = &particles;

    const Vector2 acceleration_field = { 0 , GRAVITY };

    int test_counter = 0;
    int test_limit = 1000;
    // Main game loop
    while (!WindowShouldClose() && test_counter < test_limit) {
        // Update variables here

        // Check for left click to add particle, and right click to remove last particle added
        if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) {
            Vector2 mouse_position = GetMousePosition();
            uint32_t mass = 10;
            uint32_t radius = 10;
            particlesAddParticle(p_particles, mouse_position, mass, radius);
        }

        // Vector2 test_position_i = {test_counter, test_counter};
        // particlesAddParticle(p_particles, test_position_i, 10, 10);
        // test_counter += 1;
        // fprintf(p_file, "reached test_counter %d", test_counter);

        // e.g. position += velocity
        particlesUpdate(p_particles, acceleration_field);
        
        // Check collisions
        if (p_particles->number_of_stored_particles > 0) {
            particlesCollisionsCheck(p_particles, window_width, window_height);
        }
        
        frames_counter++;

        // Draw
        BeginDrawing(); {
            ClearBackground(RAYWHITE);

            //DrawCircleV(position, (float) radius, COLOR);
            particlesDraw(p_particles);
            
            DrawFPS(10, 10);

            } 
        EndDrawing();
    }

    particlesFree(p_particles);
    CloseWindow();
    fclose(p_file);

    return 0;
}
