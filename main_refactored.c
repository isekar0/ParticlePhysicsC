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
// #include <math.h>
#include <stdint.h>
#include "raylib.h"

// Keep the ratio as 16:9 and scale it up and down
#define WINDOW_SIZE_FACTOR 50 
#define WINDOW_TITLE "Particle simulation :3"

// Physical constants
#define GRAVITY 9.81
// #define AIR_RESISTENCE
// #define ELASTICITY_PARTICLE
// #define ELASTICITY_WALL


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

struct Particles_t particlesCreate(const size_t num_particles_init) {
    struct Particles_t particles;
    
    particles.mass = calloc(num_particles_init, sizeof(*particles.mass));
    particles.radius = calloc(num_particles_init, sizeof(*particles.radius));
    particles.colors = calloc(num_particles_init, sizeof(*particles.colors));
    particles.position = calloc(num_particles_init, sizeof(*particles.position));
    particles.velocity = calloc(num_particles_init, sizeof(*particles.velocity));

    if (particles.mass != NULL && particles.radius != NULL && particles.position != NULL && particles.velocity != NULL) {
        // Fill in metadata
        particles.number_of_stored_particles = 0;
        particles.arrays_size = num_particles_init;
        return particles;
    }
    else {
        fprintf(stderr, "*particles error");
        exit(EXIT_FAILURE);
    }
} 

void particlesFree(struct Particles_t *particles) {
    free(particles->mass);
    free(particles->radius);
    free(particles->colors);
    free(particles->position);
    free(particles->velocity);
    free(particles);
}

void particlesExpandArrays(struct Particles_t *particles) {
    size_t new_size = (int) particles->arrays_size * 2;
    struct Particles_t temp_particles; 

    temp_particles.mass = (uint32_t *)realloc(particles->mass, new_size);
    if (temp_particles.mass == NULL) {
        fprintf(stderr, "Error realloc particles.mass at size %lu", particles->arrays_size);
        particlesFree(particles);
        exit(EXIT_FAILURE);
    }
    
    temp_particles.radius = (uint32_t *)realloc(particles, new_size);
    if (temp_particles.mass == NULL) {
        fprintf(stderr, "Error realloc particles.radius at size %lu", particles->arrays_size);
        particlesFree(particles);
        exit(EXIT_FAILURE);
    }
    
    temp_particles.colors = (Color *)realloc((particles->colors), new_size);
    if (temp_particles.mass == NULL) {
        fprintf(stderr, "Error realloc particles.colors at size %lu", particles->arrays_size);
        particlesFree(particles);
        exit(EXIT_FAILURE);
    }

    temp_particles.position = (Vector2 *)realloc(particles->position, new_size);
    if (temp_particles.mass == NULL) {
        fprintf(stderr, "Error realloc particles.position at size %lu", particles->arrays_size);
        particlesFree(particles);
        exit(EXIT_FAILURE);
    }
    
    temp_particles.velocity = (Vector2 *)realloc(particles, new_size);
    if (temp_particles.mass == NULL) {
        fprintf(stderr, "Error realloc particles.velocity at size %lu", particles->arrays_size);
        particlesFree(particles);
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

void particlesAddParticle(struct Particles_t *particles, const uint32_t particle_i_mass, const uint32_t particle_i_radius) {

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
    particles->position[index].x = (float) GetMousePosition().x;
    particles->position[index].y = (float) GetMousePosition().y;
    particles->velocity[index].x = 0.0f;
    particles->velocity[index].y = 0.0f;
}

void particlesUpdate(struct Particles_t *particles, const Vector2 acceleration) {
    
    for (size_t idx = 0; idx < particles->number_of_stored_particles; idx++) {
        // update the velocities
        particles->velocity[idx].x += acceleration.x;
        particles->velocity[idx].y += acceleration. y; 
        // update the positions
        particles->position[idx].x += particles->velocity[idx].x;
        particles->position[idx].y += particles->velocity[idx].y;
    }
}

// Maybe integrate the raylib function for partilce-particle collisions
// bool CheckCollisionCircles(Vector2 center1, float radius1, Vector2 center2, float radius2);  
void particlesCollisionsCheck(struct Particles_t *particles, const uint16_t window_width, const uint16_t window_height) {
    
    for (size_t idx = 0; idx < particles->number_of_stored_particles; idx++) {

            if (particles->position[idx].x <= 0) {
                particles->position[idx].x = 1;
                particles->velocity[idx].x *= -1;
            }

            if (particles->position[idx].x >= GetScreenWidth()) {
                particles->position[idx].x = GetScreenWidth() - 1;
                particles->velocity[idx].x *= -1;
            }

            if (particles->position[idx].y <= 0) {
                particles->position[idx].y = 1;
                particles->velocity[idx].y *= -1;
            }

            if (particles->position[idx].y >= GetScreenHeight()) {
                particles->position[idx].y = GetScreenHeight();
                particles->velocity[idx].y *= -1;
            }
        }
}    

void particlesDraw(struct Particles_t *particles) {
    for (size_t idx = 0; idx < particles->number_of_stored_particles; idx++) {
        DrawCircleV(particles->position[idx], (float) particles->radius[idx], MAROON);
    }
}

int main(void) {
    // Initialize the raylib window 
    const uint16_t window_width = 16 * WINDOW_SIZE_FACTOR;
    const uint16_t window_height = 9 * WINDOW_SIZE_FACTOR;

    SetConfigFlags(FLAG_MSAA_4X_HINT);
    InitWindow(window_width, window_height, WINDOW_TITLE);

    unsigned int frames_counter = 0;

    SetTargetFPS(60);

    // Initialize 
    size_t num_particles_init = 5;
    struct Particles_t particles = particlesCreate(num_particles_init);
    struct Particles_t *p_particles = &particles;

    const Vector2 acceleration_field = { 0 , GRAVITY };

    // Main game loop
    while (!WindowShouldClose()) {
        // Update variables here

        // Check for left click to add particle, and right click to remove last particle added
        if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) {
            uint32_t mass = 10;
            uint32_t radius = 10;
            particlesAddParticle(p_particles, mass, radius);
        }

        // e.g. position += velocity
        particlesUpdate(p_particles, acceleration_field);
        
        // Check collisions
        particlesCollisionsCheck(p_particles, window_width, window_height);
        
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


    // dummy test loop
    // const uint16_t window_width = 16 * WINDOW_SIZE_FACTOR;
    // const uint16_t window_height = 9 * WINDOW_SIZE_FACTOR;

    //  // test static array
    // Vector2 positions[2] = {{200.0f, 100.0f}, {250.0f, 150.0f}};
    // Vector2 velocities[2] = { {5.0f, 0.0f}, {2.0f, 3.0f} };

    // SetConfigFlags(FLAG_MSAA_4X_HINT);
    // InitWindow(window_width, window_height, WINDOW_TITLE);

    // unsigned int frames_counter = 0;

    // SetTargetFPS(60);

    // while (!WindowShouldClose()) {
    // //  Update variables here
    // //  e.g. position += velocity
    //     for (size_t i = 0; i < sizeof(positions)/sizeof(positions[0]); i++) {
    //         velocities[i].y += 9.81;

    //         positions[i].x += velocities[i].x;
    //         positions[i].y += velocities[i].y;
    //     }
        
    // //  Check collisions
    //     for (size_t i = 0; i < sizeof(positions)/sizeof(positions[0]); i++) {
    //         if (positions[i].x <= 0) {
    //             positions[i].x = 1;
    //             velocities[i].x *= -0.7;
    //         }
    //         if (positions[i].x >= window_width) {
    //             positions[i].x = window_width - 1;
    //             velocities[i].x *= -0.7;
    //         }
    //         if (positions[i].y <= 0) {
    //             positions[i].y = 1;
    //             velocities[i].x *= -0.7;
    //         }
    //         if (positions[i].y >= window_height) {
    //             positions[i].y = window_height + 1;
    //             velocities[i].y *= -0.7;
    //         }
    //     }
        
    //     frames_counter++;

    // //  Draw
    //     BeginDrawing(); {
    //         ClearBackground(RAYWHITE);

    //         //DrawCircleV(position, (float) radius, COLOR);
    //         for (size_t i = 0; i < sizeof(positions)/sizeof(positions[0]); i++) {
    //             DrawCircleV(positions[i], 10, MAROON);
    //         }
            
    //         DrawFPS(10, 10);

    //         } 
    //     EndDrawing();
    // }

    // CloseWindow();


    return 0;
}
