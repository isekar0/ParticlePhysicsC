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

// Inlucdes
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include "raylib.h"

// Keep the ratio as 16:9 and scale it up and down
#define WINDOW_SIZE_FACTOR 50 
#define WINDOW_TITLE "Particle simulation :3"
// Array of particles
#define NUM_PARICLES_INIT 10
#define NUM_PARTICLES_MAX 20;
// Physical constants
#define GRAVITY -9.81
// #define AIR_RESISTENCE
// #define ELASTICITY_PARTICLE
// #define ELASTICITY_WALL


// Consider placing structures and functions relating to the physics of particles in a separate header
// Data structures
// Vector2D_f will store all vectors other than position
struct Vector2D_f {
    float x_component_f;
    float y_component_f;
};

struct Particles_t {
    // Metadata
    uint32_t number_of_stored_particles;
    uint32_t arrays_size;
    // constants
    uint32_t *mass;
    uint32_t *radius;
    // variables
    struct Vector2 *position; // raylib 2D vector = { .x, .y }
    struct Vector2D_f *velocity;
};

struct Particles_t *particlesCreate(void) {
    const size_t num_particles = NUM_PARICLES_INIT;

    // initialize struct, allocate enough memory for num_particles number of particles, but contains no actual particles
    struct Particles_t *particles = (struct Particles_t*) malloc(
    sizeof(struct Particles_t) +
    num_particles * (sizeof(typeof(particles->mass))) + sizeof(typeof(particles->radius)) +
    sizeof(typeof(particles->position)) + sizeof(typeof(particles->velocity)));

    if (particles != NULL && particles->mass != NULL && particles->radius != NULL && particles->position != NULL && particles->velocity != NULL) {
        // Fill in metadata
        particles->number_of_stored_particles = 0;
        particles->arrays_size = num_particles;
        return particles;
    }
    else {
        fprintf(stderr, "*particles error");
        exit(EXIT_FAILURE);
    }
} 

// void particlesReallocMemory(struct Particles_t *particles) {

// }

void particlesAddParticle(struct Particles_t *particles, const uint32_t particle_i_mass, const uint32_t particle_i_radius) {
    // update number of particles stored
    particles->number_of_stored_particles += 1;
    
    // Expand size of arrays if they are full 
    if (particles->number_of_stored_particles > particles->arrays_size) {
        // particlesReallocMemory(particles);
    }



}

void particlesUpdate(struct Particles_t *particles, struct Vector2D_f acceleration) {
    for (size_t idx = 0; idx < particles->number_of_stored_particles; idx++) {
        // update the velocities
        particles->velocity[idx].x_component_f += acceleration.x_component_f;
        particles->velocity[idx].y_component_f += acceleration. y_component_f; 
        // update the positions
        particles->position[idx].x += particles->velocity[idx].x_component_f;
        particles->position[idx].y += particles->velocity[idx].y_component_f;
    }
}

// void particlesCollisionsCheck(struct Particles_t *particles, const uint16_t window_width, const uint16_t window_height) {
    
// }
    


int main(void) {
    // Initialize the raylib window 
    const uint16_t window_width = 16 * WINDOW_SIZE_FACTOR;
    const uint16_t window_height = 9 * WINDOW_SIZE_FACTOR;

    SetConfigFlags(FLAG_MSAA_4X_HINT);
    InitWindow(window_width, window_height, WINDOW_TITLE);

    unsigned int frames_counter = 0;

    SetTargetFPS(60);

    // Initialize 
    // struct Particles_t *particles = particlesCreate();

    // Main game loop
    while (!WindowShouldClose()) {
    // Update variables here
    // e.g. position += velocity

    // Check collisions

    //frames_counter++;

    // Draw
    BeginDrawing(); {
        ClearBackground(RAYWHITE);

        //DrawCircleV(position, (float) radius, COLOR);
        
        DrawFPS(10, 10);

    } EndDrawing();
    }

    CloseWindow();

    return 0;
}
