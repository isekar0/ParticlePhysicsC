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
#define PARTICLE_RADIUS 20
// #define AIR_RESISTENCE
// #define ELASTICITY_PARTICLE
// #define ELASTICITY_WALL



// Consider placing structures and functions relating to the physics of particles in a separate header
// Data structures
// Vector2D_f will store all vectors other than position
typedef struct Vector2D_f {
  float x_component_f;
  float y_component_f;
} Vector2D_f;

// Vector2D_i will store the position vector
typedef struct Vector2D_i {
  uint32_t x_component_i;
  uint32_t y_component_i;
} Vector2D_i;

typedef struct Particle {
  // Constants
  const uint32_t mass;
  const uint32_t radius; 
  // Variables to be updated per timeframe
  Vector2 position;
  Vector2D_f velocity;
  Vector2D_f acceleration;
} Particle;

void timestepUpdate(Particle *particles, const size_t num_particles) {
  for (size_t i = 0; i < num_particles; i++) {
    // update velocity
    particles[i].velocity.x_component_f += particles[i].acceleration.x_component_f;
    particles[i].velocity.y_component_f += particles[i].acceleration.y_component_f;

    // update position
    particles[i].position.x += particles[i].velocity.x_component_f;
    particles[i].position.y += particles[i].velocity.y_component_f;
  }
}

void collisionsResolve(Particle *particles, const size_t num_particles, const uint32_t window_width, const uint32_t window_height) {
  /*Check for collisions
   * Currently only works with particle-wall collisions
   * $Implment particle-particle collision checking
   */
  // check if pos_x - r < 0 || pos_x + r > width || pos_y - r < 0 || pos_y + r > height 
  const uint32_t bounds_x = window_width - PARTICLE_RADIUS; 
  const uint32_t bounds_y = window_height - PARTICLE_RADIUS;
  for (size_t idx = 0; idx < num_particles; idx++) {
    /* Iterate through particle in particles
     * Check x and y components of position and compare with wall
     * In case of collision
     * 1) Prevent out-of-bounds position
     * 2) Apply counter force and update velocity
     * */
    if (particles[idx].position.x < 0 || particles[idx].position.x > bounds_x) {
       particles[idx].velocity.x_component_f *= (float) -1.0; // * damping coefficient later
     }
    if (particles[idx].position.y < 0 || particles[idx].position.y > bounds_y) {
      particles[idx].velocity.y_component_f *= (float) -1.0;
    }
  }
}
/*
void particlesCreate(const uint32_t mass, const uint32_t radius, Particle *particles, size_t num_particles) {
 for (size_t idx = 0; idx < num_particles; idx++) {
   Vector2 idx_position = { (uint32_t) idx % (uint32_t) num_particles, (uint32_t) idx / (uint32_t) num_particles };
   Particle particle = { .mass = mass, .radius = radius, .position = idx_position, .velocity = {0, 0} };
   *particles[idx] = particle;
 }
}
*/
void particlesDraw(const Particle *particles, size_t num_particles) {
  for (size_t idx = 0; idx < num_particles; idx++) {
   DrawCircleV(particles[idx].position, particles[idx].radius, MAROON); 
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
  uint32_t num_particles = NUM_PARICLES_INIT;
  uint32_t particle_radius = PARTICLE_RADIUS;
  Particle *particles = malloc(sizeof(Particle) * num_particles);
  if (particles == NULL) {
    perror("*particles failed to initialize");
  }
  //particlesCreate(5 /*kg*/, particle_radius /*m*/, particles, num_particles);
  for (size_t idx = 0; idx < num_particles; idx++) {
    uint32_t x = (window_width / 2.0f) + particle_radius * (idx % num_particles);
    uint32_t y = (window_height / 2.0f) + particle_radius * (idx / num_particles);
    particles[idx] = (Particle[5]) { .mass = 5, .radius = particle_radius, .position = { x, y }, .velocity = {0.0f, 0.0f}, .acceleration = { 0.0f, 0.0f} };
  }
  // Main game loop
  while (!WindowShouldClose()) {
    // Update variables here
    // e.g. position += velocity
    timestepUpdate(particles, num_particles);

    // Check collisions
    collisionsResolve(particles, num_particles, window_width, window_height);

    //frames_counter++;

    // Draw
    BeginDrawing(); {
      ClearBackground(RAYWHITE);

      //DrawCircleV(position, (float) radius, COLOR);
      particlesDraw(particles, num_particles); 
      DrawFPS(10, 10);

    } EndDrawing();
  }
  
  free(particles);
  CloseWindow();

  return 0;
}
