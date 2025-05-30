#ifndef PARTICLES_ARRAYS_H
#define PARTICLES_ARRAYS_H

#include <stdlib.h>
#include <stdint.h>
#include "raylib.h"

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

struct Particles_t;

void particlesFree(struct Particles_t *particles); 
struct Particles_t particlesCreate(size_t num_particles_init); 
void particlesAddParticle(struct Particles_t *particles, const Vector2 particle_i_position, const uint32_t particle_i_mass, const uint32_t particle_i_radius, uint32_t num_particles);
void particlesRemoveParticle(struct Particles_t *particles);
void particlesDraw(struct Particles_t *particles);

#endif