#ifndef PARTICLES_SOLVER_H
#define PARTICLES_SOLVER_H

#include "../includes/structs.h"
#include "raylib.h"
#include <stdint.h>

void particlesCollisionsCheck( struct Particles_t *particles, const uint16_t window_width, const uint16_t window_height, float time_step );
void particlesUpdate( struct Particles_t *particles, const float time_step );
void particlesCollisionsCheck_Grid( struct Particles_t *particles, const uint16_t window_width, const uint16_t window_height, float time_step );

#endif