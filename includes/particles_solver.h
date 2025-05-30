#ifndef PARTICLES_SOLVER_H
#define PARTICLES_SOLVER_H

#include "../includes/structs.h"
#include "raylib.h"
#include <stdint.h>

float energyCalc(uint32_t mass, Vector2 velocity);
void particlesCollisionsCheck(struct Particles_t *particles, const uint16_t window_width, const uint16_t window_height, float time_step);
void particlesUpdate(struct Particles_t *particles, const Vector2 acceleration, const float time_step);

#endif