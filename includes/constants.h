#ifndef CONSTANTS_H
#define CONSTANTS_H

// Metadata
// Keep the ratio as 16:9 and scale it up and down
#define WINDOW_SIZE_FACTOR 75
#define WINDOW_TITLE "Particle simulation :3"

// Solver constants
#define SOLVER_TUNER 0.3f // larger values can cause this phenomenon where extra energy/momentum is introduced to the system, values greater than 0.4 result is massive impulse spikes
#define EPSILON 1e-8f     // small value for division by zero correction
#define DEPTH_CAP 0.05f
#define SLOP_FACTOR 0.01f
#define BAUMGARTE 0.2f
#define GRAVITY 9.81f // Negated so it can be used correctly (positive makes particles go down)
#define AIR_RESISTENCE 0.2
#define RESTITUTION 0.9
#define DAMPING_WALL 0.9
#define TIMESTEP_FACTOR 5 // time_step /= TIMESTEP_FACTOR
#define ENERGY_MIN 1e-4
#define CELL_SIZE_FACTOR 2 // cell size = k * max_radius; 1 <= k <= 2

// Array constants
#define NUM_INIT_SIZE 50000
#define ARRAY_GROWTH_FACTOR 2

#endif