#ifndef STRUCTS_H
#define STRUCTS_H

#include <stdint.h>
#include <stdlib.h>
#include <raylib.h>

struct Particles_t {
    // Metadata
    size_t number_of_stored_particles;
    size_t arrays_size;
    // fixed values
    uint32_t *mass;
    uint32_t *radius;
    // varying values
    Color *colors;
    Vector2 *position; // raylib 2D vector = { .x, .y }
    Vector2 *velocity;
    float *energy;
    float highest_energy;
};

struct Tuple_Vector2_t {
    Vector2 vec1, vec2;
};

#endif