#ifndef STRUCTS_H
#define STRUCTS_H

#include <raylib.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>

struct Particles_t {
    // varying values
    Color   *colors;
    Vector2 *position; // raylib 2D vector = { .x, .y }
    Vector2 *velocity;
    // fixed values
    uint32_t *mass;
    uint32_t *radius;
    // Metadata
    size_t length;
    size_t size;
};

struct Tuple_Vector2_t {
    Vector2 vec1, vec2;
};

// A helper struct for sorting by minX
struct BoundEntry {
    int   *idx;  // original particle index
    float *minX; // position[idx].x - radius[idx]
    float *maxX; // position[idx].x + radius[idx]
    size_t length;
    // size_t size;
};

#endif