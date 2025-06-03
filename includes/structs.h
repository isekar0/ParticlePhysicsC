#ifndef STRUCTS_H
#define STRUCTS_H

#include <raylib.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>

struct Particles_t {
    // varying values
    // Vector2 *velocity;
    // Vector2 *position; // raylib 2D vector = { .x, .y }
    float *v_x, *v_y;
    float *p_x, *p_y;
    Color *colors;
    // fixed values
    float *mass;
    float *radius;
    // Metadata
    size_t length;
    size_t size;
};

struct Pair {
    int i, j;
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