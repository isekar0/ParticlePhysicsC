#ifndef PARTICLES_ARRAYS_H
#define PARTICLES_ARRAYS_H

#include "raylib.h"
#include "structs.h"
#include <stdint.h>
#include <stdlib.h>

struct Particles_t;
struct BoundEntry;

void boundsFree( struct BoundEntry *bounds );
// struct BoundEntry  boundsCreate( size_t num_bounds );
// void               boundsExpand( struct BoundEntry *bounds );
void               rebuildBoundEntry( struct Particles_t *p, struct BoundEntry *BE );
void               particlesFree( struct Particles_t *particles );
struct Particles_t particlesCreate( size_t num_particles_init );
void               particlesAddParticle( struct Particles_t *particles, const Vector2 particle_i_position, const float particle_i_mass, const float particle_i_radius, uint8_t num_particles );
void               particlesRemoveParticle( struct Particles_t *particles );
Texture2D          CreateCircleTexture( int radius );
void               particlesDraw( struct Particles_t *particles, Texture2D );

#endif