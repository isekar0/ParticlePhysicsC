#include "../includes/particles_arrays.h"

#include "../includes/structs.h"
// #include "../includes/utils.h"
#include <stdio.h>

void particlesFree( struct Particles_t *particles ) {
    free( particles->mass );
    free( particles->radius );
    free( particles->colors );
    free( particles->position );
    free( particles->velocity );
    // free(particles);
}

struct Particles_t particlesCreate( size_t num_particles_init ) {
    struct Particles_t particles;
    particles.mass = calloc( num_particles_init, sizeof( *particles.mass ) );
    if ( !particles.mass ) {
        fprintf( stderr, "Failed to initialize particles" );
        particlesFree( &particles );
        exit( EXIT_FAILURE );
    }
    particles.radius = calloc( num_particles_init, sizeof( *particles.radius ) );
    if ( !particles.radius ) {
        fprintf( stderr, "Failed to initialize particles" );
        particlesFree( &particles );
        exit( EXIT_FAILURE );
    }
    particles.colors = calloc( num_particles_init, sizeof( *particles.colors ) );
    if ( !particles.colors ) {
        fprintf( stderr, "Failed to initialize particles" );
        particlesFree( &particles );
        exit( EXIT_FAILURE );
    }
    particles.position = calloc( num_particles_init, sizeof( *particles.position ) );
    if ( !particles.position ) {
        fprintf( stderr, "Failed to initialize particles" );
        particlesFree( &particles );
        exit( EXIT_FAILURE );
    }
    particles.velocity = calloc( num_particles_init, sizeof( *particles.velocity ) );
    if ( !particles.velocity ) {
        fprintf( stderr, "Failed to initialize particles" );
        particlesFree( &particles );
        exit( EXIT_FAILURE );
    }

    // Fill in metadata
    particles.number_of_stored_particles = 0;
    particles.arrays_size                = num_particles_init;
    return particles;
}

void particlesExpandArrays( struct Particles_t *particles ) {
    size_t new_size = ( (int)particles->arrays_size ) * ARRAY_GROWTH_FACTOR;
    // printf("Expand array has been called, current size is %lu and new size will be %lu", particles->arrays_size, new_size);
    struct Particles_t temp_particles = *particles;

    temp_particles.mass = realloc( temp_particles.mass, new_size * sizeof( *particles->mass ) );
    if ( !temp_particles.mass ) {
        fprintf( stderr, "Failed to initialize particles" );
        particlesFree( &temp_particles );
        exit( EXIT_FAILURE );
    }
    temp_particles.radius = realloc( temp_particles.radius, new_size * sizeof( *particles->radius ) );
    if ( !temp_particles.radius ) {
        fprintf( stderr, "Failed to initialize particles" );
        particlesFree( &temp_particles );
        exit( EXIT_FAILURE );
    }
    temp_particles.colors = realloc( temp_particles.colors, new_size * sizeof( *particles->colors ) );
    if ( !temp_particles.colors ) {
        fprintf( stderr, "Failed to initialize particles" );
        particlesFree( &temp_particles );
        exit( EXIT_FAILURE );
    }
    temp_particles.position = realloc( temp_particles.position, new_size * sizeof( *particles->position ) );
    if ( !temp_particles.position ) {
        fprintf( stderr, "Failed to initialize particles" );
        particlesFree( &temp_particles );
        exit( EXIT_FAILURE );
    }
    temp_particles.velocity = realloc( temp_particles.velocity, new_size * sizeof( *particles->velocity ) );
    if ( !temp_particles.velocity ) {
        fprintf( stderr, "Failed to initialize particles" );
        particlesFree( &temp_particles );
        particlesFree( particles ); // its joever

        exit( EXIT_FAILURE );
    }

    // once all tests pass, set particles to be the temp
    particles->arrays_size = new_size;
    particles->mass        = temp_particles.mass;
    particles->radius      = temp_particles.radius;
    particles->colors      = temp_particles.colors;
    particles->position    = temp_particles.position;
    particles->velocity    = temp_particles.velocity;
}

void particlesAddParticle( struct Particles_t *particles, const Vector2 particle_i_position, const uint32_t particle_i_mass, const uint32_t particle_i_radius, uint32_t num_particles ) {

    // Expand size of arrays if they are full
    if ( particles->number_of_stored_particles + num_particles >= particles->arrays_size ) {
        particlesExpandArrays( particles );
    }

    // update number of particles stored
    while ( num_particles-- ) {
        particles->number_of_stored_particles += 1;

        // add new particle
        size_t index                   = particles->number_of_stored_particles;
        particles->mass[ index ]       = particle_i_mass;
        particles->radius[ index ]     = particle_i_radius;
        particles->colors[ index ]     = BLUE;
        particles->position[ index ].x = (float)particle_i_position.x + particle_i_radius * ( num_particles - 1 );
        particles->position[ index ].y = (float)particle_i_position.y + particle_i_radius * ( num_particles - 1 );
        particles->velocity[ index ].x = 0.0f;
        particles->velocity[ index ].y = 0.0f;
    }
}

void particlesRemoveParticle( struct Particles_t *particles ) {
    size_t index                   = particles->number_of_stored_particles;
    particles->mass[ index ]       = 0.0f;
    particles->radius[ index ]     = 0.0f;
    particles->colors[ index ]     = RAYWHITE;
    particles->position[ index ].x = 0.0f;
    particles->position[ index ].y = 0.0f;
    particles->velocity[ index ].x = 0.0f;
    particles->velocity[ index ].y = 0.0f;

    particles->number_of_stored_particles -= 1;
}

void particlesDraw( struct Particles_t *particles ) {
    for ( size_t idx = 0; idx < particles->number_of_stored_particles; idx++ ) {
        // uint32_t m_i = particles->position[idx];
        DrawCircleV( particles->position[ idx ], (float)particles->radius[ idx ], particles->colors[ idx ] );
        // Color color_outer = TransformH( particles->colors[ idx ], 0 );
        // DrawCircleGradient( particles->position[ idx ].x, particles->position[ idx ].y, particles->radius[ idx ], particles->colors[ idx ], color_outer );
    }
}