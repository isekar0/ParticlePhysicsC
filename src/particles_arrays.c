#include "../includes/particles_arrays.h"

#include "../includes/structs.h"
// #include "../includes/utils.h"
#include <assert.h>
#include <raylib.h>
#include <stdio.h>
#include <stdlib.h>

void boundsFree( struct BoundEntry *bounds ) {
    free( bounds->idx );
    free( bounds->maxX );
    free( bounds->minX );
}

// struct BoundEntry boundsCreate( size_t num_bounds ) {
//     struct BoundEntry bounds;

//     bounds.idx = calloc( num_bounds, sizeof( *bounds.idx ) );
//     if ( !bounds.idx ) {
//         fprintf( stderr, "Failed to initalize bounds" );
//         boundsFree( &bounds );
//         exit( EXIT_FAILURE );
//     }
//     bounds.minX = calloc( num_bounds, sizeof( *bounds.minX ) );
//     if ( !bounds.minX ) {
//         fprintf( stderr, "Failed to initalize bounds" );
//         boundsFree( &bounds );
//         exit( EXIT_FAILURE );
//     }
//     bounds.maxX = calloc( num_bounds, sizeof( *bounds.idx ) );
//     if ( !bounds.maxX ) {
//         fprintf( stderr, "Failed to initalize bounds" );
//         boundsFree( &bounds );
//         exit( EXIT_FAILURE );
//     }

//     bounds.size   = num_bounds;
//     bounds.length = 0;

//     return bounds;
// }

// void boundsExpand( struct BoundEntry *bounds ) {
//     size_t new_size = ( (int)bounds->size ) * ARRAY_GROWTH_FACTOR;
//     // printf("Expand array has been called, current size is %lu and new size will be %lu", bounds->arrays_size, new_size);
//     struct BoundEntry temp_bounds = *bounds;

//     temp_bounds.idx = realloc( temp_bounds.idx, new_size * sizeof( *bounds->idx ) );
//     if ( !temp_bounds.idx ) {
//         fprintf( stderr, "Failed to initialize bounds" );
//         boundsFree( &temp_bounds );
//         exit( EXIT_FAILURE );
//     }
//     temp_bounds.maxX = realloc( temp_bounds.maxX, new_size * sizeof( *bounds->maxX ) );
//     if ( !temp_bounds.maxX ) {
//         fprintf( stderr, "Failed to initialize bounds" );
//         boundsFree( &temp_bounds );
//         exit( EXIT_FAILURE );
//     }
//     temp_bounds.minX = realloc( temp_bounds.minX, new_size * sizeof( *bounds->minX ) );
//     if ( !temp_bounds.minX ) {
//         fprintf( stderr, "Failed to initialize bounds" );
//         boundsFree( &temp_bounds );
//         boundsFree( bounds ); // its joever

//         exit( EXIT_FAILURE );
//     }

//     // once all tests pass, set bounds to be the temp
//     // bounds->size = new_size;
//     bounds->idx  = temp_bounds.idx;
//     bounds->maxX = temp_bounds.maxX;
//     bounds->minX = temp_bounds.minX;
// }

struct BoundEntry *gBE_forSort = NULL;

/* Compare “row p” vs “row q” by their minX[p] vs minX[q] */
int CompareRows( size_t p, size_t q ) {
    float xp = gBE_forSort->minX[ p ];
    float xq = gBE_forSort->minX[ q ];
    if ( xp < xq )
        return -1;
    if ( xp > xq )
        return +1;
    return 0;
}

/* Swap all three fields at index p and q: idx[p]<->idx[q], minX[p]<->minX[q], maxX[p]<->maxX[q] */
void SwapRows( size_t p, size_t q ) {
    // swap idx[p], idx[q]
    int tmpIdx            = gBE_forSort->idx[ p ];
    gBE_forSort->idx[ p ] = gBE_forSort->idx[ q ];
    gBE_forSort->idx[ q ] = tmpIdx;

    // swap minX[p], minX[q]
    float tmpMin           = gBE_forSort->minX[ p ];
    gBE_forSort->minX[ p ] = gBE_forSort->minX[ q ];
    gBE_forSort->minX[ q ] = tmpMin;

    // swap maxX[p], maxX[q]
    float tmpMax           = gBE_forSort->maxX[ p ];
    gBE_forSort->maxX[ p ] = gBE_forSort->maxX[ q ];
    gBE_forSort->maxX[ q ] = tmpMax;
}

/*
  In‐place quicksort on the range [lo..hi], sorting by minX.
  We pick the middle element as pivot, then partition.
*/
void QuickSortSOA( size_t lo, size_t hi ) {
    if ( lo >= hi )
        return;

    size_t mid   = lo + ( hi - lo ) / 2;
    float  pivot = gBE_forSort->minX[ mid ];

    size_t i = lo;
    size_t j = hi;
    while ( i <= j ) {
        // Move i forward while minX[i] < pivot
        while ( gBE_forSort->minX[ i ] < pivot ) {
            i++;
        }
        // Move j backward while minX[j] > pivot
        while ( gBE_forSort->minX[ j ] > pivot ) {
            j--;
        }
        if ( i <= j ) {
            SwapRows( i, j );
            i++;
            if ( j == 0 )
                break; // avoid underflow
            j--;
        }
    }
    // Recurse on subranges
    if ( lo < j )
        QuickSortSOA( lo, j );
    if ( i < hi )
        QuickSortSOA( i, hi );
}

/*
  Public wrapper.  Call this on a freshly‐filled BoundEntry BE.
  After it returns, BE->minX is sorted ascending, and idx/maxX follow.
*/
void SortBoundEntrySOA( struct BoundEntry *BE ) {
    if ( BE == NULL || BE->length < 2 )
        return;
    gBE_forSort = BE;
    QuickSortSOA( 0, BE->length - 1 );
    gBE_forSort = NULL;
}

// void SortBoundEntrySOA( struct BoundEntry *BE );

void rebuildBoundEntry( struct Particles_t *p, struct BoundEntry *BE ) {
    size_t N = p->length;

    // 1) If the existing arrays are the wrong length, free + re‐allocate:
    if ( BE->length != N ) {
        // Free old arrays (if any):
        free( BE->idx );
        free( BE->minX );
        free( BE->maxX );

        // Allocate new arrays of length N:
        if ( N > 0 ) {
            BE->idx  = malloc( N * sizeof( int ) );
            BE->minX = malloc( N * sizeof( float ) );
            BE->maxX = malloc( N * sizeof( float ) );
            assert( BE->idx != NULL );
            assert( BE->minX != NULL );
            assert( BE->maxX != NULL );
        } else {
            BE->idx  = NULL;
            BE->minX = NULL;
            BE->maxX = NULL;
        }
        BE->length = N;
    }

    // 2) Fill the SoA in “unsorted” order (slot a → particle a):
    //    minX[a] = p->position[a].x − p->radius[a], etc.
    for ( size_t a = 0; a < N; a++ ) {
        BE->idx[ a ]  = (int)a;
        BE->minX[ a ] = p->position[ a ].x - (float)p->radius[ a ];
        BE->maxX[ a ] = p->position[ a ].x + (float)p->radius[ a ];
    }

    // 3) Sort the three parallel arrays so that minX[] is ascending.
    SortBoundEntrySOA( BE );
}

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
    particles.length = 0;
    particles.size   = num_particles_init;
    return particles;
}

void particlesExpandArrays( struct Particles_t *particles ) {
    size_t new_size = ( (int)particles->size ) * ARRAY_GROWTH_FACTOR;
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
    particles->size     = new_size;
    particles->mass     = temp_particles.mass;
    particles->radius   = temp_particles.radius;
    particles->colors   = temp_particles.colors;
    particles->position = temp_particles.position;
    particles->velocity = temp_particles.velocity;
}

void particlesAddParticle( struct Particles_t *particles, const Vector2 particle_i_position, const uint32_t particle_i_mass, const uint32_t particle_i_radius, uint32_t num_particles ) {

    // Expand size of arrays if they are full
    if ( particles->length + num_particles >= particles->size ) {
        particlesExpandArrays( particles );
    }

    // update number of particles stored
    while ( num_particles-- ) {
        particles->length += 1;

        // add new particle
        size_t index                   = particles->length;
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
    size_t index                   = particles->length;
    particles->mass[ index ]       = 0.0f;
    particles->radius[ index ]     = 0.0f;
    particles->colors[ index ]     = RAYWHITE;
    particles->position[ index ].x = 0.0f;
    particles->position[ index ].y = 0.0f;
    particles->velocity[ index ].x = 0.0f;
    particles->velocity[ index ].y = 0.0f;

    particles->length -= 1;
}

Texture2D CreateCircleTexture( int radius ) {
    int diameter = radius * 2;

    // 1) Create an Image in CPU memory, fully transparent
    Image circleImg = GenImageColor( diameter, diameter, (Color){ 0, 0, 0, 0 } );

    // 2) Draw a white circle (filled) into that image
    //    ImageDrawCircle(&img, centerX, centerY, radius, color);
    ImageDrawCircle( &circleImg, radius, radius, radius, WHITE );

    // 3) Upload that Image to the GPU as a Texture2D
    Texture2D circleTex = LoadTextureFromImage( circleImg );

    // 4) Free the CPU‐side Image; the GPU texture remains
    UnloadImage( circleImg );

    return circleTex;
}

static int       circleBaseRad = 12;
static Texture2D circleTex     = { 0 };

// void particlesDraw( struct Particles_t *particles ) {
//     for ( size_t idx = 0; idx < particles->length; idx++ ) {
//         DrawCircleV( particles->position[ idx ], (float)particles->radius[ idx ], particles->colors[ idx ] );
//     }
// }

void particlesDraw( struct Particles_t *particles, Texture2D circleTex ) {

    for ( size_t i = 0; i < particles->length; i++ ) {
        Vector2 pos  = particles->position[ i ];
        float   rad  = particles->radius[ i ];
        Color   tint = particles->colors[ i ];

        // 1) Compute scale factor so that a circle of radius `circleBaseRad`
        //    becomes radius `rad` on screen:
        float scale = rad / (float)circleBaseRad;

        // 2) After scaling, the (unscaled) texture is 2*circleBaseRad × 2*circleBaseRad.
        //    On screen it will become (2*circleBaseRad * scale)² = (2*rad)² in size.
        //    We want the center of the scaled texture to coincide with pos (the particle center).
        //    DrawTextureEx draws the top‐left corner at `drawPos`, so:
        float   scaledDiameter = ( (float)circleBaseRad * 2.0f ) * scale;
        Vector2 drawPos        = (Vector2){
            pos.x - ( scaledDiameter * 0.5f ),
            pos.y - ( scaledDiameter * 0.5f ) };

        // 3) Draw the circle texture at drawPos, with no rotation, at the computed scale,
        //    tinted by `tint`. The original texture was opaque‐white, so tinting changes its color.
        DrawTextureEx( circleTex, drawPos, 0.0f, scale, tint );
    }

    // Draw other UI, FPS, etc.
}