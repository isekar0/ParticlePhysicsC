#include "../includes/particles_solver.h"

#include "../includes/particles_arrays.h"
#include "../includes/structs.h"
#include "../includes/utils.h"
#include <assert.h>
#include <math.h>
#include <raylib.h>
#include <stdio.h>

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
#define NUM_INIT_SIZE 20
#define ARRAY_GROWTH_FACTOR 2

// https://stackoverflow.com/a/35212639
// https://en.wikipedia.org/wiki/Elastic_collision
// OPTIMIZE
struct Tuple_Vector2_t particlesResolveVelocity( uint32_t m1, uint32_t m2, Vector2 v1, Vector2 v2, Vector2 *p1, Vector2 *p2, uint32_t r1, uint32_t r2, float time_step ) {

    // Euclidean distance between the two particles
    float euclid_distance_squared = ( p1->x - p2->x ) * ( p1->x - p2->x ) + ( p1->y - p2->y ) * ( p1->y - p2->y );
    euclid_distance_squared       = fmaxf( euclid_distance_squared, EPSILON ); // in case it is too close to zero
    float euclid_distance         = sqrtf( euclid_distance_squared );
    // inverse of each mass
    float m1_inverse = 1.0f / m1;
    float m2_inverse = 1.0f / m2;

    // the contact normal between the two particles
    Vector2 contact_normal;
    if ( euclid_distance < EPSILON ) {
        euclid_distance  = 1;
        contact_normal.x = 1.0f;
        contact_normal.y = 0.0f;
    } else {
        contact_normal.x = ( p1->x - p2->x ) / euclid_distance;
        contact_normal.y = ( p1->y - p2->y ) / euclid_distance;
    }

    // positional correction if there is some penetration
    float velocity_correction_term = 0; // set to non zero if there is penetration
    float penetration_depth_corrected;
    float penetration_depth = r1 + r2 - euclid_distance;

    if ( penetration_depth > 0.0f ) {

        float penetration_depth_max = DEPTH_CAP * fminf( r1, r2 );
        penetration_depth_corrected = penetration_depth * SOLVER_TUNER;
        penetration_depth_corrected = fmaxf( -penetration_depth_max, fminf( penetration_depth_corrected, penetration_depth_max ) );

        float p1_factor = ( m1_inverse / ( m1_inverse + m2_inverse ) ) * penetration_depth_corrected;
        float p2_factor = ( m2_inverse / ( m1_inverse + m2_inverse ) ) * penetration_depth_corrected;

        // correct positions
        p1->x += p1_factor * contact_normal.x;
        p1->y += p1_factor * contact_normal.y;
        p2->x -= p2_factor * contact_normal.x;
        p2->y -= p2_factor * contact_normal.y;

        // recompute
        euclid_distance  = sqrtf( fmaxf( ( p1->x - p2->x ) * ( p1->x - p2->x ) + ( p1->y - p2->y ) * ( p1->y - p2->y ), EPSILON ) );
        contact_normal.x = ( p1->x - p2->x ) / euclid_distance;
        contact_normal.y = ( p1->y - p2->y ) / euclid_distance;

        float residual_depth = penetration_depth - penetration_depth_corrected;

        // set velocity_correction_term to a non zero value if the remaining depth fulfills the condition
        float slop = SLOP_FACTOR * fmaxf( r1, r2 );
        if ( residual_depth > slop ) {
            velocity_correction_term = BAUMGARTE * ( residual_depth - slop ) / time_step;
            // velocity_correction_term = ( BAUMGARTE / ( m1 + m2 ) ) * residual_depth;
        }
    }

    // Component of relative velocity between the particles will also be the same
    Vector2 relative_velocity = { v1.x - v2.x, v1.y - v2.y };
    float   velocity_normal   = relative_velocity.x * contact_normal.x + relative_velocity.y * contact_normal.y;

    if ( velocity_normal > 0.0f && velocity_correction_term == 0.0f ) {
        struct Tuple_Vector2_t old_velocities = { v1, v2 };
        return old_velocities;
    }

    float velocity_normal_desired = -RESTITUTION * velocity_normal + velocity_correction_term;

    float scalar_impulse = ( velocity_normal_desired - velocity_normal ) / ( m1_inverse + m2_inverse );

    Vector2 vector_impulse = { scalar_impulse * contact_normal.x, scalar_impulse * contact_normal.y };

    Vector2 v1_prime = { v1.x + vector_impulse.x * m1_inverse, v1.y + vector_impulse.y * m1_inverse };
    Vector2 v2_prime = { v2.x - vector_impulse.x * m2_inverse, v2.y - vector_impulse.y * m2_inverse };

    struct Tuple_Vector2_t new_velocities = { v1_prime, v2_prime };

    return new_velocities;
}

// Maybe integrate the raylib function for partilce-particle collisions
// bool CheckCollisionCircles(Vector2 center1, float radius1, Vector2 center2, float radius2);
void particlesCollisionsCheck( struct Particles_t *particles, const uint16_t window_width, const uint16_t window_height, float time_step ) {

    for ( size_t idx = 0; idx < particles->number_of_stored_particles; idx++ ) {

        if ( particles->position[ idx ].x <= 0 + particles->radius[ idx ] ) {
            particles->position[ idx ].x = particles->radius[ idx ];
            particles->velocity[ idx ].x *= -1 * DAMPING_WALL;
        }

        if ( particles->position[ idx ].x >= window_width - particles->radius[ idx ] ) {
            particles->position[ idx ].x = window_width - particles->radius[ idx ];
            particles->velocity[ idx ].x *= -1 * DAMPING_WALL;
        }

        if ( particles->position[ idx ].y < 0 + particles->radius[ idx ] ) {
            particles->position[ idx ].y = particles->radius[ idx ];
            particles->velocity[ idx ].y *= -1 * DAMPING_WALL;
            // particles->velocity[idx].y += particles->mass[idx] * GRAVITY;
        }

        if ( particles->position[ idx ].y >= window_height - particles->radius[ idx ] ) {
            particles->position[ idx ].y = window_height - particles->radius[ idx ];
            particles->velocity[ idx ].y *= -1 * DAMPING_WALL;
        }
    }

    for ( size_t idx = 0; idx < particles->number_of_stored_particles - 1; idx++ ) {
        Vector2  position_idx = particles->position[ idx ];
        uint32_t radius_idx   = particles->radius[ idx ];
        uint32_t mass_idx     = particles->mass[ idx ];
        Vector2  velocity_idx = particles->velocity[ idx ];
        for ( size_t jdx = idx + 1; jdx < particles->number_of_stored_particles; jdx++ ) {
            if ( CheckCollisionCircles( position_idx, radius_idx, particles->position[ jdx ], particles->radius[ jdx ] ) ) {
                Vector2  position_jdx = particles->position[ jdx ];
                uint32_t mass_jdx     = particles->mass[ jdx ];
                uint32_t radius_jdx   = particles->radius[ jdx ];
                Vector2  velocity_jdx = particles->velocity[ jdx ];
                // solve for final velocities
                struct Tuple_Vector2_t v1_v2_new = particlesResolveVelocity( mass_idx, mass_jdx, velocity_idx, velocity_jdx, &position_idx, &position_jdx, radius_idx, radius_jdx, time_step );
                particles->velocity[ idx ]       = v1_v2_new.vec1;
                particles->velocity[ jdx ]       = v1_v2_new.vec2;
            }
        }
    }
}

void particlesUpdate( struct Particles_t *particles, const float time_step ) {
    // float highest_energy_iteration = 0;
    for ( size_t idx = 0; idx < particles->number_of_stored_particles; idx++ ) {
        // update the positions
        particles->position[ idx ].x += 0.5 * 0 * time_step * time_step + particles->velocity[ idx ].x * time_step;
        particles->position[ idx ].y += 0.5 * GRAVITY * time_step * time_step + particles->velocity[ idx ].y * time_step;

        // update the velocities
        // particles->velocity[ idx ].x += 0 * time_step * AIR_RESISTENCE;
        // particles->velocity[ idx ].y += GRAVITY * time_step * AIR_RESISTENCE;
        float speed_idx = particles->velocity[ idx ].x * particles->velocity[ idx ].x + particles->velocity[ idx ].y * particles->velocity[ idx ].y;
        // float sqrt_speed_idx = sqrtf( speed_idx );

        // Vector2 drag_velocity = { -AIR_RESISTENCE * particles->velocity[ idx ].x * sqrt_speed_idx, -AIR_RESISTENCE * particles->velocity[ idx ].y * sqrt_speed_idx };

        particles->velocity[ idx ].x *= ( 1 - AIR_RESISTENCE * time_step );
        particles->velocity[ idx ].y *= ( 1 - AIR_RESISTENCE * time_step );

        particles->velocity[ idx ].x += GRAVITY * time_step;
        particles->velocity[ idx ].y += GRAVITY * time_step;

        // update the color
        float temp = ( particles->mass[ idx ] * speed_idx + ENERGY_MIN ) / ( 1000000 - ENERGY_MIN );

        temp           = clamp_f( temp, 0, 1 );
        float nth_root = 3.0f;
        temp           = pow( temp, 1.0f / nth_root );

        particles->colors[ idx ]   = TransformH( RED, 240 * ( 1 - temp ) );
        particles->colors[ idx ].a = 255;
    }
}

// A helper struct for sorting by minX
typedef struct {
    int   idx;  // original particle index
    float minX; // position[idx].x - radius[idx]
    float maxX; // position[idx].x + radius[idx]
} BoundEntry;

// Compare two BoundEntry by their minX
static int CompareBoundEntry( const void *a, const void *b ) {
    const BoundEntry *A = (const BoundEntry *)a;
    const BoundEntry *B = (const BoundEntry *)b;
    if ( A->minX < B->minX )
        return -1;
    if ( A->minX > B->minX )
        return +1;
    return 0;
}

void particlesCollisionsCheck_SweepPrune( struct Particles_t *particles,
                                          const uint16_t      window_width,
                                          const uint16_t      window_height,
                                          float               time_step ) {
    size_t N = particles->number_of_stored_particles;

    // 1) Wall‐bounce / boundary checks (unchanged from before)
    for ( size_t idx = 0; idx < N; idx++ ) {
        float x = particles->position[ idx ].x;
        float y = particles->position[ idx ].y;
        float r = (float)particles->radius[ idx ];

        if ( x <= r ) {
            particles->position[ idx ].x = r;
            particles->velocity[ idx ].x *= -1.0f * DAMPING_WALL;
        }
        if ( x >= (float)window_width - r ) {
            particles->position[ idx ].x = (float)window_width - r;
            particles->velocity[ idx ].x *= -1.0f * DAMPING_WALL;
        }
        if ( y <= r ) {
            particles->position[ idx ].y = r;
            particles->velocity[ idx ].y *= -1.0f * DAMPING_WALL;
        }
        if ( y >= (float)window_height - r ) {
            particles->position[ idx ].y = (float)window_height - r;
            particles->velocity[ idx ].y *= -1.0f * DAMPING_WALL;
        }
    }

    // If there are 0 or 1 particles, there are no pairwise collisions:
    if ( N < 2 )
        return;

    // 2) Build and sort the “bounds” array by minX
    BoundEntry *bounds = malloc( N * sizeof( BoundEntry ) );
    if ( !bounds )
        return; // out of memory—just bail.
    for ( size_t i = 0; i < N; i++ ) {
        float cx         = particles->position[ i ].x;
        float r          = (float)particles->radius[ i ];
        bounds[ i ].idx  = (int)i;
        bounds[ i ].minX = cx - r;
        bounds[ i ].maxX = cx + r;
    }
    qsort( bounds, N, sizeof( BoundEntry ), CompareBoundEntry );

    // 3) “Sweep” through the sorted list
    // For each entry a = 0..N-2, compare against b = a+1.. until minX[b] > maxX[a].
    for ( size_t a = 0; a < N - 1; a++ ) {
        int   i     = bounds[ a ].idx;
        float maxXi = bounds[ a ].maxX;

        // HOTSPOT
        // Look at all subsequent entries whose minX ≤ maxXi
        for ( size_t b = a + 1; b < N; b++ ) {
            if ( bounds[ b ].minX > maxXi ) {
                // The sorted‐by‐minX property guarantees no further x‐overlaps
                break;
            }

            int j = bounds[ b ].idx;
            // At this point, the intervals [minX[i], maxX[i]] and [minX[j], maxX[j]] do overlap in x.
            // We now do the full circle‐circle test (y‐coordinate + radii).
            Vector2 pi = particles->position[ i ];
            Vector2 pj = particles->position[ j ];
            float   ri = (float)particles->radius[ i ];
            float   rj = (float)particles->radius[ j ];

            if ( CheckCollisionCircles( pi, ri, pj, rj ) ) {
                // **CALL your impulse solver** exactly as before.
                uint32_t m1 = particles->mass[ i ];
                uint32_t m2 = particles->mass[ j ];
                Vector2  v1 = particles->velocity[ i ];
                Vector2  v2 = particles->velocity[ j ];

                // Here, pass &particles->position[i] and &particles->position[j]
                // so that positional correction is actually written back:
                struct Tuple_Vector2_t new_vs =
                    particlesResolveVelocity( m1, m2, v1, v2,
                                              &particles->position[ i ],
                                              &particles->position[ j ],
                                              (uint32_t)ri, (uint32_t)rj,
                                              time_step );

                particles->velocity[ i ] = new_vs.vec1;
                particles->velocity[ j ] = new_vs.vec2;
            }
        }
    }

    free( bounds );
}
