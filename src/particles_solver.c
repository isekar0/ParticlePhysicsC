#include "../includes/particles_solver.h"

#include "../includes/constants.h"
#include "../includes/particles_arrays.h"
#include "../includes/structs.h"
#include "../includes/utils.h"
#include <assert.h>
#include <math.h>
#include <raylib.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

// https://stackoverflow.com/a/35212639
// https://en.wikipedia.org/wiki/Elastic_collision
void particlesUpdate( struct Particles_t *particles, const float time_step ) {
    // float highest_energy_iteration = 0;
    for ( size_t idx = 0; idx < particles->length; idx++ ) {
        // update the positions
        particles->p_x[ idx ] += 0.5 * 0 * time_step * time_step + particles->v_x[ idx ] * time_step;
        particles->p_y[ idx ] += 0.5 * GRAVITY * time_step * time_step + particles->v_y[ idx ] * time_step;

        // update the velocities
        particles->v_x[ idx ] *= ( 1 - AIR_RESISTENCE * time_step );
        particles->v_y[ idx ] *= ( 1 - AIR_RESISTENCE * time_step );

        particles->v_x[ idx ] += GRAVITY * time_step;
        particles->v_y[ idx ] += GRAVITY * time_step;

        // update the color
        float speed_idx = particles->v_x[ idx ] * particles->v_x[ idx ] + particles->v_y[ idx ] * particles->v_y[ idx ];
        float temp      = ( particles->mass[ idx ] * speed_idx + ENERGY_MIN ) / ( 1000000 - ENERGY_MIN );

        temp           = clamp_f( temp, 0, 1 );
        float nth_root = 3.0f;
        temp           = pow( temp, 1.0f / nth_root );

        particles->colors[ idx ]   = TransformH( RED, 240 * ( 1 - temp ) );
        particles->colors[ idx ].a = 255;

        // Check Highest values
    }
}

typedef struct {
    int   i, j;     // indices of the two colliding particles
    float distance; // Euclidean distance
    float normal_x; // contact normal (unit)
    float normal_y;
    float penetration_depth;     // raw depth
    float penetration_corrected; // clamped depth
    float velocity_correction;   // Baumgarte term, if any
    float scalar_impulse;        // final impulse magnitude
    float impulse_vec_x;         // impulse vector = scalar_impulse * normal
    float impulse_vec_y;
} CollisionIntermediate;

void particlesResolveVelocity( struct Particles_t *particles, struct Pair *collisions, size_t collisionCount, float time_step ) {
    size_t N = particles->length;
    if ( collisionCount == 0 || N < 2 )
        return;

    // 1) Allocate one big array of intermediates (heap):
    CollisionIntermediate *ci =
        (CollisionIntermediate *)aligned_alloc(
            32, collisionCount * sizeof( CollisionIntermediate ) );
    if ( !ci )
        return; // out of memory

    // 2) Fill in geometry + Baumgarte data for each collision:
    for ( size_t idx = 0; idx < collisionCount; idx++ ) {
        int i = collisions[ idx ].i;
        int j = collisions[ idx ].j;
        // Defensive bounds check:
        if ( i < 0 || i >= (int)N || j < 0 || j >= (int)N ) {
            // Skip invalid pair; mark a dummy entry. We'll never use it in pass (3).
            ci[ idx ].i = -1;
            ci[ idx ].j = -1;
            continue;
        }
        ci[ idx ].i = collisions[ idx ].i;
        ci[ idx ].j = collisions[ idx ].j;

        // --- A) Compute dx, dy, dist, normal ---
        float dx           = particles->p_x[ i ] - particles->p_x[ j ];
        float dy           = particles->p_y[ i ] - particles->p_y[ j ];
        float dist2        = dx * dx + dy * dy;
        dist2              = fmaxf( dist2, EPSILON );
        float dist         = sqrtf( dist2 );
        ci[ idx ].distance = dist;

        if ( dist < EPSILON ) {
            ci[ idx ].normal_x = 1.0f;
            ci[ idx ].normal_y = 0.0f;
        } else {
            ci[ idx ].normal_x = dx / dist;
            ci[ idx ].normal_y = dy / dist;
        }

        // --- B) Compute penetration depths + correction ---
        float ri                    = (float)particles->radius[ i ];
        float rj                    = (float)particles->radius[ j ];
        float rawDepth              = ( ri + rj ) - dist;
        ci[ idx ].penetration_depth = rawDepth;

        if ( rawDepth > 0.0f ) {
            float maxDepth                  = DEPTH_CAP * fminf( ri, rj );
            float corrected                 = rawDepth * SOLVER_TUNER;
            corrected                       = fmaxf( -maxDepth, fminf( corrected, maxDepth ) );
            ci[ idx ].penetration_corrected = corrected;
        } else {
            ci[ idx ].penetration_corrected = 0.0f;
        }

        // --- C) Compute Baumgarte / “velocity correction” term ---
        if ( rawDepth > 0.0f ) {
            float residual = rawDepth - ci[ idx ].penetration_corrected;
            float slopVal  = SLOP_FACTOR * fmaxf( ri, rj );
            if ( residual > slopVal ) {
                ci[ idx ].velocity_correction =
                    BAUMGARTE * ( residual - slopVal ) / time_step;
            } else {
                ci[ idx ].velocity_correction = 0.0f;
            }
        } else {
            ci[ idx ].velocity_correction = 0.0f;
        }
    }

    // 3) Now do the actual impulse/positional correction pass:
    for ( size_t idx = 0; idx < collisionCount; idx++ ) {
        // If we flagged i = -1 because of an out‐of‐bounds pair, skip it:
        int i = ci[ idx ].i;
        int j = ci[ idx ].j;
        if ( i < 0 || j < 0 )
            continue;

        float   dist = ci[ idx ].distance;
        float   n_x  = ci[ idx ].normal_x;
        float   n_y  = ci[ idx ].normal_y;
        uint8_t m1   = particles->mass[ i ];
        uint8_t m2   = particles->mass[ j ];

        // --- A) Positional correction (move out of penetration) ---
        if ( ci[ idx ].penetration_depth > 0.0f ) {
            float invM1    = 1.0f / (float)m1;
            float invM2    = 1.0f / (float)m2;
            float totalInv = invM1 + invM2; // = (m1+m2)/(m1*m2)
            float p1factor = ( invM1 / totalInv ) * ci[ idx ].penetration_corrected;
            float p2factor = ( invM2 / totalInv ) * ci[ idx ].penetration_corrected;

            particles->p_x[ i ] += p1factor * n_x;
            particles->p_y[ i ] += p1factor * n_y;
            particles->p_x[ j ] -= p2factor * n_x;
            particles->p_y[ j ] -= p2factor * n_y;
        }

        // --- B) Compute current relative velocity along normal ---
        float v1_x  = particles->v_x[ i ];
        float v1_y  = particles->v_y[ i ];
        float v2_x  = particles->v_x[ j ];
        float v2_y  = particles->v_y[ j ];
        float rel_x = v1_x - v2_x;
        float rel_y = v1_y - v2_y;

        float velAlongNormal = rel_x * n_x + rel_y * n_y;

        // If they are separating (velAlongNormal > 0) and no Baumgarte push,
        // skip the impulse.
        if ( velAlongNormal > 0.0f && ci[ idx ].velocity_correction == 0.0f ) {
            continue;
        }

        // --- C) Compute desired normal velocity and impulse scalar ---
        float vNormDesired = -RESTITUTION * velAlongNormal + ci[ idx ].velocity_correction;

        float invM1 = 1.0f / (float)m1;
        float invM2 = 1.0f / (float)m2;
        float denom = invM1 + invM2; // = (m1 + m2)/(m1*m2)
        float jmag  = ( vNormDesired - velAlongNormal ) / denom;

        ci[ idx ].scalar_impulse = jmag;
        ci[ idx ].impulse_vec_x  = jmag * n_x;
        ci[ idx ].impulse_vec_y  = jmag * n_y;

        // --- D) Apply velocity change to both particles ---
        particles->v_x[ i ] += invM1 * ci[ idx ].impulse_vec_x;
        particles->v_y[ i ] += invM1 * ci[ idx ].impulse_vec_y;
        particles->v_x[ j ] -= invM2 * ci[ idx ].impulse_vec_x;
        particles->v_y[ j ] -= invM2 * ci[ idx ].impulse_vec_y;
    }

    // 4) Clean up:
    free( ci );
}

void particlesCollisionsCheck_Grid(
    struct Particles_t *particles,
    const uint16_t      window_width,
    const uint16_t      window_height,
    float               time_step ) {

    size_t N = particles->length;
    if ( N < 2 )
        return;

    for ( size_t idx = 0; idx < N; idx++ ) {
        float x = particles->p_x[ idx ];
        float y = particles->p_y[ idx ];
        float r = (float)particles->radius[ idx ];

        if ( x <= r ) {
            particles->p_x[ idx ] = r;
            particles->v_x[ idx ] *= -DAMPING_WALL;
        }
        if ( x >= window_width - r ) {
            particles->p_x[ idx ] = window_width - r;
            particles->v_x[ idx ] *= -DAMPING_WALL;
        }
        if ( y <= r ) {
            particles->p_y[ idx ] = r;
            particles->v_y[ idx ] *= -DAMPING_WALL;
        }
        if ( y >= window_height - r ) {
            particles->p_y[ idx ] = window_height - r;
            particles->v_y[ idx ] *= -DAMPING_WALL;
        }
    }

    // 1) Find the maximum particle radius (O(N))
    float maxRadius = 6.0f;
    for ( size_t i = 0; i < N; i++ ) {
        float r = (float)particles->radius[ i ];
        if ( r > maxRadius )
            maxRadius = r;
    }
    // If all radii are zero (unlikely), bail out
    if ( maxRadius <= 0.0f )
        return;

    // 2) Choose cellSize = 2 * maxRadius.
    //    That guarantees any two colliding particles must lie
    //    in the same or neighbor cell.
    float cellSize = 2.0f * maxRadius;

    // 3) Compute number of cells in x and y directions
    int cols      = (int)ceilf( (float)window_width / cellSize );
    int rows      = (int)ceilf( (float)window_height / cellSize );
    int cellCount = cols * rows;

    // 4) Allocate & initialize the "head of list" per cell.
    //    We'll use an array of length `cellCount`, each entry is
    //    the index of the first particle in that bucket (or -1 if empty).
    int *cellHeads   = (int *)aligned_alloc( 32, cellCount * sizeof( int ) );
    int *nextIndices = (int *)aligned_alloc( 32, N * sizeof( int ) );
    if ( !cellHeads || !nextIndices ) {
        // Allocation failure: clean up and return
        if ( cellHeads )
            free( cellHeads );
        if ( nextIndices )
            free( nextIndices );
        return;
    }
    // Initialize all heads to -1, and all nextIndices to -1
    for ( int c = 0; c < cellCount; c++ )
        cellHeads[ c ] = -1;
    for ( size_t i = 0; i < N; i++ )
        nextIndices[ i ] = -1;

    // 5) Insert each particle into its grid cell:
    //    Compute (cx, cy) = floor(x/cellSize), floor(y/cellSize).
    //    Cell index = cy*cols + cx. Then push into that bucket.
    for ( int i = 0; i < (int)N; i++ ) {
        float px = particles->p_x[ i ];
        float py = particles->p_y[ i ];

        // Clamp into [0..width) or [0..height). If out of bounds,
        // you can either skip or clamp to valid cell. We assume the
        // boundary‐bounce code has already kept positions valid.
        int cx = (int)( px / cellSize );
        int cy = (int)( py / cellSize );
        if ( cx < 0 )
            cx = 0;
        if ( cy < 0 )
            cy = 0;
        if ( cx >= cols )
            cx = cols - 1;
        if ( cy >= rows )
            cy = rows - 1;

        int bucketIndex = cy * cols + cx;
        // Chain‐in at the front
        nextIndices[ i ]         = cellHeads[ bucketIndex ];
        cellHeads[ bucketIndex ] = i;
    }

    // 6) For each particle i, check only the 3×3 neighborhood
    //    (its own cell plus eight neighbors). We only compare
    //    to j > i to avoid double‐testing pairs.

    size_t capacity = N * 4; // heuristic: expect “O(N)” collisions in many fluid cases
    if ( capacity < 16 )
        capacity = 16; // never go below 16
    struct Pair *collisions = (struct Pair *)malloc( capacity * sizeof( struct Pair ) );
    if ( !collisions ) {
        free( cellHeads );
        free( nextIndices );
        return;
    }
    size_t collisionCount = 0;
    for ( int i = 0; i < (int)N; i++ ) {
        float px = particles->p_x[ i ];
        float py = particles->p_y[ i ];
        float ri = (float)particles->radius[ i ];
        // uint8_t m1 = particles->mass[ i ];
        // Vector2  v1 = particles->velocity[ i ];

        // Recompute the cell of i (slightly redundant but simpler)
        int cx = (int)( px / cellSize );
        int cy = (int)( py / cellSize );
        if ( cx < 0 )
            cx = 0;
        if ( cy < 0 )
            cy = 0;
        if ( cx >= cols )
            cx = cols - 1;
        if ( cy >= rows )
            cy = rows - 1;

        // Loop over neighbor cells dx = -1,0,+1; dy = -1,0,+1
        for ( int dy = -1; dy <= +1; dy++ ) {
            int ny = cy + dy;
            if ( ny < 0 || ny >= rows )
                continue;
            for ( int dx = -1; dx <= +1; dx++ ) {
                int nx = cx + dx;
                if ( nx < 0 || nx >= cols )
                    continue;

                int neighborBucket = ny * cols + nx;
                // Traverse the linked list in neighborBucket
                for ( int j = cellHeads[ neighborBucket ]; j != -1; j = nextIndices[ j ] ) {
                    // Only test j > i to avoid double checks and skip self
                    if ( j <= i )
                        continue;

                    float px_j = particles->p_x[ j ];
                    float py_j = particles->p_y[ j ];
                    float rj   = (float)particles->radius[ j ];

                    // Quick AABB‐style reject (optional, but x/y separation check is cheap):
                    float dx_ij = px - px_j;
                    if ( dx_ij > ( ri + rj ) || dx_ij < -( ri + rj ) )
                        continue;
                    float dy_ij = py - py_j;
                    if ( dy_ij > ( ri + rj ) || dy_ij < -( ri + rj ) )
                        continue;

                    // Precise circle‐circle test:
                    float dist2 = dx_ij * dx_ij + dy_ij * dy_ij;
                    float rsum2 = ( ri + rj ) * ( ri + rj );
                    if ( dist2 <= rsum2 ) {

                        if ( collisionCount >= capacity ) {
                            // need more room → double the capacity
                            size_t       newCap = capacity * 2;
                            struct Pair *tmp    = (struct Pair *)realloc( collisions, newCap * sizeof( struct Pair ) );
                            if ( !tmp ) {
                                // failed to expand; bail early
                                goto cleanup;
                            }
                            collisions = tmp;
                            capacity   = newCap;
                        }

                        collisions[ collisionCount ].i = i;
                        collisions[ collisionCount ].j = j;
                        collisionCount++;
                    }
                }
            }
        }
    }

    // particlesResolveVelocity_NEW( particles, collisions, collisionCount, time_step );
    particlesResolveVelocity( particles, collisions, collisionCount, time_step );

// 7) Free our temporary arrays
cleanup:
    free( collisions );
    free( cellHeads );
    free( nextIndices );
}
