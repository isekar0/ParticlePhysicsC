#include "../includes/particles_solver.h"

#include "../includes/particles_arrays.h"
#include "../includes/structs.h"
#include "../includes/utils.h"
#include <assert.h>
#include <math.h>
#include <raylib.h>
#include <stddef.h>
#include <stdint.h>
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
#define NUM_INIT_SIZE 60
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

void particlesUpdate( struct Particles_t *particles, const float time_step ) {
    // float highest_energy_iteration = 0;
    for ( size_t idx = 0; idx < particles->length; idx++ ) {
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

void particlesCollisionsCheck_SweepPrune( struct Particles_t *particles,
                                          struct BoundEntry  *BE,
                                          const uint16_t      window_width,
                                          const uint16_t      window_height,
                                          float               time_step ) {
    size_t N = particles->length;

    // 1) Wall‐bounce / boundary checks (unchanged from before)
    for ( size_t idx = 0; idx < N; idx++ ) {
        float x = particles->position[ idx ].x;
        float y = particles->position[ idx ].y;
        float r = (float)particles->radius[ idx ];

        if ( x <= r ) {
            particles->position[ idx ].x = r;
            particles->velocity[ idx ].x *= -DAMPING_WALL;
        }
        if ( x >= window_width - r ) {
            particles->position[ idx ].x = window_width - r;
            particles->velocity[ idx ].x *= -DAMPING_WALL;
        }
        if ( y <= r ) {
            particles->position[ idx ].y = r;
            particles->velocity[ idx ].y *= -DAMPING_WALL;
        }
        if ( y >= window_height - r ) {
            particles->position[ idx ].y = window_height - r;
            particles->velocity[ idx ].y *= -DAMPING_WALL;
        }
    }

    // If there are 0 or 1 particles, there are no pairwise collisions:
    if ( N < 2 )
        return;

    // --- 2) Sweep-and-prune over sorted slots 0..N-1 ---
    for ( size_t a = 0; a < N - 1; a++ ) {
        // “a” is the pivot slot.  Retrieve its particle index i:
        int      i     = BE->idx[ a ];
        float    maxXi = BE->maxX[ a ]; // right edge of particle i’s x-interval
        float    px_i  = particles->position[ i ].x;
        float    py_i  = particles->position[ i ].y;
        float    ri    = (float)particles->radius[ i ];
        uint32_t m1    = particles->mass[ i ];

        // Compare against b = a+1.. until BE->minX[b] > maxXi
        for ( size_t b = a + 1; b < N; b++ ) {
            // If the next slot’s left edge exceeds pivot’s right edge, we’re done
            if ( BE->minX[ b ] > maxXi ) {
                break;
            }

            // Otherwise, perform a full 2D circle-circle test
            int   j    = BE->idx[ b ]; // the actual particle index at slot b
            float px_j = particles->position[ j ].x;
            float py_j = particles->position[ j ].y;
            float rj   = (float)particles->radius[ j ];

            // dx² + dy² <= (ri + rj)² ?
            float rsum = ri + rj;
            float dx   = px_i - px_j;

            // if ( dx > rsum || dx < -rsum ) {
            //     break;
            // }

            float dy = py_i - py_j;

            // if ( dy > rsum || dy < -rsum ) {
            //     break;
            // }

            float dist2 = dx * dx + dy * dy;
            if ( dist2 <= ( rsum * rsum ) ) {
                // They overlap—resolve via your existing impulse solver.
                uint32_t m2 = particles->mass[ j ];
                Vector2  v1 = particles->velocity[ i ];
                Vector2  v2 = particles->velocity[ j ];

                // Pass pointers so that positional correction actually updates particles->position[...]:
                struct Tuple_Vector2_t new_vs =
                    particlesResolveVelocity(
                        m1, m2,
                        v1, v2,
                        &particles->position[ i ],
                        &particles->position[ j ],
                        (uint32_t)ri, (uint32_t)rj,
                        time_step );

                particles->velocity[ i ] = new_vs.vec1;
                particles->velocity[ j ] = new_vs.vec2;
            }
        }
    }
}

struct Pair {
    int i, j;
};

typedef struct {
    int     i, j;                  // indices of the two colliding particles
    float   distance;              // Euclidean distance
    Vector2 normal;                // contact normal (unit)
    float   penetration_depth;     // raw depth
    float   penetration_corrected; // clamped depth
    float   velocity_correction;   // Baumgarte term, if any
    float   scalar_impulse;        // final impulse magnitude
    Vector2 impulse_vec;           // impulse vector = scalar_impulse * normal
} CollisionIntermediate;

void particlesResolveVelocity_NEW(
    struct Particles_t *particles,
    struct Pair        *collisions,
    size_t              collisionCount,
    float               time_step ) {

    size_t N = particles->length;
    if ( collisionCount == 0 || N < 2 )
        return;

    // 1) Allocate one big array of intermediates (heap):
    CollisionIntermediate *ci =
        (CollisionIntermediate *)malloc(
            collisionCount * sizeof( CollisionIntermediate ) );
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
        float dx           = particles->position[ i ].x - particles->position[ j ].x;
        float dy           = particles->position[ i ].y - particles->position[ j ].y;
        float dist2        = dx * dx + dy * dy;
        dist2              = fmaxf( dist2, EPSILON );
        float dist         = sqrtf( dist2 );
        ci[ idx ].distance = dist;

        if ( dist < EPSILON ) {
            ci[ idx ].normal.x = 1.0f;
            ci[ idx ].normal.y = 0.0f;
        } else {
            ci[ idx ].normal.x = dx / dist;
            ci[ idx ].normal.y = dy / dist;
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

        float    dist = ci[ idx ].distance;
        Vector2  n    = ci[ idx ].normal;
        uint32_t m1   = particles->mass[ i ];
        uint32_t m2   = particles->mass[ j ];

        // --- A) Positional correction (move out of penetration) ---
        if ( ci[ idx ].penetration_depth > 0.0f ) {
            float invM1    = 1.0f / (float)m1;
            float invM2    = 1.0f / (float)m2;
            float totalInv = invM1 + invM2; // = (m1+m2)/(m1*m2)
            float p1factor = ( invM1 / totalInv ) * ci[ idx ].penetration_corrected;
            float p2factor = ( invM2 / totalInv ) * ci[ idx ].penetration_corrected;

            particles->position[ i ].x += p1factor * n.x;
            particles->position[ i ].y += p1factor * n.y;
            particles->position[ j ].x -= p2factor * n.x;
            particles->position[ j ].y -= p2factor * n.y;
        }

        // --- B) Compute current relative velocity along normal ---
        Vector2 v1             = particles->velocity[ i ];
        Vector2 v2             = particles->velocity[ j ];
        Vector2 rel            = { v1.x - v2.x, v1.y - v2.y };
        float   velAlongNormal = rel.x * n.x + rel.y * n.y;

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
        ci[ idx ].impulse_vec.x  = jmag * n.x;
        ci[ idx ].impulse_vec.y  = jmag * n.y;

        // --- D) Apply velocity change to both particles ---
        particles->velocity[ i ].x += invM1 * ci[ idx ].impulse_vec.x;
        particles->velocity[ i ].y += invM1 * ci[ idx ].impulse_vec.y;
        particles->velocity[ j ].x -= invM2 * ci[ idx ].impulse_vec.x;
        particles->velocity[ j ].y -= invM2 * ci[ idx ].impulse_vec.y;
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
        float x = particles->position[ idx ].x;
        float y = particles->position[ idx ].y;
        float r = (float)particles->radius[ idx ];

        if ( x <= r ) {
            particles->position[ idx ].x = r;
            particles->velocity[ idx ].x *= -DAMPING_WALL;
        }
        if ( x >= window_width - r ) {
            particles->position[ idx ].x = window_width - r;
            particles->velocity[ idx ].x *= -DAMPING_WALL;
        }
        if ( y <= r ) {
            particles->position[ idx ].y = r;
            particles->velocity[ idx ].y *= -DAMPING_WALL;
        }
        if ( y >= window_height - r ) {
            particles->position[ idx ].y = window_height - r;
            particles->velocity[ idx ].y *= -DAMPING_WALL;
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
    int *cellHeads   = (int *)malloc( cellCount * sizeof( int ) );
    int *nextIndices = (int *)malloc( N * sizeof( int ) );
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
        float px = particles->position[ i ].x;
        float py = particles->position[ i ].y;

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
        float px = particles->position[ i ].x;
        float py = particles->position[ i ].y;
        float ri = (float)particles->radius[ i ];
        // uint32_t m1 = particles->mass[ i ];
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

                    float px_j = particles->position[ j ].x;
                    float py_j = particles->position[ j ].y;
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

    particlesResolveVelocity_NEW( particles, collisions, collisionCount, time_step );

// 7) Free our temporary arrays
cleanup:
    free( collisions );
    free( cellHeads );
    free( nextIndices );
}
