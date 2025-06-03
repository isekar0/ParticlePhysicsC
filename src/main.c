/*
 * =====================================================================================
 *
 *       Filename:  main.c
 *
 *    Description: Basic fluid/particle simulation in C using Raylib for rendering
 *
 *        Version:  1.1
 *        Created:  05/11/2025 03:59:13 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Martin Hanna
 *   Organization:  N/A
 *
 * =====================================================================================
 */

// Includes
#include "raylib.h"
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <time.h>

#include "../includes/constants.h"
#include "../includes/particles_arrays.h"
#include "../includes/particles_solver.h"
#include "../includes/structs.h"

int main( void ) {
    // Debug
    int      test_counter       = 0;
    int      test_limit         = 4500;
    uint32_t test_num_particles = 10;

    // random seed
    srand( time( NULL ) );
    SetRandomSeed( time( NULL ) );

    // Initialize the raylib window
    const uint16_t window_width  = 16 * WINDOW_SIZE_FACTOR;
    const uint16_t window_height = 9 * WINDOW_SIZE_FACTOR;

    SetConfigFlags( FLAG_MSAA_4X_HINT );
    InitWindow( window_width, window_height, WINDOW_TITLE );

    unsigned int frames_counter = 0;

    SetTargetFPS( 144 );

    // Initialize
    size_t              num_particles_init = NUM_INIT_SIZE; // WILL CRASH IF NOT MULTIPLE OF 20
    struct Particles_t  particles          = particlesCreate( num_particles_init );
    struct Particles_t *p_particles        = &particles;

    static Texture2D circleTex = { 0 };
    circleTex                  = CreateCircleTexture( 6 );

    float previous_time = GetFrameTime();
    // Main game loop
    while ( !WindowShouldClose() && test_counter < test_limit ) {
        // Update variables here

        // Check for left click to add particle, and right click to remove last particle added
        if ( IsMouseButtonPressed( MOUSE_BUTTON_LEFT ) ) {
            Vector2  mouse_position = GetMousePosition();
            uint32_t mass           = GetRandomValue( 5, 10 );
            uint32_t radius         = GetRandomValue( 5, 15 );
            particlesAddParticle( p_particles, mouse_position, mass, radius, 1 );
        }
        if ( IsMouseButtonPressed( MOUSE_BUTTON_RIGHT ) && p_particles->length > 0 ) {
            particlesRemoveParticle( p_particles );
        }

        // Debug
        if ( test_counter < test_limit ) {
            Vector2  test_position_i = { GetRandomValue( 0, window_width ), GetRandomValue( 0, window_height ) };
            uint32_t m_i             = GetRandomValue( 5, 10 );
            uint32_t r_i             = GetRandomValue( 2, 8 );
            particlesAddParticle( p_particles, test_position_i, m_i, r_i, test_num_particles );
        }
        // e.g. position += velocity
        float current_time = GetFrameTime();
        float delta_time   = current_time - previous_time;

        delta_time = fminf( delta_time, 1.0f / 60 );
        // Check collisions
        if ( p_particles->length > 0 ) {
            particlesUpdate( p_particles, delta_time );
            particlesCollisionsCheck_Grid(
                p_particles,
                window_width, window_height,
                delta_time );
        }

        frames_counter++;

        BeginDrawing();
        ClearBackground( RAYWHITE );

        particlesDraw( p_particles, circleTex );

        // DrawFPS(10, 10);
        DrawText( TextFormat( "FPS: %d", GetFPS() ), 10, 10, 20, BLACK );
        DrawText( TextFormat( "No. Particles: %lu", particles.length ), window_width - 200, 10, 20, BLACK );

        EndDrawing();

        test_counter += 1;
    }

    // free and shutdown
    particlesFree( p_particles );
    CloseWindow();

    return 0;
}