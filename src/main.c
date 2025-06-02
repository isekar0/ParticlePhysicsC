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
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
// #include <errno.h>
#include "raylib.h"
#include <sys/types.h>

#include "../includes/particles_arrays.h"
#include "../includes/particles_solver.h"
#include "../includes/structs.h"
// #include "../includes/utils.h"

#if defined( PLATFORM_DESKTOP )
#define GLSL_VERSION 330
#else // PLATFORM_ANDROID, PLATFORM_WEB
#define GLSL_VERSION 100
#endif

// Metadata
// Keep the ratio as 16:9 and scale it up and down
#define WINDOW_SIZE_FACTOR 75
#define WINDOW_TITLE "Particle simulation :3"

int main( void ) {
    // Debug
    int      test_counter       = 0;
    int      test_limit         = 1000;
    uint32_t test_num_particles = 40;

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

    // struct BoundEntry  bounds   = boundsCreate( num_particles_init );
    // struct BoundEntry *p_bounds = &bounds;

    // struct BoundEntry bounds = { NULL, NULL, NULL, 0 };

    static Texture2D circleTex = { 0 };
    circleTex                  = CreateCircleTexture( 6 );

    // const Vector2 acceleration_field = { 0 , GRAVITY };

    // shaders
    //  Shader blur_effect = {0};
    //  blur_effect = LoadShader("shaders/base.vs", "shaders/blur.fs");

    // int locRes  = GetShaderLocation(blur_effect, "resolution");
    // int locDir  = GetShaderLocation(blur_effect, "direction");
    // RenderTexture2D sceneRT = LoadRenderTexture(window_width, window_height);
    // RenderTexture2D hpassRT = LoadRenderTexture(window_width, window_height);
    // Vector2 res = { (float)window_width, (float)window_height };
    // SetShaderValue(blur_effect, locRes, &res, SHADER_UNIFORM_VEC2);

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
            // Vector2 test_position_i = { 500, 500 };
            Vector2  test_position_i = { ( window_width / 2.0f ) + GetRandomValue( -100, 100 ), ( window_width / 2.0f ) + GetRandomValue( -100, 100 ) };
            uint32_t m_i             = GetRandomValue( 5, 10 );
            uint32_t r_i             = GetRandomValue( 3, 6 );
            particlesAddParticle( p_particles, test_position_i, m_i, r_i, test_num_particles );
        }
        // e.g. position += velocity
        float current_time = GetFrameTime();
        float delta_time   = current_time - previous_time;
        // float delta_time = 1.0f / GetFPS();

        // delta_time *= 10;
        // previous_time = current_time;

        delta_time = fminf( delta_time, 1.0f / 60 );
        // Check collisions
        if ( p_particles->length > 0 ) {
            // particlesCollisionsCheck( p_particles, window_width, window_height, delta_time );
            particlesUpdate( p_particles, delta_time );
            // rebuildBoundEntry( p_particles, &bounds );
            // particlesCollisionsCheck_SweepPrune( p_particles, &bounds, window_width, window_height, delta_time );
            particlesCollisionsCheck_Grid(
                p_particles,
                window_width, window_height,
                delta_time );
        }

        frames_counter++;

        // Draw
        // BeginTextureMode(sceneRT);      // all draws now go into sceneRT
        //     ClearBackground(BLANK);     // or transparent black
        //     particlesDraw(p_particles); // uses DrawCircle(), no shader active
        // EndTextureMode();

        // BeginTextureMode(hpassRT);
        //     ClearBackground(BLANK);

        //     // tell shader “horizontal”
        //     float dirH[2] = { 1.0f, 0.0f };
        //     SetShaderValue(blur_effect, locDir, dirH, SHADER_UNIFORM_VEC2);

        //     BeginShaderMode(blur_effect);
        //         // draw the off-screen colour buffer as a texture
        //         DrawTextureRec(sceneRT.texture,
        //                     (Rectangle){0,0, (float)sceneRT.texture.width,
        //                                         (float)-sceneRT.texture.height},   // -height flips Y
        //                     (Vector2){0,0}, WHITE);
        //     EndShaderMode();
        // EndTextureMode();

        BeginDrawing();
        ClearBackground( RAYWHITE );

        // float dirV[2] = { 0.0f, 1.0f };
        // SetShaderValue(blur_effect, locDir, dirV, SHADER_UNIFORM_VEC2);

        // BeginShaderMode(blur_effect);
        //     DrawTextureRec(hpassRT.texture,
        //                 (Rectangle){0,0, (float)hpassRT.texture.width,
        //                                     (float)-hpassRT.texture.height},
        //                 (Vector2){0,0}, WHITE);
        // EndShaderMode();

        particlesDraw( p_particles, circleTex );

        // DrawFPS(10, 10);
        DrawText( TextFormat( "FPS: %d", GetFPS() ), 10, 10, 20, BLACK );
        DrawText( TextFormat( "No. Particles: %lu", particles.length ), window_width - 200, 10, 20, BLACK );

        EndDrawing();

        test_counter += 1;
    }

    // free and shutdown
    // UnloadRenderTexture(sceneRT);
    // UnloadRenderTexture(hpassRT);
    // UnloadShader(blur_effect);
    // boundsFree( &bounds );
    particlesFree( p_particles );
    CloseWindow();

    return 0;
}