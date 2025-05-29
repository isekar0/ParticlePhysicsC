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
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <time.h>
// #include <errno.h>
#include <sys/types.h>
#include "raylib.h"

#if defined(PLATFORM_DESKTOP)
    #define GLSL_VERSION            330
#else   // PLATFORM_ANDROID, PLATFORM_WEB
    #define GLSL_VERSION            100
#endif

#include <string.h> // for log files

// Metadata
// Keep the ratio as 16:9 and scale it up and down
#define WINDOW_SIZE_FACTOR 75 
#define WINDOW_TITLE "Particle simulation :3"

// Solver constants
#define SOLVER_TUNER 0.3f // larger values can cause this phenomenon where extra energy/momentum is introduced to the system, values greater than 0.4 result is massive impulse spikes
#define EPSILON 1e-8f // small value for division by zero correction
#define DEPTH_CAP 0.1f
#define SLOP_FACTOR 0.005f
#define BAUMGARTE 0.2f
#define GRAVITY 9.81 // Negated so it can be used correctly (positive makes particles go down)
#define AIR_RESISTENCE 0.98
#define RESTITUTION 0.90
#define DAMPING_WALL 0.90
#define TIMESTEP_FACTOR 5 // time_step /= TIMESTEP_FACTOR
#define ENERGY_MIN 1e-4

// Testing
#define NUM_INIT_SIZE 20
#define ARRAY_GROWTH_FACTOR 2

// Consider placing structures and functions relating to the physics of particles in a separate header
// Data structures
struct Particles_t {
    // Metadata
    size_t number_of_stored_particles;
    size_t arrays_size;
    // fixed values
    uint32_t *mass;
    uint32_t *radius;
    // varying values
    Color *colors;
    Vector2 *position; // raylib 2D vector = { .x, .y }
    Vector2 *velocity;
    float *energy;
    float highest_energy;
};

struct Tuple_Vector2_t {
    Vector2 vec1, vec2;
};

unsigned char clamp_color(float intensity); // 0 to 255
float clamp_f(float value, float min, float max);
Color TransformH(const Color in, const float fHue);
void particlesFree(struct Particles_t *particles); 
struct Particles_t particlesCreate(size_t num_particles_init); 
void particlesAddParticle(struct Particles_t *particles, const Vector2 particle_i_position, const uint32_t particle_i_mass, const uint32_t particle_i_radius);
void particlesRemoveParticle(struct Particles_t *particles);
void particlesCollisionsCheck(struct Particles_t *particles, const uint16_t window_width, const uint16_t window_height, float time_step);
void particlesUpdate(struct Particles_t *particles, const Vector2 acceleration, const float time_step);
void particlesDraw(struct Particles_t *particles);

int main(void) {
    // Debug
    FILE *p_file_last_log = fopen("log.txt", "w");
    if (!p_file_last_log) {
        fprintf(stderr, "Cannot open log file");
        exit(EXIT_FAILURE);
    }
    FILE *p_file_benchmarks = fopen("benchmarks.txt", "a+");
    if (!p_file_benchmarks) {
        fprintf(stderr, "Cannot open benchmarks file");
        exit(EXIT_FAILURE);
    }

    int test_counter = 0;
    int test_limit = 5000;
    // performance
    float running_time_per_particle = 0;
    float current_mean_time_per_particle = 0;

    // random seed
    srand(time(NULL));
    SetRandomSeed(time(NULL));

    // Initialize the raylib window 
    const uint16_t window_width = 16 * WINDOW_SIZE_FACTOR;
    const uint16_t window_height = 9 * WINDOW_SIZE_FACTOR;

    SetConfigFlags(FLAG_MSAA_4X_HINT);
    InitWindow(window_width, window_height, WINDOW_TITLE);

    unsigned int frames_counter = 0;

    SetTargetFPS(144);

    // Initialize 
    size_t num_particles_init = NUM_INIT_SIZE; // WILL CRASH IF NOT MULTIPLE OF 20
    struct Particles_t particles = particlesCreate(num_particles_init);
    struct Particles_t *p_particles = &particles;

    const Vector2 acceleration_field = { 0 , GRAVITY };

    //shaders
    Shader blur_effect = {0};
    blur_effect = LoadShader("shaders/base.vs", "shaders/blur.fs");

    int locRes  = GetShaderLocation(blur_effect, "resolution");
    int locDir  = GetShaderLocation(blur_effect, "direction");
    RenderTexture2D sceneRT = LoadRenderTexture(window_width, window_height);
    RenderTexture2D hpassRT = LoadRenderTexture(window_width, window_height);
    Vector2 res = { (float)window_width, (float)window_height };
    SetShaderValue(blur_effect, locRes, &res, SHADER_UNIFORM_VEC2);


    float last_time = GetFrameTime();
    // Main game loop
    while (!WindowShouldClose() && test_counter < test_limit) {
        // Update variables here
        
        // Check for left click to add particle, and right click to remove last particle added
        if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) {
            Vector2 mouse_position = GetMousePosition();
            uint32_t mass = GetRandomValue(5, 10);
            uint32_t radius = GetRandomValue(5, 15);
            particlesAddParticle(p_particles, mouse_position, mass, radius);
        }
        if (IsMouseButtonPressed(MOUSE_BUTTON_RIGHT) && p_particles->number_of_stored_particles > 0) {
            particlesRemoveParticle(p_particles);
        }

        // Debug
        Vector2 test_position_i = {(window_width / 2.0f) + GetRandomValue(-100, 100), (window_width / 2.0f) + GetRandomValue(-100, 100)};
        uint32_t m_i = GetRandomValue(5, 10);
        uint32_t r_i = GetRandomValue(5, 15);
        particlesAddParticle(p_particles, test_position_i, m_i, r_i);
        test_counter += 1;
        fprintf(p_file_last_log, " number of particles %d, FPS %d\n", test_counter, GetFPS());

        // e.g. position += velocity
        float current_time = GetFrameTime();
        float delta_time = current_time - last_time;
        running_time_per_particle += delta_time;
        delta_time *= TIMESTEP_FACTOR;
        // particlesUpdate(p_particles, acceleration_field, delta_time, window_width, window_height);
        
        // Check collisions
        if (p_particles->number_of_stored_particles > 0) {
            particlesCollisionsCheck(p_particles, window_width, window_height, delta_time);
            particlesUpdate(p_particles, acceleration_field, delta_time);

        }
        
        frames_counter++;

        // Draw
        // BeginDrawing(); {
            // ClearBackground(RAYWHITE);

            // particlesDraw(p_particles);

            BeginTextureMode(sceneRT);      // all draws now go into sceneRT
                ClearBackground(BLANK);     // or transparent black
                particlesDraw(p_particles); // uses DrawCircle(), no shader active
            EndTextureMode();

            BeginTextureMode(hpassRT);
                ClearBackground(BLANK);

                // tell shader “horizontal”
                float dirH[2] = { 1.0f, 0.0f };
                SetShaderValue(blur_effect, locDir, dirH, SHADER_UNIFORM_VEC2);

                BeginShaderMode(blur_effect);
                    // draw the off-screen colour buffer as a texture
                    DrawTextureRec(sceneRT.texture,
                                (Rectangle){0,0, (float)sceneRT.texture.width,
                                                    (float)-sceneRT.texture.height},   // -height flips Y
                                (Vector2){0,0}, WHITE);
                EndShaderMode();
            EndTextureMode();


            BeginDrawing();
                ClearBackground(RAYWHITE);

                float dirV[2] = { 0.0f, 1.0f };
                SetShaderValue(blur_effect, locDir, dirV, SHADER_UNIFORM_VEC2);

                BeginShaderMode(blur_effect);
                    DrawTextureRec(hpassRT.texture,
                                (Rectangle){0,0, (float)hpassRT.texture.width,
                                                    (float)-hpassRT.texture.height},
                                (Vector2){0,0}, WHITE);
                EndShaderMode();

                // DrawFPS(10, 10);
                DrawText(TextFormat("FPS: %d", GetFPS()), 10, 10, 20, LIGHTGRAY);
                DrawText(TextFormat("No. Particles: %lu", particles.number_of_stored_particles), window_width - 200, 10, 20, LIGHTGRAY);

            EndDrawing();

            
            

        // } 
        // EndDrawing();
    }

    // print out some benchmarks
    rewind(p_file_benchmarks);                /* go to beginning for reading */

    float last_mean_time_per_particle = 0.0f; /* initialise */

    /* Read file line by line, always keeping the last value we see */
    char line[1024];
    while (fgets(line, sizeof line, p_file_benchmarks)) {
        float tmp;
        if (sscanf(line,
                "Current mean time per particle = %f",  /* exact template */
                &tmp) == 1)
        {
            last_mean_time_per_particle = tmp;
        }
    }

    /* Now switch to append-mode writing */
    fseek(p_file_benchmarks, 0, SEEK_END);

    current_mean_time_per_particle = running_time_per_particle / particles.number_of_stored_particles;
    float percentage_difference = 2 * (last_mean_time_per_particle - current_mean_time_per_particle) / (last_mean_time_per_particle + current_mean_time_per_particle);
    fprintf(p_file_benchmarks,"Previous mean time per particle = %.6f\n", last_mean_time_per_particle);
    fprintf(p_file_benchmarks, "Percentage difference (+ve if improved) = %.2f%%\n", percentage_difference * 100);
    fprintf(p_file_benchmarks,"Current mean time per particle = %.6f\n",current_mean_time_per_particle);

    // free and shutdown
    UnloadRenderTexture(sceneRT);
    UnloadRenderTexture(hpassRT);
    UnloadShader(blur_effect);
    particlesFree(p_particles);
    CloseWindow();
    fclose(p_file_last_log);
    fclose(p_file_benchmarks);

    return 0;
}

unsigned char clamp(float intensity) {
    if (intensity < 0) {
        return 0;
    }
    if (intensity > 255) {
        return 255;
    }

    return (unsigned char) intensity;
}

float clamp_f(float val, float min, float max) {
    if (val < min) {
        return min;
    }
    if (val > max) {
        return max;
    }
    
    return val;
}

Color TransformH(const Color in, const float fHue) {
    Color out;
    const float cosA = cos(fHue*3.14159265f/180); //convert degrees to radians
    const float sinA = sin(fHue*3.14159265f/180); //convert degrees to radians
    const float sqrt_1_3 = sqrtf(1.0f/3.0f);
    const float one_minus_cosA = 1.0f - cosA;
    const float one_third = 1.0f / 3.0f;
    //calculate the rotation matrix, only depends on Hue
    float matrix[3][3] = {{cosA + one_minus_cosA * one_third, one_third* one_minus_cosA - sqrt_1_3 * sinA, one_third* one_minus_cosA + sqrt_1_3 * sinA},
        {one_third * one_minus_cosA + sqrt_1_3 * sinA, cosA + one_third* one_minus_cosA, one_third * one_minus_cosA - sqrt_1_3 * sinA},
        {one_third * one_minus_cosA - sqrt_1_3 * sinA, one_third * one_minus_cosA + sqrt_1_3 * sinA, cosA + one_third * one_minus_cosA}};
    //Use the rotation matrix to convert the RGB directly
    out.r = clamp(in.r*matrix[0][0] + in.g*matrix[0][1] + in.b*matrix[0][2]);
    out.g = clamp(in.r*matrix[1][0] + in.g*matrix[1][1] + in.b*matrix[1][2]);
    out.b = clamp(in.r*matrix[2][0] + in.g*matrix[2][1] + in.b*matrix[2][2]);
    return out;
}

void particlesFree(struct Particles_t *particles) {
    free(particles->mass);
    free(particles->radius);
    free(particles->colors);
    free(particles->position);
    free(particles->velocity);
    free(particles->energy);
    // free(particles);
}

struct Particles_t particlesCreate(size_t num_particles_init) {
    struct Particles_t particles;
    particles.mass = calloc(num_particles_init, sizeof(*particles.mass));
    if (!particles.mass) {
        fprintf(stderr, "Failed to initialize particles");
        particlesFree(&particles);
        exit(EXIT_FAILURE);
    }
    particles.radius = calloc(num_particles_init, sizeof(*particles.radius));
    if (!particles.radius) {
        fprintf(stderr, "Failed to initialize particles");
        particlesFree(&particles);
        exit(EXIT_FAILURE);
    }
    particles.colors = calloc(num_particles_init, sizeof(*particles.colors));
    if (!particles.colors) {
        fprintf(stderr, "Failed to initialize particles");
        particlesFree(&particles);
        exit(EXIT_FAILURE);
    }
    particles.position = calloc(num_particles_init, sizeof(*particles.position));
    if (!particles.position) {
        fprintf(stderr, "Failed to initialize particles");
        particlesFree(&particles);
        exit(EXIT_FAILURE);
    }
    particles.velocity = calloc(num_particles_init, sizeof(*particles.velocity));
    if (!particles.velocity) {
        fprintf(stderr, "Failed to initialize particles");
        particlesFree(&particles);
        exit(EXIT_FAILURE);
    }
    particles.energy = calloc(num_particles_init, sizeof(*particles.energy));
    if (!particles.energy) {
        fprintf(stderr, "Failed to initialize particles");
        particlesFree(&particles);
        exit(EXIT_FAILURE);
    }

    // Fill in metadata
    particles.number_of_stored_particles = 0;
    particles.arrays_size = num_particles_init;
    particles.highest_energy = 0;
    return particles;
} 

void particlesExpandArrays(struct Particles_t *particles) {
    size_t new_size = ((int) particles->arrays_size ) * ARRAY_GROWTH_FACTOR;
    // printf("Expand array has been called, current size is %lu and new size will be %lu", particles->arrays_size, new_size);
    struct Particles_t temp_particles = *particles; 

    temp_particles.mass = realloc(temp_particles.mass, new_size * sizeof(*particles->mass));
    if (!temp_particles.mass) {
        fprintf(stderr, "Failed to initialize particles");
        particlesFree(&temp_particles);
        exit(EXIT_FAILURE);
    }
    temp_particles.radius = realloc(temp_particles.radius, new_size * sizeof(*particles->radius));
    if (!temp_particles.radius) {
        fprintf(stderr, "Failed to initialize particles");
        particlesFree(&temp_particles);
        exit(EXIT_FAILURE);
    }
    temp_particles.colors = realloc(temp_particles.colors, new_size * sizeof(*particles->colors));
    if (!temp_particles.colors) {
        fprintf(stderr, "Failed to initialize particles");
        particlesFree(&temp_particles);
        exit(EXIT_FAILURE);
    }
    temp_particles.position = realloc(temp_particles.position, new_size * sizeof(*particles->position));
    if (!temp_particles.position) {
        fprintf(stderr, "Failed to initialize particles");
        particlesFree(&temp_particles);
        exit(EXIT_FAILURE);
    }
    temp_particles.velocity = realloc(temp_particles.velocity, new_size * sizeof(*particles->velocity));
    if (!temp_particles.velocity) {
        fprintf(stderr, "Failed to initialize particles");
        particlesFree(&temp_particles);
        exit(EXIT_FAILURE);
    }
     temp_particles.energy = realloc(temp_particles.energy, new_size * sizeof(*particles->energy));
    if (!temp_particles.energy) {
        fprintf(stderr, "Failed to initialize particles");
        particlesFree(&temp_particles);
        particlesFree(particles); // its joever
        exit(EXIT_FAILURE);
    }

    // once all tests pass, set particles to be the temp
    particles->arrays_size = new_size;
    particles->mass = temp_particles.mass;
    particles->radius = temp_particles.radius;
    particles->colors = temp_particles.colors;
    particles->position = temp_particles.position;
    particles->velocity = temp_particles.velocity;
    particles->energy = temp_particles.energy;
}

float energyCalc(uint32_t mass, Vector2 velocity) {
    return 0.5 * mass * (velocity.x * velocity.x + velocity.y * velocity.y);
}

void particlesAddParticle(struct Particles_t *particles, const Vector2 particle_i_position, const uint32_t particle_i_mass, const uint32_t particle_i_radius) {

    // Expand size of arrays if they are full 
    if (particles->number_of_stored_particles >= particles->arrays_size) {
        particlesExpandArrays(particles);
    }
    
    // update number of particles stored
    particles->number_of_stored_particles += 1;

    // add new particle
    size_t index = particles->number_of_stored_particles;
    particles->mass[index] = particle_i_mass;
    particles->radius[index] = particle_i_radius;
    particles->colors[index] = BLUE;
    particles->position[index].x = (float) particle_i_position.x;
    particles->position[index].y = (float) particle_i_position.y;
    particles->velocity[index].x = 0.0f;
    particles->velocity[index].y = 0.0f;
    particles->energy[index] = 0.0f;
}

void particlesRemoveParticle(struct Particles_t *particles) {
    size_t index = particles->number_of_stored_particles;
    particles->mass[index] = 0.0f;
    particles->radius[index] = 0.0f;
    particles->colors[index] = RAYWHITE;
    particles->position[index].x = 0.0f;
    particles->position[index].y = 0.0f;
    particles->velocity[index].x = 0.0f;
    particles->velocity[index].y = 0.0f;
    particles->energy[index] = 0.0f;

    particles->number_of_stored_particles -= 1;
}

// https://stackoverflow.com/a/35212639
// https://en.wikipedia.org/wiki/Elastic_collision
struct Tuple_Vector2_t particlesResolveVelocity(uint32_t m1, uint32_t m2, Vector2 v1, Vector2 v2, Vector2 p1, Vector2 p2, uint32_t r1, uint32_t r2, float time_step) {
    
    // Euclidean distance between the two particles
    float euclid_distance_squared = (p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y);
    euclid_distance_squared = fmaxf(euclid_distance_squared, EPSILON); // in case it is too close to zero
    float euclid_distance = sqrtf(euclid_distance_squared);
    // inverse of each mass
    float m1_inverse = 1.0f / m1;
    float m2_inverse = 1.0f / m2;

    // the contact normal between the two particles
    Vector2 contact_normal;
    if (euclid_distance < EPSILON) {
        euclid_distance = 1;
        // default values, change to be random values on the unit circle (random but still of length 1)
        
        contact_normal.x = 1.0f;
        contact_normal.y = 0.0f;
        // int random_angle = (float) rand() / (float) (RAND_MAX / (2 * PI));
        // contact_normal.x = 1.0f / cos(random_angle);
        // contact_normal.y = 1.0f / sin(random_angle);
    }
    else {
        contact_normal.x = (p1.x - p2.x) / euclid_distance;
        contact_normal.y = (p1.y - p2.y) / euclid_distance;
    }
    
    // positional correction if there is some penetration
    float velocity_correction_term = 0; // set to non zero if there is penetration
    float penetration_depth_corrected;
    float penetration_depth = r1 + r2 - euclid_distance;

    if (penetration_depth > 0.0f) {

        float penetration_depth_max = DEPTH_CAP * fminf(r1, r2);
        penetration_depth_corrected = penetration_depth * SOLVER_TUNER;
        penetration_depth_corrected = fmaxf(-penetration_depth_max, fminf(penetration_depth_corrected, penetration_depth_max));

        float p1_factor = (m1_inverse / (m1_inverse + m2_inverse)) * penetration_depth_corrected;
        float p2_factor = (m2_inverse / (m1_inverse + m2_inverse)) * penetration_depth_corrected;

        // correct positions
        p1.x += p1_factor * contact_normal.x;
        p1.y += p1_factor * contact_normal.y;
        p2.x -= p2_factor * contact_normal.x;
        p2.y -= p2_factor * contact_normal.y;

        // recompute
        euclid_distance = sqrtf(fmaxf((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y), EPSILON));
        contact_normal.x = (p1.x - p2.x) / euclid_distance;
        contact_normal.y = (p1.y - p2.y) / euclid_distance;

        float residual_depth = penetration_depth - penetration_depth_corrected;

        // set velocity_correction_term to a non zero value if the remaining depth fulfills the condition
        float slop = SLOP_FACTOR * fmaxf(r1, r2);
        if (residual_depth > slop) {
            velocity_correction_term = BAUMGARTE * (residual_depth - slop) / time_step;
        }
    }
    
    // Component of relative velocity between the particles will also be the same
    Vector2 relative_velocity = { v1.x - v2.x, v1.y - v2.y };
    float velocity_normal = relative_velocity.x * contact_normal.x + relative_velocity.y * contact_normal.y;

    if (velocity_normal > 0.0f && velocity_correction_term == 0.0f) {
        struct Tuple_Vector2_t old_velocities = { v1, v2 };
        return old_velocities;
    }

    float velocity_normal_desired = -RESTITUTION * velocity_normal + velocity_correction_term;

    float scalar_impulse = (velocity_normal_desired - velocity_normal) / (m1_inverse + m2_inverse);

    Vector2 vector_impulse = { scalar_impulse * contact_normal.x, scalar_impulse * contact_normal.y };

    Vector2 v1_prime = {v1.x + vector_impulse.x * m1_inverse, v1.y + vector_impulse.y * m1_inverse };
    Vector2 v2_prime = { v2.x - vector_impulse.x * m2_inverse, v2.y - vector_impulse.y * m2_inverse };
    
    struct Tuple_Vector2_t new_velocities = { v1_prime, v2_prime };

    return new_velocities;
}

// Maybe integrate the raylib function for partilce-particle collisions
// bool CheckCollisionCircles(Vector2 center1, float radius1, Vector2 center2, float radius2);  
void particlesCollisionsCheck(struct Particles_t *particles, const uint16_t window_width, const uint16_t window_height, float time_step) {
    
    for (size_t idx = 0; idx < particles->number_of_stored_particles; idx++) {

        if (particles->position[idx].x <= 0) {
            particles->position[idx].x = 1;
            particles->velocity[idx].x *= -1 * DAMPING_WALL;
        }

        if (particles->position[idx].x >= window_width) {
            particles->position[idx].x = window_width - 1;
            particles->velocity[idx].x *= -1 * DAMPING_WALL;
        }

        if (particles->position[idx].y < 0) {
            particles->position[idx].y = 1;
            particles->velocity[idx].y *= -1 * DAMPING_WALL;
            // particles->velocity[idx].y += particles->mass[idx] * GRAVITY;
        }

        if (particles->position[idx].y >= window_height) {
            particles->position[idx].y = window_height - 1;
            particles->velocity[idx].y *= -1 * DAMPING_WALL;
        }
    }

    for (size_t idx = 0; idx < particles->number_of_stored_particles - 1; idx++) {
        Vector2 position_idx = particles->position[idx];
        uint32_t radius_idx = particles->radius[idx];
        uint32_t mass_idx = particles->mass[idx];
        Vector2 velocity_idx = particles->velocity[idx];
        for (size_t jdx = idx + 1; jdx < particles->number_of_stored_particles; jdx++) {
            if (CheckCollisionCircles(position_idx, radius_idx, particles->position[jdx], particles->radius[jdx])) {
                Vector2 position_jdx = particles->position[jdx];
                uint32_t mass_jdx = particles->mass[jdx];
                uint32_t radius_jdx = particles->radius[jdx];
                Vector2 velocity_jdx = particles->velocity[jdx];
                // solve for final velocities
                struct Tuple_Vector2_t v1_v2_new = particlesResolveVelocity(mass_idx, mass_jdx, velocity_idx, velocity_jdx, position_idx, position_jdx, radius_idx, radius_jdx, time_step);
                particles->velocity[idx] = v1_v2_new.vec1;
                particles->velocity[jdx] = v1_v2_new.vec2;
            }
        }
    }
}    

void particlesUpdate(struct Particles_t *particles, const Vector2 acceleration, const float time_step) {
    float highest_energy_iteration = 0;
    for (size_t idx = 0; idx < particles->number_of_stored_particles; idx++) {
        // update the positions
        particles->position[idx].x += 0.5 * acceleration.x * time_step * time_step + particles->velocity[idx].x * time_step;
        particles->position[idx].y += 0.5 * acceleration.y * time_step * time_step + particles->velocity[idx].y * time_step;
    
        // update the velocities
        particles->velocity[idx].x += acceleration.x * time_step * AIR_RESISTENCE; 
        particles->velocity[idx].y += acceleration.y * time_step * AIR_RESISTENCE;

        // update the energy
        particles->energy[idx] = energyCalc(particles->mass[idx], particles->velocity[idx]);

        float energy_idx = particles->energy[idx];
        if (energy_idx > highest_energy_iteration) {
            highest_energy_iteration = energy_idx;
        }
        
    }

    particles->highest_energy = fmaxf(particles->highest_energy * 0.95, highest_energy_iteration);
    particles->highest_energy = fmaxf(particles->highest_energy, 50);

    // update the colors
    for (size_t idx = 0; idx < particles->number_of_stored_particles; idx++) {
        float temp = (particles->energy[idx] + ENERGY_MIN) / (particles->highest_energy - ENERGY_MIN); 
        temp = clamp_f(temp, 0, 1);
        float nth_root = 5.0f;
        temp = pow(temp, 1.0f / nth_root);

        // const float TAU = 0.05f;
        // float temp_shoulder = (temp - TAU) / (1 - TAU);
        
        // float temp_biased = clamp_f(temp_shoulder, 0.0f, 1.0f); 
        
        particles->colors[idx] = TransformH(RED, 360 * (1 - temp) - 60);
        particles->colors[idx].a = 255;
    }
    
}

void particlesDraw(struct Particles_t *particles) {
    for (size_t idx = 0; idx < particles->number_of_stored_particles; idx++) {
        // uint32_t m_i = particles->position[idx];
        // DrawCircleV(particles->position[idx], (float) particles->radius[idx], particles->colors[idx]);
        Color color_outer = TransformH(particles->colors[idx], 10);
        DrawCircleGradient(particles->position[idx].x, particles->position[idx].y, particles->radius[idx] * 2, particles->colors[idx], color_outer);
    }
}