#include "../includes/utils.h"

#include <math.h>

#include "raylib.h"

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