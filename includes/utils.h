#ifndef UTILS_H
#define UTILS_H

#include <raylib.h>

unsigned char clamp_color( float intensity ); // 0 to 255
float         clamp_f( float value, float min, float max );
Color         TransformH( const Color in, const float fHue );

#endif