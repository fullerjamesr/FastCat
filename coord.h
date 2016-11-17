#pragma once
#ifndef FASTCAT_COORD_H
#define FASTCAT_COORD_H

#include <math.h>
#include <stdbool.h>

static const double coord_golden_ratio = 1.618033988749895;
static const double coord_pi = 3.14159265358979323846;

typedef struct coord
{
    double x;
    double y;
    double z;
} coord;

static const coord UNITIALIZED_COORD = {0.0, 0.0, 0.0};

static inline double distance_squared(const coord c1, const coord c2)
{
    return ((c1.x - c2.x)*(c1.x - c2.x) + (c1.y - c2.y)*(c1.y - c2.y) + (c1.z - c2.z)*(c1.z - c2.z));
}

static inline double distance(const coord c1, const coord c2)
{
    return sqrt(distance_squared(c1,c2));
}

static inline double dot_product(const coord a, const coord b)
{
    return (a.x*b.x + a.y*b.y + a.z*b.z);
}

static inline coord vec_add(const coord c1, const coord c2)
{
    coord ret = {c1.x + c2.x, c1.y + c2.y, c1.z + c2.z};
    return ret;
}

static inline coord vec_subtract(const coord c1, const coord c2)
{
    coord ret = {c1.x - c2.x, c1.y - c2.y, c1.z - c2.z};
    return ret;
}

static inline coord vec_multiply(const coord c, double m)
{
    coord ret = {c.x * m, c.y * m, c.z * m};
    return ret;
}

static inline bool closeish(const coord c1, const coord c2, const double cutoff)
{
    coord dist = vec_subtract(c1, c2);
    return ((fabs(dist.x) < cutoff) && (fabs(dist.y) < cutoff) && (fabs(dist.z) < cutoff));
}

static inline double sphere_volume(double r)
{
    return 4.0 / 3.0 * coord_pi * r * r * r;
}

size_t golden_spiral_distribution(coord * vectors, size_t num_vectors);

#endif //FASTCAT_COORD_H