#ifndef NUMBERTYPE_H
#define NUMBERTYPE_H

#include <limits>

//foundation of doubles
static bool nearlyEqual(double x, double y)
{
	const double eps = 0.0000000001;
	return (y - eps <= x) && (x <= y + eps);
}
static bool nearlyInClosedRange(double x, double range_min, double range_max)
{
	const double eps = 0.0000000001;
	return (range_min - eps <= x) && (x <= range_max + eps);
}
static bool strictlyInOpenRange(double x, double range_min, double range_max)
{
	const double eps = 0.0000000001;
	return (range_min + eps < x) && (x < range_max - eps);
}

#endif // NUMBERTYPE_H
