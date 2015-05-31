#ifndef NUMBERTYPE_H
#define NUMBERTYPE_H

#include <limits>

//foundation of doubles
static bool nearlyEqual(double x, double y)
{
	const double eps = 0.0000000000001;
	double ratio = x / y;
    return (ratio < 1+2*eps) && (ratio > 1-2*eps);
}
static bool nearlyInRange(double x, double range_min, double range_max)
{
	const double eps = 0.0000000000001;
	double range_epsmin = std::min(range_min * (1 - eps), range_min * (1 + eps));
	double range_epsmax = std::max(range_max * (1 - eps), range_max * (1 + eps));
	return (range_epsmin <= x) && (x <= range_epsmax);
}

#endif // NUMBERTYPE_H
