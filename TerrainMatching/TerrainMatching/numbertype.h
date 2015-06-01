#ifndef NUMBERTYPE_H
#define NUMBERTYPE_H

#include <limits>

//foundation of doubles
static bool nearlyEqual(double x, double y)
{
	const double eps = 0.000000000001;
	std::pair<double, double> range_eps = std::minmax(y * (1 + eps), y * (1 - eps));
	return (range_eps.first <= x) && (x <= range_eps.second);
}
static bool nearlyInClosedRange(double x, double range_min, double range_max)
{
	const double eps = 0.000000000001;
	double range_epsmin = std::min(range_min * (1 - eps), range_min * (1 + eps));
	double range_epsmax = std::max(range_max * (1 - eps), range_max * (1 + eps));
	return (range_epsmin <= x) && (x <= range_epsmax);
}
static bool strictlyInOpenRange(double x, double range_min, double range_max)
{
	const double eps = 0.000000000001;
	double range_epsmin = std::max(range_min * (1 - eps), range_min * (1 + eps));
	double range_epsmax = std::min(range_max * (1 - eps), range_max * (1 + eps));
	return (range_epsmin < x) && (x < range_epsmax);
}

#endif // NUMBERTYPE_H
