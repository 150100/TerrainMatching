#ifndef NUMBERTYPE_H
#define NUMBERTYPE_H

#include <limits>

//foundation of doubles
static bool nearlyEqual(double x, double y)
{
    const double eps = std::numeric_limits<double>::epsilon();
    return (x/y < 1+2*eps) && (y/x < 1+2*eps);
}


#endif // NUMBERTYPE_H
