#include "plane.h"


Point::Point()
{
    x = y = z = 0.0;
}

Point::Point(double a, double b, double c)
{
    x = a;
    y = b;
    z = c;
}

//Point::Point(Number_type a, Number_type b, Number_type c)
//{
//    x = a.to_double();
//    y = b.to_double();
//    z = c.to_double();
//}


Line::Line()
{
    a = 0;
    b = 0;
}

Line::Line(double x, double y)
{
    a = x;
    b = y;
}


Plane::Plane()
{
    a=b=c=d=0;
}

Plane::Plane(double a_, double b_, double c_, double d_)
{
    a = a_; b=b_; c=c_; d=d_;
}

void Plane::setPlane_VTpair(Point &a_1, Point &a_2, Point &a_3, Point &b_1, int indicator)
{
    //a_2 - a_1 vector
    Point vector_1(a_2.x-a_1.x, a_2.y-a_1.y, a_2.z-a_1.z);

    //a_3 - a_1 vector
    Point vector_2(a_3.x-a_1.x, a_3.y-a_1.y, a_3.z-a_1.z);

    //normal vectors of the plane
    double normal_x, normal_y, normal_z;

    //cross product
    normal_x = vector_1.y * vector_2.z - vector_1.z * vector_2.y;
    normal_y = -vector_1.x * vector_2.z + vector_1.z * vector_2.x;
    normal_z = vector_1.x * vector_2.y - vector_1.y * vector_2.x;

    // Plane equation : normal_x(x-ref_x) + normal_y(y-ref_y) + normal_z(z-vertical_dist) = 0
    // (b_1-a_1)*indicator is translation vector that makes zero-distance.
    a = normal_x;
    b = normal_y;
    c = normal_z;
    d = indicator * (normal_x * (b_1.x - a_1.x) + normal_y * (b_1.y - a_1.y) + normal_z * (b_1.z - a_1.z));
}

void Plane::setPlane_EEpair(Point &a_1, Point &a_2, Point &b_1, Point &b_2, int indicator)
{
    //a_2 - a_1 vector
    Point vector_1(a_2.x-a_1.x, a_2.y-a_1.y, a_2.z-a_1.z);

    //b_2 - b_1 vector
    Point vector_2(b_2.x-b_1.x, b_2.y-b_1.y, b_2.z-b_1.z);

    //평면의 normal vector
    double normal_x, normal_y, normal_z;

    //cross product
    normal_x = vector_1.y * vector_2.z - vector_1.z * vector_2.y;
    normal_y = -vector_1.x * vector_2.z + vector_1.z * vector_2.x;
    normal_z = vector_1.x * vector_2.y - vector_1.y * vector_2.x;

    // Plane equation : normal_x(x-tr.x) + normal_y(y-tr.y) + normal_z(z-tr.z) = 0
    // where tr is a translation vector that makes zero-distance.
    // (b_1-a_1)*indicator is translation vector that makes zero-distance.
    a = normal_x;
    b = normal_y;
    c = normal_z;
    d = indicator * (normal_x * (b_1.x - a_1.x) + normal_y * (b_1.y - a_1.y) + normal_z * (b_1.z - a_1.z));
}
