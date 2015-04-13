//
//  plane.h
//  3D minimum extent finder
//
//  Created by 윤상덕 on 2014. 6. 28..
//  Copyright (c) 2014년 GA. All rights reserved.
//

#ifndef PLANE_H
#define PLANE_H


#include <iostream>
//#include <gmp.h>


//점 class
class Point
{
public:
    double x,y,z;
    Point();
    Point(double a, double b, double c);
    //Point(Number_type a, Number_type b, Number_type c);
};


class Line
{
public:
    //y = ax + b
    double a,b;
    Line();
    Line(double x, double y);

    bool operator==(Line &l) {return l.a == a && l.b == b;}

    void print() {std::cout << a << "x + " << b << " = 0" << std::endl;}
};


class Plane
{
public:
    double a, b, c, d; // ax + by + cz = d

    Plane();
    Plane(double a_, double b_, double c_, double d_);

    bool operator==(Plane &p) {return p.a == a && p.b == b && p.c == c && p.d == d;}

    // indicator=1 when a_123 is from patch, indicator=-1 when a_123 is from domain.
    void setPlane_VTpair(Point &a_1, Point &a_2, Point &a_3, Point &b_1, int indicator); // a_123:triangle, b_1:vertex
    void setPlane_EEpair(Point &a_1, Point &a_2, Point &b_1, Point &b_2, int indicator); // a_12:edge, b_12:edge

    //void print() {std::cout << a.to_double() << "x + " << b.to_double() << "y + " << c.to_double() << "z = " << d.to_double() << std::endl;}
    void print() {std::cout << a << "x + " << b << "y + " << c << "z = " << d << std::endl;}
};


#endif // PLANE_H
