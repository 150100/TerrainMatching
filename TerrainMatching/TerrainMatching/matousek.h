//
//  matousek.h
//  3D minimum extent finder
//
//  Created by 윤상덕 on 2014. 6. 28..
//  Copyright (c) 2014년 GA. All rights reserved.
//

#ifndef MATOUSEK_H
#define MATOUSEK_H


#include <iostream>
#include <algorithm>
#include <list>
#include <vector>
#include <math.h>
#include "plane.h"

//basis computation에서 사용하는 리턴타입. location 이 upper envelope과 lower envelope의 vertical distance가 최소가 되는점 (x,y,z)이고 distance가 그때의 거리. (점에서 envelope까지의 거리)
class BasisResult
{
public:
    Point location;
    double distance;
    bool all_parallel;

    BasisResult();
    BasisResult(Point a, double d);
};


class Matousek
{
private:
    std::vector<Line> constraints_line; // WARNING : not containing basis
    std::vector<Plane> constraints_plane; // WARNING : not containing basis

public:

    //마투섹 논문의 SUBEX_lp와 동일
    // WARNING : m_plane_list_vector does not contain any basis plane.
    std::vector<Plane> subex_lp(std::vector<Plane> initial_basis, BasisResult &br);

    //Basis computation. basis set에 대해서 upper, lower envelope까지의 거리가 최소인 점과 그때의 거리를 리턴.
    //이를 이용해서 새 베이시스 계산 가능.
    BasisResult basis_computation(std::vector<Plane> &basis, std::vector<Plane> &res_basis);

    //1을 리턴하면 test plane이 basis를 violate, 0이면 반대.
    bool violation_test(BasisResult current_basis, Plane test);


    //위에서 쓴 것들과 완전히 동일한 함수의 2d버전.
    // WARNING : m_plane_list_vector does not contain any basis plane.
    std::vector<Line> subex_lp_2d(std::vector<Line> initial_basis, BasisResult &br);

    BasisResult basis_computation_2d(std::vector<Line> &basis, std::vector<Line> &res_basis);

    bool violation_test_2d(BasisResult current_basisresult, Line test);


    //edge의 양 끝점을 a,b로 함. 둘 중에서 x좌표가 더 작은것을 a로 함!!! 중요함 이거!!!
    std::vector<Line> plane_to_line(std::vector<Plane> &planes, Point a, Point b);

    //constrant를 담고있을 vector
    void insert_constraints_plane(std::vector<Plane> constraints);
    void insert_constraints_line(std::vector<Line> constraints);
};

#endif // MATOUSEK_H
