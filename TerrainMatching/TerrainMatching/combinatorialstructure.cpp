#include "combinatorialstructure.h"

BasisResult CombinatorialStructure::solveLP(std::vector<Point> boundary_pts)
{
    BasisResult br(Point(), std::numeric_limits<double>::infinity());
    for (std::vector<Point>::iterator bit = boundary_pts.begin(); bit != boundary_pts.end(); ++bit) {
        double z_min(std::numeric_limits<double>::infinity()), z_max(-std::numeric_limits<double>::infinity());
        for (std::list<Plane>::iterator pit = m_plane_list.begin(); pit != m_plane_list.end(); ++pit) {
            double z_p = (pit->d - (pit->a * bit->x) - (pit->b * bit->y)) / pit->c;
            if (z_p < z_min) z_min = z_p;
            if (z_p > z_max) z_max = z_p;
        }

        if (z_min == std::numeric_limits<double>::infinity())
            throw cpp::Exception("z_min not updated. (CombinatorialStructure::solveLP)");
        if (z_max == -std::numeric_limits<double>::infinity())
            throw cpp::Exception("z_max not updated. (CombinatorialStructure::solveLP)");

        double z_diff = (z_max - z_min) / 2.0;
        if (z_diff <= 0)
            throw cpp::Exception("z_diff should be positive. (CombinatorialStructure::solveLP)");

        double z_mid = z_min + z_diff;
        if (z_diff < br.distance) {
            br.location = Point(bit->x, bit->y, z_mid);
            br.distance = z_diff;
        }
    }
    return br;


//    std::vector<Plane> m_plane_list_vector(m_plane_list.begin(), m_plane_list.end());
//    std::random_shuffle(m_plane_list_vector.begin(), m_plane_list_vector.end());

//    std::vector<Plane> m_all_plane_list_vector = m_plane_list_vector; // maintaining all planes for boundary condition

//    std::vector<Plane> m_basis_list_vector;
//    m_basis_list_vector.push_back(m_plane_list_vector[0]);
//    m_plane_list_vector.erase(m_plane_list_vector.begin());

//    unsigned int cur_number = 1;
//    unsigned int pl_index = 1;
//    while (cur_number < 4 && pl_index < m_plane_list_vector.size()) {
//        Plane *p_p1, *p_p2;
//        p_p2 = &m_plane_list_vector[pl_index];
//        bool has_parallel = false;
//        for (unsigned int bl_index = 0; bl_index < m_basis_list_vector.size(); ++bl_index)
//        {
//            p_p1 = &m_basis_list_vector[bl_index];
//            if ( equal(p_p1->a * p_p2->b, p_p1->b * p_p2->a) &&
//                 equal(p_p1->b * p_p2->c, p_p1->c * p_p2->b) &&
//                 equal(p_p1->a * p_p2->c, p_p1->c * p_p2->a))
//            {
//                has_parallel = true;
//                break;
//            }
//        }

//        if (!has_parallel) { // m_plane_list_vector does not contain any basis plane.
//            m_basis_list_vector.push_back(*p_p2);
//            m_plane_list_vector.erase(m_plane_list_vector.begin() + pl_index);
//            cur_number++;
//        }
//        else
//            pl_index++;
//    }

////    std::cout << "Solve start: number=" << m_plane_list_vector.size() << " " << m_basis_list_vector.size() << std::endl;

////    std::cout << std::endl << "[plane]\n";
////    for (std::vector<Plane>::iterator pit = m_plane_list_vector.begin(); pit != m_plane_list_vector.end(); ++pit) {
////        pit->print();
////    }
////    std::cout << std::endl;

////    std::cout << std::endl << "[basis]\n";
////    for (std::vector<Plane>::iterator pit = m_basis_list_vector.begin(); pit != m_basis_list_vector.end(); ++pit) {
////        pit->print();
////    }
////    std::cout << std::endl;

//    Matousek m;
//    m.insert_constraints_plane(m_plane_list_vector);
//    std::vector<Plane> temp_basis;
//    BasisResult br = m.basis_computation(m_basis_list_vector, temp_basis);

//    bool outside = false; // when the solution is outside of the cell, true.

//    if (!br.all_parallel) {

////        CGAL_assertion(temp_basis.size() == m_basis_list_vector.size());
////        std::vector<Plane>::iterator pit1 = temp_basis.begin();
////        std::vector<Plane>::iterator pit2 = m_basis_list_vector.begin();
////        while (pit1 != temp_basis.end()) {
////            CGAL_assertion(*pit1 == *pit2);
////            ++pit1;
////            ++pit2;
////        }

//        std::vector<Plane> res_basis = m.subex_lp(m_basis_list_vector, br);
////        for (std::vector<Plane>::iterator bit = res_basis.begin(); bit != res_basis.end(); ++bit)
////            bit->print();

//        // check the solution is outside the cell
//        for (std::vector<Point>::iterator bit = boundary_pts.begin(); bit != boundary_pts.end(); ++bit) {
//            std::vector<Point>::iterator prev_bit;
//            if (bit == boundary_pts.begin()) {prev_bit = boundary_pts.end(); --prev_bit;}
//            else {prev_bit = bit; --prev_bit;}

//            Point b1(*prev_bit), b2(*bit);
//            if ((b2.x-b1.x)*(br.location.y-b1.y)-(b2.y-b1.y)*(br.location.x-b1.x) < 0) {
//                outside = true; // right side of b1b2 is outside
//                break;
//            }
//        }
//    }

//    // if the solution is outside the boundary
//    if (outside || br.all_parallel) {
//        //std::cout << "Solution is outside the cell." << std::endl;
//        bool first = true;
//        for (std::vector<Point>::iterator bit = boundary_pts.begin(); bit != boundary_pts.end(); ++bit) {
//            //std::cout << "(" << bit->x << "," << bit->y << "," << bit->z << ") \n";
//            std::vector<Point>::iterator prev_bit;
//            if (bit == boundary_pts.begin()) {prev_bit = boundary_pts.end(); --prev_bit;}
//            else {prev_bit = bit; --prev_bit;}

//            Matousek m_l;

//            // two boundary points
//            Point b1(*prev_bit), b2(*bit);
//            if (b1.x > b2.x) {Point temp = b1; b1 = b2; b2 = temp;}
//            //else if (b1.x == b2.x) std::cout << "Error: Same x-coordinate. " << b1.x.to_double() << std::endl;
//            else if (b1.x == b2.x) std::cout << "Error: Same x-coordinate. " << b1.x << std::endl;

//            // basis and constraints
//            std::vector<Line> line_vector, basis_vector;

//            // all planes
//            line_vector = m_l.plane_to_line(m_all_plane_list_vector, b1, b2);

////            for (std::vector<Line>::iterator lit = line_vector.begin(); lit != line_vector.end(); ++lit) {
////                std::cout << lit->a << "x + " << lit->b << " = 0" << std::endl;
////            }
////            std::cout << std::endl;


////            std::cout << std::endl << "[all line]\n";
////            for (std::vector<Line>::iterator pit = line_vector.begin(); pit != line_vector.end(); ++pit) {
////                pit->print();
////            }
////            std::cout << std::endl;

//            basis_vector.push_back(line_vector[0]);
//            line_vector.erase(line_vector.begin());
//            unsigned int cur_number = 1;
//            unsigned int pl_index = 1;
//            while (cur_number < 3 && pl_index < line_vector.size()) {
//                Line *p_l1, *p_l2;
//                p_l2 = &line_vector[pl_index];
//                bool has_parallel = false;
//                for (unsigned int bl_index = 0; bl_index < basis_vector.size(); ++bl_index)
//                {
//                    p_l1 = &basis_vector[bl_index];
//                    if (equal(p_l1->a, p_l2->a)) {
//                        has_parallel = true;
//                        break;
//                    }
//                }

//                if (!has_parallel) { // line_vector does not contain any basis line.
//                    basis_vector.push_back(*p_l2);
//                    line_vector.erase(line_vector.begin() + pl_index);
//                    ++cur_number;
//                }
//                else
//                    ++pl_index;
//            }


////            std::cout << std::endl << "[line]\n";
////            for (std::vector<Line>::iterator pit = line_vector.begin(); pit != line_vector.end(); ++pit) {
////                pit->print();
////            }
////            std::cout << std::endl;

////            std::cout << std::endl << "[basis]\n";
////            for (std::vector<Line>::iterator pit = basis_vector.begin(); pit != basis_vector.end(); ++pit) {
////                pit->print();
////            }
////            std::cout << std::endl;

//            // compute local solution
//            m_l.insert_constraints_line(line_vector);
//            std::vector<Line> temp_basis_line;
//            BasisResult br_h = m_l.basis_computation_2d(basis_vector, temp_basis_line);

////            CGAL_assertion(temp_basis_line.size() == basis_vector.size());
////            std::vector<Line>::iterator pit1 = temp_basis_line.begin();
////            std::vector<Line>::iterator pit2 = basis_vector.begin();
////            while (pit1 != temp_basis_line.end()) {
////                CGAL_assertion(*pit1 == *pit2);
////                ++pit1;
////                ++pit2;
////            }

//            if (!br.all_parallel)
//                std::vector<Line> res_basis_line = m_l.subex_lp_2d(basis_vector, br_h); // b_cur is in coordinate of plane "h" in basis_computation_2d


//            // if the solution is outside the boundary, pick solution as a point
//            double b1b2_abs = sqrt((b1.x-b2.x)*(b1.x-b2.x) + (b1.y-b2.y)*(b1.y-b2.y)); // | (segment) b1 b2 |
//            BasisResult br_cur;

//            if (br_h.all_parallel || br_h.location.x < 0 || br_h.location.x > b1b2_abs) {
//                //std::cout << " Solution is outside the line boundary." << std::endl;
//                double min_z_b1(std::numeric_limits<double>::infinity()), max_z_b1(-std::numeric_limits<double>::infinity());
//                double min_z_b2(std::numeric_limits<double>::infinity()), max_z_b2(-std::numeric_limits<double>::infinity());

//                for (std::vector<Line>::iterator lit = line_vector.begin(); lit != line_vector.end(); ++lit) {
//                    double z_b1 = lit->b;
//                    if (z_b1 < min_z_b1) min_z_b1 = z_b1;
//                    if (z_b1 > max_z_b1) max_z_b1 = z_b1;

//                    double z_b2 = lit->a * b1b2_abs + lit->b;
//                    if (z_b2 < min_z_b2) min_z_b2 = z_b2;
//                    if (z_b2 > max_z_b2) max_z_b2 = z_b2;
//                }

//                double b1_dist = (max_z_b1-min_z_b1)/2.0;
//                double b2_dist = (max_z_b2-min_z_b2)/2.0;
//                if (b1_dist < b2_dist)
//                    br_cur = BasisResult(Point(b1.x, b1.y, (max_z_b1-min_z_b1)/2.0), b1_dist);
//                else
//                    br_cur = BasisResult(Point(b1.x, b1.y, (max_z_b2-min_z_b2)/2.0), b2_dist);

//                //std::cout << " (" << br_cur.location.x << "," << br_cur.location.y << "," << br_cur.location.z << ") dist=" << br_cur.distance << std::endl;
//            }
//            else {
//                //std::cout << " Solution is inside the line boundary." << std::endl;
//                CGAL_assertion(br_h.location.y == 0);
//                double ratio = br_h.location.x / b1b2_abs;
//                CGAL_assertion(ratio >= 0 && ratio <= 1);
//                br_cur = BasisResult( Point((b2.x-b1.x)*ratio+b1.x, (b2.y-b1.y)*ratio+b1.y, br_h.location.z), br_h.distance );
//            }

//            if (first) {first = false; br = br_cur;}
//            else if (br_cur.distance < br.distance) br = br_cur;
//        }
//    }
//    else {
//        //std::cout << "Solution is inside the cell." << std::endl;
//        inside_num++;
//    }

////    std::cout << "Solve end: ("
////              << br.location.x << ","
////              << br.location.y << ","
////              << br.location.z << ") dist="
////              << br.distance << std::endl;

////    char c;
////    std::cin >> c;

    //return br;
}

void CombinatorialStructure::printPlanes()
{
    std::list<Plane>::iterator it;
    for (it = m_plane_list.begin(); it != m_plane_list.end(); ++it) {
        it->print();
    }
}


Plane EdgeEdgePair::makePlane()
{
    Point &a_1 = edge1.he->getOrigin()->getData().p;
    Point &a_2 = edge1.he->getTwin()->getOrigin()->getData().p;
    Point &b_1 = edge2.he->getOrigin()->getData().p;
    Point &b_2 = edge2.he->getTwin()->getOrigin()->getData().p;

    Plane p;
    if (edge1_is_from_patch)
       p.setPlane_EEpair(a_1, a_2, b_1, b_2, 1);
    else
       p.setPlane_EEpair(a_1, a_2, b_1, b_2, -1);

    return p;
}


Plane VertexTrianglePair::makePlane()
{
    HalfEdge *boundary_he = pTriangle->getBoundary();
    HalfEdge *first_boundary_he = boundary_he;
    Point &a_1 = boundary_he->getOrigin()->getData().p;
    boundary_he = boundary_he->getNext();
    Point &a_2 = boundary_he->getOrigin()->getData().p;
    boundary_he = boundary_he->getNext();
    Point &a_3 = boundary_he->getOrigin()->getData().p;
    if (first_boundary_he != boundary_he->getNext()) throw cpp::Exception("Face is not a triangle.");
    Point &b_1 = pVertex->getData().p;

    Plane p;
    if (tri_is_from_patch)
       p.setPlane_VTpair(a_1, a_2, a_3, b_1, 1);
    else
       p.setPlane_VTpair(a_1, a_2, a_3, b_1, -1);

    return p;
}

