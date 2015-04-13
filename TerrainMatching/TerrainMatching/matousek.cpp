#include "matousek.h"
#include "numbertype.h"

//#define MATOUSEK_DEBUG

BasisResult::BasisResult()
{
    location=Point();
    distance=0;
    all_parallel=false;
}

BasisResult::BasisResult(Point a, double d)
{
    location = Point(a.x, a.y, a.z);
    distance = d;
    all_parallel = false;
}

void Matousek::insert_constraints_plane(std::vector<Plane> constraints)
{
    constraints_plane = constraints;
//    for (unsigned int i=0; i<constraints.size(); i++)
//    {
//        constraints_plane.push_back(constraints[i]);
//    }
}

void Matousek::insert_constraints_line(std::vector<Line> constraints)
{
    constraints_line = constraints;
//    for (unsigned int i=0; i<constraints.size(); i++)
//    {
//        constraints_line.push_back(constraints[i]);
//    }
}


std::vector<Plane> Matousek::subex_lp(std::vector<Plane> initial_basis, BasisResult &br)
{
#ifdef MATOUSEK_DEBUG
    std::cout << "subex_lp start\n";
    std::cout << "==constraints\n";
    for (unsigned int i=0; i<constraints_plane.size(); i++) {
        constraints_plane[i].print();
    }
    std::cout << "==basis\n";
    for (unsigned int i=0; i<initial_basis.size(); i++) {
        initial_basis[i].print();
    }
#endif

    if (initial_basis.size() < 3) {

#ifdef MATOUSEK_DEBUG
        CGAL_assertion(br.all_parallel);
        std::cout << "subex_lp end. All Parallel !\n";
#endif

        return initial_basis;
    }

    //constraints size가 0이면 그대로 종료
    if (constraints_plane.size() == 0) {

#ifdef MATOUSEK_DEBUG
        std::cout << br.location.x << "," << br.location.y << "," << br.location.z << " " << br.distance << std::endl;
        std::cout << "subex_lp end 0\n";
#endif

        return initial_basis;
    }

#ifdef MATOUSEK_DEBUG
    std::cout<<" # of constraits: "<<constraints_plane.size()<<", # of basis: "<<initial_basis.size()<<std::endl;
#endif

    //constraints에서 random하게 하나를 고름 $$$
    int rand_number = rand()%constraints_plane.size();

    Plane h = constraints_plane[rand_number];

    //골라서 뺌
    constraints_plane.erase(constraints_plane.begin()+rand_number);

    //그걸 뺀 상태로 recursive call.
    // sub_res_basis : res_basis of the sub-problem (with constraints except h).
    std::vector<Plane> sub_res_basis = subex_lp(initial_basis, br);

    if (br.all_parallel) return initial_basis; // all lines are parallel (search boundary lines)

#ifdef MATOUSEK_DEBUG
        std::cout<<"h was : ";
        h.print();
#endif

    //violation test에서 h가 violated
    if (violation_test(br, h))
    {

#ifdef MATOUSEK_DEBUG
        std::cout<<"Violated"<<std::endl;
#endif

        //h를 추가해서...
        sub_res_basis.push_back(h);

        //h를 추가한 상태로 새 basis계산
        std::vector<Plane> new_res_basis;
        br = basis_computation(sub_res_basis, new_res_basis);

#ifdef MATOUSEK_DEBUG
        bool contain = false; // h should be contained in res_basis
        for (std::vector<Plane>::iterator pit = new_res_basis.begin(); pit != new_res_basis.end(); ++pit) {
            pit->print();
            h.print();
            if (*pit == h) {contain = true; break;}
        }
        std::cout<<"Contain = " << contain << std::endl;
        //CGAL_assertion(contain);
        std::cout<<"Size of new result basis: " << new_res_basis.size() << std::endl;
        std::cout<<"Current minimum distance: " << br.distance << std::endl;
#endif

#ifdef MATOUSEK_DEBUG
        std::cout << "subex_lp_2d end 2\n";
#endif

        //정제된 basis로 다시 recursive call
        return subex_lp(new_res_basis, br);
    }
    else
    {
        //violation안된 경우.
        //h 복구
        constraints_plane.push_back(h);

#ifdef MATOUSEK_DEBUG
        std::cout << "subex_lp_2d end 3\n";
#endif

        //리턴
        return sub_res_basis;
    }
}


BasisResult Matousek::basis_computation(std::vector<Plane> &basis, std::vector<Plane> &res_basis)
{
#ifdef MATOUSEK_DEBUG
    std::cout << "\tbc start. basis size=" << basis.size() << std::endl;
    //CGAL_assertion(basis.size() >= 3);

//    for (std::vector<Plane>::iterator pit = basis.begin(); pit != basis.end(); ++pit) {
//        std::cout << pit->a << "," << pit->b;
//    }
//    std::cout << std::endl;
#endif

    std::vector<Line> line_list;
    std::vector<std::pair<unsigned int, unsigned int> > plane_pair_list; // for each point, remember the intersecting plane indexes of basis

    std::vector<Point> point_list;
    std::vector<std::pair<unsigned int, unsigned int> > line_pair_list; // for each point, remember the intersecting line indexes of basis

    const int basis_size = basis.size();

    //주어진 평면들의 교선들의 xy-평면으로의 projection을 구함.
    for (int i=0; i<basis_size-1; i++)
    {
        for (int j=i+1; j<basis_size; j++)
        {
            //두 평면이 평행한 경우는 계산하지 않음
            if (!nearlyEqual(basis[i].a*basis[j].b, basis[i].b*basis[j].a) ||
                !nearlyEqual(basis[i].a*basis[j].c, basis[i].c*basis[j].a) ||
                !nearlyEqual(basis[i].b*basis[j].c, basis[i].c*basis[j].b))
            {
                //computation;
                double a = (basis[i].c*basis[j].a - basis[i].a*basis[j].c)/(basis[i].b*basis[j].c - basis[i].c*basis[j].b);
                double b = (basis[i].d*basis[j].c - basis[i].c*basis[j].d)/(basis[i].b*basis[j].c - basis[i].c*basis[j].b);
                line_list.push_back(Line(a, b));

                std::pair<unsigned int, unsigned int> plane_pair_entry(i, j);
                plane_pair_list.push_back(plane_pair_entry);

#ifdef MATOUSEK_DEBUG
                Line(a, b).print();
#endif
            }
        }
    }

    //구한 교선들에서 교점을 구함.
    const int line_list_size = line_list.size();
    //int count =0;
    for (int i=0; i<line_list_size-1; i++)
    {
        for (int j=i+1; j<line_list_size; j++)
        {
            //std::cout << line_list[i].a << " " << line_list[j].a << std::endl;
            //두 직선이 평행하지 않은 경우면 계산
            if (!nearlyEqual(line_list[i].a, line_list[j].a))
            {
                //계산
                double x = (line_list[i].b - line_list[j].b)/(line_list[j].a - line_list[i].a);
                double y = line_list[i].a*x+line_list[i].b;
                point_list.push_back(Point(x, y, 0.0));

                std::pair<unsigned int, unsigned int> line_pair_entry(i, j);
                line_pair_list.push_back(line_pair_entry);
            }
        }
    }

    Point temp_location;
    double temp_distance=std::numeric_limits<double>::infinity();

    if (point_list.size() == 0) { // the planes in basis are all parallel  or  the lines are all parallel.
        BasisResult b(temp_location, temp_distance);
        b.all_parallel = true;

#ifdef MATOUSEK_DEBUG
    std::cout << "\tbc end parallel\n";
#endif

        return b;
    }

    std::vector<unsigned int> cur_res_basis_idx;

    for (unsigned int i=0; i<point_list.size(); i++)
    {

#ifdef MATOUSEK_DEBUG
        std::cout << "== " <<point_list[i].x << ", " << point_list[i].y << ", " << point_list[i].z << std::endl;
#endif

        std::pair<unsigned int, unsigned int> line_pair_entry(line_pair_list[i]);

        double temp_max=-std::numeric_limits<double>::infinity();
        double temp_min=std::numeric_limits<double>::infinity();

        unsigned int j_max=0; // j that determines temp_max
        unsigned int j_min=0; // j that determines temp_min

        for (unsigned int j=0; j<basis.size(); j++)
        {
            //각 평면의 주어진 점에서의 z좌표를 구함. 제일 높은걸 temp_max에, 제일 낮은걸 temp_min에 저장. 나중에 다 계산 뒤에 temp_max-temp_min으로 주어진 점에서의 upper envolope과 lower envelope의 vertical distance를 구함.
            double plane_z = (basis[j].d - basis[j].a*point_list[i].x - basis[j].b*point_list[i].y)/basis[j].c;
            if (temp_max < plane_z) {
                temp_max = plane_z;
                j_max = j;

#ifdef MATOUSEK_DEBUG
//                std::cout << "max update : temp_max=" << temp_max << ", j_max=" << j << std::endl;
#endif
            }
            if (temp_min > plane_z) {
                temp_min = plane_z;
                j_min = j;

#ifdef MATOUSEK_DEBUG
//                std::cout << "min update : temp_min=" << temp_min << ", j_min=" << j << std::endl;
#endif
            }
        }

        double distance_i = (temp_max - temp_min)/2.0;

#ifdef MATOUSEK_DEBUG
        CGAL_assertion(temp_max >= temp_min);
        std::cout << "temp_distance=" << temp_distance << ", distance_i=" << distance_i << std::endl;
#endif

        if (temp_distance > distance_i)
        {
            std::pair<unsigned int, unsigned int> line1_plane_pair_entry = plane_pair_list[line_pair_entry.first];
            std::pair<unsigned int, unsigned int> line2_plane_pair_entry = plane_pair_list[line_pair_entry.second];

            // collect related planes
            cur_res_basis_idx.clear();
            cur_res_basis_idx.push_back(line1_plane_pair_entry.first);
            cur_res_basis_idx.push_back(line1_plane_pair_entry.second);
            cur_res_basis_idx.push_back(line2_plane_pair_entry.first);
            cur_res_basis_idx.push_back(line2_plane_pair_entry.second);
            cur_res_basis_idx.push_back(j_max);
            cur_res_basis_idx.push_back(j_min);

#ifdef MATOUSEK_DEBUG
            for (std::vector<unsigned int>::iterator piit = cur_res_basis_idx.begin(); piit != cur_res_basis_idx.end(); ++piit)
                std::cout << *piit << " ";
            std::cout << distance_i << std::endl;
            line_list[line_pair_entry.first].print();
            line_list[line_pair_entry.second].print();
#endif

            // remove duplicates
            std::sort(cur_res_basis_idx.begin(), cur_res_basis_idx.end());
            std::vector<unsigned int>::iterator last = std::unique(cur_res_basis_idx.begin(), cur_res_basis_idx.end());
            cur_res_basis_idx.erase(last, cur_res_basis_idx.end());

            temp_distance = distance_i;
            temp_location = Point(point_list[i].x, point_list[i].y, (temp_max + temp_min)/2.0);
        }
    }

    //CGAL_assertion(cur_res_basis_idx.size() >= 4);

    // plane index to plane
    for (std::vector<unsigned int>::iterator piit = cur_res_basis_idx.begin(); piit != cur_res_basis_idx.end(); ++piit) {
        res_basis.push_back(basis[*piit]);

#ifdef MATOUSEK_DEBUG
        basis[*piit].print();
#endif
    }

    // result basis info
    BasisResult b(temp_location, temp_distance);
    b.all_parallel = false;

#ifdef MATOUSEK_DEBUG
    std::cout << "\tbc end. dist=" << temp_distance << ", loc=" << temp_location.x << "," << temp_location.y << "," << temp_location.z << "\n";
#endif

    return b;
}

bool Matousek::violation_test(BasisResult current_basisresult, Plane test)
{
    double test_z_in_res = (test.d - current_basisresult.location.x * test.a - current_basisresult.location.y * test.b)/test.c;
    double vertical_dist = test_z_in_res - current_basisresult.location.z;

#ifdef MATOUSEK_DEBUG
    std::cout << "-- test_z_in_res=" << test_z_in_res << ", vectical_dist=" << vertical_dist << std::endl;
#endif

    return (vertical_dist > current_basisresult.distance ||
            -vertical_dist > current_basisresult.distance);
}


std::vector<Line> Matousek::subex_lp_2d(std::vector<Line> initial_basis, BasisResult &br)
{
#ifdef MATOUSEK_DEBUG
    std::cout << "subex_lp_2d start\n";
    std::cout << "==constraints\n";
    for (unsigned int i=0; i<constraints_line.size(); i++) {
        constraints_line[i].print();
    }
    std::cout << "==basis\n";
    for (unsigned int i=0; i<initial_basis.size(); i++) {
        initial_basis[i].print();
    }
#endif

    if (initial_basis.size() < 2) {

#ifdef MATOUSEK_DEBUG
        CGAL_assertion(br.all_parallel);
        std::cout << "subex_lp end. All Parallel !\n";
#endif

        return initial_basis;
    }

    //constraints가 없으면 그대로 종료 initial basis로 basis계산하고 종료.
    if (constraints_line.size() == 0) {

#ifdef MATOUSEK_DEBUG
        std::cout << "subex_lp_2d end 0\n";
#endif

        return initial_basis;
    }

#ifdef MATOUSEK_DEBUG
    std::cout<<" # of constraits: "<<constraints_line.size()<<", # of basis: "<<initial_basis.size()<<std::endl;
#endif

    //constraints에서 random하게 하나를 고름
    int rand_number = rand()%constraints_line.size();

    Line h = constraints_line[rand_number];

    //골라서 뺌
    constraints_line.erase(constraints_line.begin()+rand_number);

    //그걸 뺀 상태로 recursive call.
    std::vector<Line> sub_res_basis = subex_lp_2d(initial_basis, br);

    if (br.all_parallel) return initial_basis; // all lines are parallel (search boundary pts)

#ifdef MATOUSEK_DEBUG
    if (initial_basis.size() < 2) {
        char c;
        std::cin >> c;
    }
#endif

    //violation test에서 h가 violated
    if (violation_test_2d(br, h))
    {

#ifdef MATOUSEK_DEBUG
        std::cout<<"Violated"<<std::endl;
#endif

        //h를 추가해서...
        sub_res_basis.push_back(h);

        //cout<<h.a<<" "<<h.b<<endl;

        //h를 추가한 상태로 새 basis계산
        std::vector<Line> new_res_basis;
        br = basis_computation_2d(sub_res_basis, new_res_basis);

#ifdef MATOUSEK_DEBUG
        bool contain = false; // h should be contained in res_basis
        for (std::vector<Line>::iterator pit = new_res_basis.begin(); pit != new_res_basis.end(); ++pit) {
            if (*pit == h) {contain = true; break;}
        }
        //CGAL_assertion(contain);
        std::cout<<"Contain = " << contain << std::endl;
        std::cout<<"Size of new result basis: " << new_res_basis.size() <<std::endl;
        std::cout<<"Current minimum distance: " << br.distance << std::endl;
#endif

#ifdef MATOUSEK_DEBUG
        std::cout << "subex_lp_2d end 2\n";
#endif

        //h를 복구하고 정제된 basis로 다시 recursive call
        return subex_lp_2d(new_res_basis, br);
    }
    else
    {
        //violation안된 경우. h를 복구하고 리턴
        constraints_line.push_back(h);

#ifdef MATOUSEK_DEBUG
        std::cout << "subex_lp_2d end 3\n";
#endif

        return sub_res_basis;
    }
}


BasisResult Matousek::basis_computation_2d(std::vector<Line> &basis, std::vector<Line> &res_basis)
{
#ifdef MATOUSEK_DEBUG
    std::cout << "\tbc2d start\n";
    //CGAL_assertion(basis.size() >= 3);

    for (std::vector<Line>::iterator lit = basis.begin(); lit != basis.end(); ++lit) {
        std::cout << lit->a << " & ";
    }
    std::cout << std::endl;
#endif

    std::vector<Point> point_list;
    std::vector<std::pair<unsigned int, unsigned int> > basis_pair_list; // for each point, remember the intersecting line indexes of basis

    //주어진 basis의 size
    const int basis_size = basis.size();

    for (int i=0; i<basis_size-1; i++)
    {
        for (int j=i+1; j<basis_size; j++)
        {
            //두 직선이 평행하지 않은 경우면 계산
            if (!nearlyEqual(basis[i].a, basis[j].a))
            {
                //계산
                double x = (basis[i].b - basis[j].b)/(basis[j].a - basis[i].a);
                point_list.push_back(Point(x, 0.0, 0.0));

                std::pair<unsigned int, unsigned int> basis_pair_entry(i, j);
                basis_pair_list.push_back(basis_pair_entry);
                //std::cout<<"point"<<" x="<<point_list.back().x<<std::endl;
            }
            else{
                //std::cout<<"parallel"<<std::endl;
            }

        }
    }

    Point temp_location;
    double temp_distance=std::numeric_limits<double>::infinity();

    if (point_list.size() == 0) { // the lines in basis are all parallel
        BasisResult b(temp_location, temp_distance);
        b.all_parallel = true;
        return b;
    }

    std::vector<unsigned int> cur_res_basis_idx;

    for (unsigned int i=0; i<point_list.size(); i++)
    {
        std::pair<unsigned int, unsigned int> basis_pair_entry(basis_pair_list[i]); // two lines corresponding to the point

        double temp_max=-std::numeric_limits<double>::infinity();
        double temp_min=std::numeric_limits<double>::infinity();

        unsigned int j_max=0; // j that determines temp_max
        unsigned int j_min=0; // j that determines temp_min

        for (unsigned int j=0; j<basis.size(); j++)
        {
            //각 직선의 주어진 점에서의 y좌표를 구함. 제일 높은걸 temp_max에, 제일 낮은걸 temp_min에 저장. 나중에 다 계산 뒤에 temp_max-temp_min으로 주어진 점에서의 upper envolope과 lower envelope의 vertical distance를 구함.
            double line_z = basis[j].a * point_list[i].x + basis[j].b;
            if (temp_max < line_z) {
                temp_max = line_z;
                j_max = j;
            }
            if (temp_min > line_z) {
                temp_min = line_z;
                j_min = j;
            }
        }

        double distance_i = (temp_max - temp_min)/2.0;

        if (temp_distance > distance_i)
        {
            // collect related planes
            cur_res_basis_idx.clear();
            cur_res_basis_idx.push_back(basis_pair_entry.first);
            cur_res_basis_idx.push_back(basis_pair_entry.second);
            cur_res_basis_idx.push_back(j_max);
            cur_res_basis_idx.push_back(j_min);

            // remove duplicates
            std::sort(cur_res_basis_idx.begin(), cur_res_basis_idx.end());
            std::vector<unsigned int>::iterator last = std::unique(cur_res_basis_idx.begin(), cur_res_basis_idx.end());
            cur_res_basis_idx.erase(last, cur_res_basis_idx.end());

            temp_distance = distance_i;
            temp_location = Point(point_list[i].x, 0.0, (temp_max + temp_min)/2.0);
        }
    }

    //CGAL_assertion(cur_res_basis_idx.size() >= 3);

    // line index to line
    for (std::vector<unsigned int>::iterator liit = cur_res_basis_idx.begin(); liit != cur_res_basis_idx.end(); ++liit) {
        res_basis.push_back(basis[*liit]);

#ifdef MATOUSEK_DEBUG
        basis[*liit].print();
#endif
    };

    // result basis info
    BasisResult b(temp_location, temp_distance);
    b.all_parallel = false;

#ifdef MATOUSEK_DEBUG
    std::cout << "\tbc2d end " << b.all_parallel << " dist=" << temp_distance << "\n";
//    char c;
//    if (cur_res_basis_idx.size() < 3) std::cin >> c;
#endif

    return b;
}

bool Matousek::violation_test_2d(BasisResult current_basisresult, Line test)
{
    double test_z_in_res = test.a*current_basisresult.location.x + test.b;
    double vertical_dist = test_z_in_res - current_basisresult.location.z;

#ifdef MATOUSEK_DEBUG
    std::cout << "-- test_z_in_res=" << test_z_in_res << ", vectical_dist=" << vertical_dist << std::endl;
#endif

    return (vertical_dist > current_basisresult.distance ||
            -vertical_dist > current_basisresult.distance);
}

std::vector<Line> Matousek::plane_to_line(std::vector<Plane> &planes, Point a, Point b)
{
    //두 점 a,b를 지나나고 z축에 평행한 평면h와 planes들의 교선들을 구한다.
    //구한 교선은 h상에서 정의되어있다고 가정하고, a를 h평면상의 원점으로 가정함.

    //두 점 a, b를 지나고 z축에 평행한 평면을 만든다.
    Plane h = Plane(a.y-b.y, b.x-a.x, 0.0, a.x*(a.y-b.y)+a.y*(b.x-a.x));

    std::vector<Line> intersecting_lines;

    //planes에 들어있는 각각의 평면과 h와의 교선을 구함.
    for(unsigned int i=0;i<planes.size();i++)
    {
        //h가 z축과 평행한 녀석이고 평행한 건 들어오지 않는다고 가정할 수 있으므로 체크는 안하는것으로..
        /*
        if (planes[i].a/h.a != planes[i].b/h.b ||
            planes[i].a/h.a != planes[i].c/h.c ||
            planes[i].b/h.b != planes[i].c/h.c)
        {
         */
            //computation

            //기울기
            //3차원 공간상에서 직선의 기울기는
            //x : y : z = planes[i].c * h.b : - planes[i].c*h.a :h.a*planes[i].b - planes[i].a*h.b
            //따라서 h상에서 직선의 기울기는 z / sqrt(x^2 + y^2)가 됨.

            //sqrt를 사용하기 위해서 잠시 mpf_t로 만들었다가 다시 변환.
        double squared = planes[i].c*planes[i].c*(h.a*h.a + h.b*h.b);

//        mpf_t temp;
//        mpf_init2(temp, 10000);
//        mpf_set_q(temp, squared.mpq());

//        mpf_t sqrt;
//        mpf_init2(sqrt, 10000);
//        mpf_sqrt(sqrt, temp);

//        mpq_t temp2;
//        mpq_init(temp2);
//        mpq_set_f(temp2, sqrt);

//        double denom(temp2);

        double denom = sqrt(squared);

        //test
        //temp2.print();
        //cout<<endl;


        //mpf_clears(temp,sqrt);

            //double ratio = (h.a*planes[i].b - planes[i].a*h.b)/((planes[i].c*planes[i].c*(h.a*h.a + h.b*h.b)).sqrt());
            double ratio = (h.a*planes[i].b - planes[i].a*h.b)/denom;

        //우리가 구하는 평면은 x가 증가하는 방향, 혹은 x가 같은 경우 y가 증가하는 방향을 양의 방향이라고 정의한다.

        if(planes[i].c*h.b != 0.0)
        {
            if (planes[i].c*h.b < 0.0)// && h.a*planes[i].b - planes[i].a*h.b > 0.0)
            //if (planes[i].c*h.b * h.a< 0.0)
            {
                ratio = - ratio;
            }

        //}
        }
        else
        {
            if(- planes[i].c*h.a < 0.0)// && h.a*planes[i].b - planes[i].a*h.b > 0.0)
            {
                ratio = -ratio;
            }
        }

        intersecting_lines.push_back(Line(ratio, (planes[i].d - planes[i].a * a.x - planes[i].b * a.y)/planes[i].c));
    }


    return intersecting_lines;
}
