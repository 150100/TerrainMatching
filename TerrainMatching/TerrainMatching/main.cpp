#include <iostream>
#include <iomanip>
#include <ctime>

#include "terrain.h"
#include "translationspacesubdivision.h"


int main(int argc, char **argv)
{
    char ch;

    // ------------- Terrain setting ---------------
    std::cerr << std::setprecision(std::numeric_limits<double>::digits10);
    std::cerr << "Domain loading..." << std::endl;

    // two vectors to hold point coordinates and triangle vertex indices
	std::string Tname, Pname;
    bool part_selection = false;

    // test : arguments
    if (argc > 3) throw cpp::Exception("Invalid arguments.");

	// file name
	if (argc == 1) {
		std::cerr << "Enter filename of T : ";
		std::cin >> Tname;
		std::cerr << "Enter filename of P : ";
		std::cin >> Pname;
	}
	else if (argc == 2) {
		part_selection = true;
		Tname = argv[1];
	}
	else {
		Tname = argv[1];
		Pname = argv[2];
	}

    // build a terrain T
    Terrain T;
    T.loadData(Tname.c_str());

    std::cerr << "Domain information" << std::endl;
    std::cerr << "V = " << T.number_of_vertices() << std::endl;
    std::cerr << "E = " << T.number_of_edges() << std::endl;
    std::cerr << "F = " << T.number_of_faces() << std::endl;


    // ------------- Patch setting ---------------
    std::cerr << "Patch loading..." << std::endl;
		
    Terrain P;
	if (part_selection) {
		double part_x_min, part_x_max, part_y_min, part_y_max;
		std::cerr << "range of x (0~1) (ex- 0.3 0.4) :";
		std::cin >> part_x_min >> part_x_max;
		std::cerr << "range of y (0~1) (ex- 0.3 0.4) :";
		std::cin >> part_y_min >> part_y_max;
		T.createGetTrimTerrain(part_x_min, part_x_max, part_y_min, part_y_max, &P);
	}
    else
        P.loadData(Pname.c_str());

    std::cerr << "Patch information" << std::endl;
    std::cerr << "V = " << P.number_of_vertices() << std::endl;
    std::cerr << "E = " << P.number_of_edges() << std::endl;
    std::cerr << "F = " << P.number_of_faces() << std::endl;



    // ------------- TSS setting ---------------

    std::cerr << "Subdividing Translation Space..." << std::endl;

    // TSS construction
	clock_t clock_TSSstart = clock();

	TranslationSpaceSubdivision *p_tss;
	try {
		p_tss = new TranslationSpaceSubdivision(&T, &P); // with construction
	}
	catch (std::exception& e) {
		std::cerr << "An exception occurred: " << e.what() << std::endl;
		std::cin.sync();
		std::cin.get();
		return 1;
	}
	clock_t clock_TSSend = clock();

	TranslationSpaceSubdivision &tss = *p_tss;

    std::cerr << "Translation Space Subdivision information" << std::endl;
    std::cerr << "V = " << tss.arr.number_of_vertices() << std::endl;
	std::cerr << "E = " << tss.arr.number_of_edges() << std::endl;
	std::cerr << "F = " << tss.arr.number_of_faces() << std::endl;
	std::cerr << "time elapsed = " << (double)(clock_TSSend - clock_TSSstart) / (double)CLOCKS_PER_SEC << std::endl;


    // ------------- Traverse structure ---------------

    std::cerr << "Traversal start (s to start) : ";
    std::cin >> ch;
	if (ch != 's') throw cpp::Exception("Wrong command.");

    unsigned num = 0;
    BasisResult global_optimal(Point(0,0,0),std::numeric_limits<double>::infinity()),
                local_optimal(Point(0,0,0),std::numeric_limits<double>::infinity());
    bool first = true;

    //
	clock_t clock_IntervalStart = clock();
	do {
		//std::cin >> ch;
		++num;
		if (num % 1 == 0) {
			clock_t clock_IntervalEnd = clock();
			std::cerr << "\n== Cell " << num << " ==\n";
			std::cerr << "Time elapsed : " << (double)(clock_IntervalEnd - clock_IntervalStart) / (double)CLOCKS_PER_SEC << '\n';
			clock_IntervalStart = clock();
			tss.printCSPlanes();
			//std::cerr << ".........." << std::endl;
		}

        bool success = tss.solveLP(local_optimal);

        if (success) {
            if (local_optimal.distance < 0) {
                std::cerr << "ERROR: distance < 0" << std::endl;
                std::cerr << "Cell " << num << std::endl;
                std::cerr << "distance=" << local_optimal.distance << std::endl;
                //tss.printCSPlanes();
                char c;
                std::cin >> c;
                continue;
            }
            else if (local_optimal.distance == 0) {
                std::cerr << "WARNING: distance == 0\n";
            }

            if (first) {
                first = false;
                global_optimal = local_optimal;
            }
            else if (local_optimal.distance < global_optimal.distance)
                global_optimal = local_optimal;
        }
    }
    while (tss.advance() == tss.DFS_ADVANCED);

	// result
    std::cerr << "the number of cells traversed = " << num << std::endl;

    std::cerr << "\nGlobal Solution: ("
                << global_optimal.location.x << ","
                << global_optimal.location.y << ","
                << global_optimal.location.z << ") dist="
                << global_optimal.distance << std::endl;
    std::cerr << "Inside number: " << tss.getInsideNum() << std::endl;

	std::cin.sync();
	std::cin.get();
	return 0;
}
