#include <iostream>
#include <iomanip>
#include <ctime>

#include "common.h"
#include "terrain.h"
#include "translationspacesubdivision.h"


int main(int argc, char **argv)
{
    char ch;

    // ------------- Terrain setting ---------------
    //std::cerr << std::setprecision(std::numeric_limits<double>::digits10);
	std::cerr << std::setprecision(6);
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
    TerrainWithGrids T;
    T.loadData(Tname.c_str());

	//T.sortVerticesAndUpdate();

#ifdef DEBUG
    std::cerr << "Domain information" << std::endl;
	T.print(std::cerr);
#endif


    // ------------- Patch setting ---------------
    std::cerr << "Patch loading..." << std::endl;
		
	TerrainWithGrids P;
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

	//P.sortVerticesAndUpdate();

#ifdef DEBUG
    std::cerr << "Patch information" << std::endl;
	P.print(std::cerr);
#endif


	// ------------- Traverse structure ---------------

	//std::cerr << "Traversal start (s to start) : ";
	//std::cin >> ch;
	//if (ch != 's') throw cpp::Exception("Wrong command.");

	std::cerr << "Traversal start.\n";

    // TSS construction
#ifdef DEBUG
	unsigned int numTotalVertices = 0;
	unsigned int numTotalEdges = 0;
	unsigned int numTotalFaces = 0;
	clock_t clock_TotalInterval = 0;

	clock_t clock_TSSstart = clock();
#endif

	double maxEdgeLength = std::max(T.edgeLength_max, P.edgeLength_max);
	T.makeGrids(maxEdgeLength);
	P.makeGrids(maxEdgeLength);

	TranslationSpaceSubdivision *p_tss;
	try {
		p_tss = new TranslationSpaceSubdivision(&T, &P);
	}
	catch (std::exception& e) {
		std::cerr << "An exception occurred: " << e.what() << std::endl;
		std::cin.sync();
		std::cin.get();
		return 1;
	}
#ifdef DEBUG
	clock_t clock_TSSend = clock();
#endif

	TranslationSpaceSubdivision &tss = *p_tss;

#ifdef DEBUG
	numTotalVertices += tss.arr.number_of_vertices();
	numTotalEdges += tss.arr.number_of_edges();
	numTotalFaces += tss.arr.number_of_faces();
	clock_TotalInterval += clock_TSSend - clock_TSSstart;
#endif

 //   std::cerr << "Translation Space Subdivision information" << std::endl;
 //   std::cerr << "V = " << tss.arr.number_of_vertices() << std::endl;
	//std::cerr << "E = " << tss.arr.number_of_edges() << std::endl;
	//std::cerr << "F = " << tss.arr.number_of_faces() << std::endl;
	//std::cerr << "time elapsed = " << (double)() / (double)CLOCKS_PER_SEC << std::endl;


#ifdef DEBUG
    unsigned int numCellTraversed = 0;
	clock_t clock_TraversalStart = clock();
#endif
    BasisResult global_optimal(Point(0,0,0),std::numeric_limits<double>::infinity()),
                local_optimal(Point(0,0,0),std::numeric_limits<double>::infinity());
    bool first = true;

    //
#ifdef DEBUG
	clock_t clock_IntervalStart = clock();
#endif

	try {
		do {
			std::cerr << "\n<< Next Grid >>\n\n";
			do {
				//std::cin >> ch;
#ifdef DEBUG
				++numCellTraversed;
				if (numCellTraversed % 1 == 0) {
					clock_t clock_IntervalEnd = clock();
					std::cerr << "== Cell " << numCellTraversed << " ==\n";
					std::cerr << "Time elapsed : " << (double)(clock_IntervalEnd - clock_IntervalStart) / (double)CLOCKS_PER_SEC << '\n';
					tss.printCSPlanes();
					std::cerr << "====\n";
					clock_IntervalStart = clock();
					//std::cerr << ".........." << std::endl;
				}
#endif
				bool success = tss.solveLP(local_optimal);

				if (success) {
#ifdef DEBUG
					if (local_optimal.distance < 0) {
						std::cerr << "ERROR: distance < 0" << std::endl;
						std::cerr << "Cell " << numCellTraversed << std::endl;
						std::cerr << "distance=" << local_optimal.distance << std::endl;
						//tss.printCSPlanes();
						char c;
						std::cin >> c;
						continue;
					}
					else if (local_optimal.distance == 0) {
						std::cerr << "WARNING: distance == 0\n";
					}
#endif
					if (first) {
						first = false;
						global_optimal = local_optimal;
					}
					else if (local_optimal.distance < global_optimal.distance)
						global_optimal = local_optimal;
				}
			} while (tss.advanceDFS() == DFS_ADVANCED);
		} while (tss.advanceGridCell() == GridCellSearch_ADVANCED);

	}
	catch (std::exception& e) {
		std::cerr << "An exception occurred: " << e.what() << std::endl;
		std::cin.sync();
		std::cin.get();
		return 1;
	}

	clock_t clock_TraversalEnd = clock();

	// result
#ifdef DEBUG
	std::cerr << "Translation Space Subdivision information" << std::endl;
	std::cerr << "V = " << numTotalVertices << std::endl;
	std::cerr << "E = " << numTotalEdges << std::endl;
	std::cerr << "F = " << numTotalFaces << std::endl;
	std::cerr << "time elapsed = " << (double)(clock_TotalInterval) / (double)CLOCKS_PER_SEC << std::endl;

    std::cerr << "the number of cells traversed = " << numCellTraversed << '\n';
#endif

    std::cerr << "\nGlobal Solution: ("
                << global_optimal.location.x << ","
                << global_optimal.location.y << ","
                << global_optimal.location.z << ") dist="
                << global_optimal.distance << std::endl;
	std::cerr << "Inside number: " << tss.getInsideNum() << '\n';
#ifdef DEBUG
	std::cerr << "Traversal time elapsed : " << (double)(clock_TraversalEnd - clock_TraversalStart) / (double)CLOCKS_PER_SEC << '\n';
#endif

	std::cin.sync();
	std::cin.get();
	return 0;
}
