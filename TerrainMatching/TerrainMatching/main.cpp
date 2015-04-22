#include <iostream>
#include <iomanip>

#include "terrain.h"
#include "translationspacesubdivision.h"


int main(int argc, char **argv)
{
    try {

    //    QGraphicsScene scene;
    //    scene.setSceneRect(0,0,1024,800);
    ////    scene.addRect(QRectF(0,0,1024,800));
    ////    scene.addLine(QLineF(0,0,1024,800));
    ////    scene.addLine(QLineF(0,800,1024,0));

    //    // setting Qt view
    //    QGraphicsView* view = new QGraphicsView(&scene);
    //    CGAL::Qt::GraphicsViewNavigation navigation;
    //    view->installEventFilter(&navigation);
    //    view->viewport()->installEventFilter(&navigation);
    //    view->setRenderHint(QPainter::Antialiasing);

        // setting CGAL timer
        char ch;

        // ------------- Terrain setting ---------------
        std::cerr << std::setprecision(std::numeric_limits<double>::digits10);
        std::cerr << "Domain loading..." << std::endl;

        // two vectors to hold point coordinates and triangle vertex indices
        std::vector<double> T_coords;
        std::vector<unsigned int> T_tris;
        std::vector<double> P_coords;
        std::vector<unsigned int> P_tris;
        bool part_selection = false;

        // test : arguments
        if (argc > 3) throw cpp::Exception("Invalid arguments.");

		// file name
		std::string Tname, Pname;
		if (argc == 1) {
			std::cerr << "Enter filename of T : ";
			std::cin >> Tname;
			std::cerr << "Enter filename of P : ";
			std::cin >> Pname;
		}
		else {
			Tname = argv[1];
			Pname = argv[2];
		}

        // build a terrain T
        Terrain T;
        if (argc == 1)
            T.loadData(T_coords, T_tris);
        else
            T.loadData(argv[1]);

        std::cerr << "Domain information" << std::endl;
        std::cerr << "V = " << T.number_of_vertices() << std::endl;
        std::cerr << "E = " << T.number_of_edges() << std::endl;
        std::cerr << "F = " << T.number_of_faces() << std::endl;

  //      std::cerr << "c to continue: ";
  //      std::cin >> ch;
		//if (ch != 'c') throw cpp::Exception("Wrong command.");



        // ------------- Patch setting ---------------
        std::cerr << "Patch loading..." << std::endl;

        double part_x_min, part_x_max, part_y_min, part_y_max;

        if (part_selection) {
            P_coords = T_coords;
            P_tris = T_tris;

            std::cerr << "range of x (0~1) (ex- 0.3 0.4) :";
            std::cin >> part_x_min >> part_x_max;
            std::cerr << "range of y (0~1) (ex- 0.3 0.4) :";
            std::cin >> part_y_min >> part_y_max;
        }

        Terrain P;
        if (argc == 1)
            P.loadData(P_coords, P_tris);
//        else if (part_selection)
//            P = T.createGetTrimTerrain(part_x_min, part_x_max, part_y_min, part_y_max);
        else
            P.loadData(argv[2]);

        std::cerr << "Patch information" << std::endl;
        std::cerr << "V = " << P.number_of_vertices() << std::endl;
        std::cerr << "E = " << P.number_of_edges() << std::endl;
        std::cerr << "F = " << P.number_of_faces() << std::endl;

  //      std::cerr << "c to continue: ";
  //      std::cin >> ch;
		//if (ch != 'c') throw cpp::Exception("Wrong command.");



        // ------------- TSS setting ---------------

        std::cerr << "Subdividing Translation Space..." << std::endl;

        // TSS construction
        TranslationSpaceSubdivision tss(&T, &P); // with construction

        std::cerr << "Translation Space Subdivision information" << std::endl;
        std::cerr << "V = " << tss.arr.number_of_vertices() << std::endl;
		std::cerr << "E = " << tss.arr.number_of_edges() << std::endl;
		std::cerr << "F = " << tss.arr.number_of_faces() << std::endl;



        // ------------- Traverse structure ---------------

        std::cerr << "Traversal start (s to start) : ";
        std::cin >> ch;
		if (ch != 's') throw cpp::Exception("Wrong command.");

        unsigned num = 0;
        BasisResult global_optimal(Point(0,0,0),std::numeric_limits<double>::infinity()),
                    local_optimal(Point(0,0,0),std::numeric_limits<double>::infinity());
        bool first = true;

        //
        do {
            //std::cin >> ch;
            ++num;
            if (num % 1 == 0) std::cerr << "\n== Cell " << num << " ==" << std::endl;
            //tss.printCSPlanes();
            //std::cerr << ".........." << std::endl;

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

    catch (std::exception& e)
    {
        std::cerr << "An exception occurred: " << e.what() << std::endl;

		std::cin.sync();
		std::cin.get();
        return 1;
    }

}
