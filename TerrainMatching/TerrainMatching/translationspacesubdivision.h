#ifndef TRANSLATIONSPACESUBDIVISION_H
#define TRANSLATIONSPACESUBDIVISION_H

#include <cstddef>
#include <iostream>
#include <stack>

#include "arrangement.h"
#include "combinatorialstructure.h"
#include "terrain.h"

enum DFSState { DFS_ADVANCED, DFS_END, DFS_FAILED };

class TranslationSpaceSubdivision
{
	typedef Arrangement::TerrainHalfEdge TerrainHalfEdge;
	typedef Arrangement::TerrainVertex TerrainVertex;
	typedef Arrangement::TerrainFace TerrainFace;

    // members
	std::stack<Arrangement::HalfEdge*> m_he_stack; // stack for DFS
	std::stack<Arrangement::HalfEdge*> m_he_path; // path from the root
    CombinatorialStructure m_CS;
    std::vector<bool *> m_checked_plist; // pointers to checked

    // private function
    /**
     * @brief switch_VTpair Switch a vertex-triangle pair to new one according to the movement.
     * @param he A halfedge that v moves over.
     * @param v A vertex that crosses he.
     * @param he_is_from_patch Whether he is from patch.
     */
    void switch_VTpair(TerrainHalfEdge *he, TerrainVertex *v, bool he_is_from_patch);
    void switch_EEpairs(TerrainHalfEdge *he, TerrainVertex *v, bool he_is_from_patch);
//    void remove_crossing_EEpairs(TerrainHalfEdge &he, TerrainVertex &v, bool he_is_from_patch);
//    void insert_EEpairs_without_checked(TerrainHalfEdge &he, TerrainVertex &v, bool he_is_from_patch);
	void update_CS(Arrangement::HalfEdge* he);
	void update_CS(Arrangement::EdgeData::Source &ex_data);

public:
	// Constructor. Make sure that t1 and t2 are sored.
	TranslationSpaceSubdivision(TerrainWithGrids *_t1, TerrainWithGrids *_t2);

	// arrangement
	Arrangement arr; 

	// domain
	TerrainWithGrids *t1;

	// patch
	TerrainWithGrids *t2;

	// Initialization function. Call whenever a new arrangement structure is made.
    void init();

	// Go to next cell through DFS.
    DFSState advanceDFS();

	// Go to next grid-cell of the arrangement. It makes a new arrangement structure.
	inline GridCellSearchState advanceGridCell() { 
		GridCellSearchState state = arr.advanceGridCell(); 
		init(); 
		return state; 
	}

	// Solve LP-type problem to get an optimal solution under the current arrangement.
    bool solveLP(BasisResult &b);

    inline int getInsideNum() {return m_CS.inside_num;}
    inline void printCSPlanes() { m_CS.printPlanes(); }
	inline int numCSPlanes() { return m_CS.m_plane_list.size(); }
};

#endif // TRANSLATIONSPACESUBDIVISION_H
