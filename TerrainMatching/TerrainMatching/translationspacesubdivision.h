#ifndef TRANSLATIONSPACESUBDIVISION_H
#define TRANSLATIONSPACESUBDIVISION_H

#include <cstddef>
#include <iostream>
#include <stack>

#include "arrangement.h"
#include "combinatorialstructure.h"
#include "terrain.h"


class TranslationSpaceSubdivision : public Arrangement
{
    // members
    std::stack<HalfEdge*> m_he_stack; // stack for DFS
    std::stack<HalfEdge*> m_he_path; // path from the root
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
    void update_CS(HalfEdge* he);
    void update_CS(EdgeData::Source &ex_data);

public: 
    enum DFSState { DFS_ADVANCED, DFS_END, DFS_FAILED };

    TranslationSpaceSubdivision(Terrain *_t1, Terrain *_t2);

    Terrain *t1; // domain
    Terrain *t2; // patch

    void init();

    DFSState advance();

    bool solveLP(BasisResult &b);

    int getInsideNum() {return m_CS.inside_num;}
    void printCSPlanes() {m_CS.printPlanes();}
};

#endif // TRANSLATIONSPACESUBDIVISION_H
