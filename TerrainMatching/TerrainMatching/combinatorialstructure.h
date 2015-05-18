#ifndef COMBINATORIALSTRUCTURE_H
#define COMBINATORIALSTRUCTURE_H

#include "terrain.h"
#include "matousek.h"
#include <list>

///

class CombinatorialPair
{
protected:
    bool m_inserted;
    std::list<Plane>::iterator it_plane; // kind of pointer to a list element

public:
    CombinatorialPair() {m_inserted = false;}
    virtual ~CombinatorialPair() {}

    /**
     * @brief makePlane Make a plane from the combinatorial pair.
     * @return 3D plane object.
     */
    virtual Plane makePlane() = 0;

    bool is_inserted() {return m_inserted;}

    void insertPlane(std::list<Plane> &plane_list) {
        if (m_inserted) {
            std::cout << "Error: already inserted." << std::endl;
        }
        else {
            plane_list.push_front(makePlane());
            m_inserted = true;
            it_plane = plane_list.begin();
        }
    }
    void removePlane(std::list<Plane> &plane_list) {
        if (!m_inserted) {
            std::cout << "Error: removing non-inserted data." << std::endl;
        }
        else {
            plane_list.erase(it_plane);
            m_inserted = false;
        }
    }

    /**
     * @brief print Print pair information.
     */
    virtual void print() = 0;
};

class EdgeEdgePair : public CombinatorialPair
{
public:
    typedef Terrain::TerrainMesh::Vertex TerrainVertex;
    typedef Terrain::TerrainMesh::HalfEdge TerrainHalfEdge;
	typedef Terrain::EdgeData TerrainEdgeData;

    class Edge
    {
    public:
        Edge() {}
        Edge(Edge &_e) {*this = _e;}
        Edge(TerrainHalfEdge *_he) {he = _he;}

        TerrainHalfEdge *he;

        bool sameEdge(Edge &_e) {
            return  (_e.he->getOrigin() == he->getOrigin()
                     && _e.he->getNext()->getOrigin() == he->getNext()->getOrigin())
                 || (_e.he->getOrigin() == he->getNext()->getOrigin()
                     && _e.he->getNext()->getOrigin() == he->getOrigin());
        }
        bool sameEdge(TerrainHalfEdge *_he) {
            return  (_he->getOrigin() == he->getOrigin()
                     && _he->getNext()->getOrigin() == he->getNext()->getOrigin())
                 || (_he->getOrigin() == he->getNext()->getOrigin()
                     && _he->getNext()->getOrigin() == he->getOrigin());
        }
        bool containsVertex(TerrainVertex *_v) {
            return he->getOrigin() == _v || he->getNext()->getOrigin() == _v;
        }
    };

    Edge edge1, edge2;
    bool edge1_is_from_patch;

    EdgeEdgePair(TerrainHalfEdge *_he1, TerrainHalfEdge *_he2, bool _he1_is_from_patch)
        : edge1(_he1), edge2(_he2), edge1_is_from_patch(_he1_is_from_patch) {}

    virtual Plane makePlane();

    virtual void print() {
        std::cout << "edge-edge pair, ck: " << edge1_is_from_patch;
    }
};

class VertexTrianglePair : public CombinatorialPair
{
public:
    typedef Terrain::TerrainMesh::Vertex Vertex;
    typedef Terrain::TerrainMesh::HalfEdge HalfEdge;
    typedef Terrain::TerrainMesh::Face Triangle;

    Vertex* pVertex;
    Triangle* pTriangle;
    bool tri_is_from_patch;

    VertexTrianglePair(Vertex* _v, Triangle* _t, bool _tri_is_from_patch)
        : pVertex(_v), pTriangle(_t), tri_is_from_patch(_tri_is_from_patch) {}

    virtual Plane makePlane();

    virtual void print() {
        std::cout << "vertex-triangle pair, ck: " << tri_is_from_patch;
    }
};

class CombinatorialStructure
{
public:
    unsigned int inside_num;

    typedef Terrain::TerrainMesh::Vertex TerrainVertex;
    typedef Terrain::TerrainMesh::Face TerrainFace;
    typedef Terrain::TerrainMesh::HalfEdge TerrainHalfEdge;
    typedef Terrain::EdgeData TerrainEdge;

    CombinatorialStructure() {inside_num = 0;}

    std::list<Plane> m_plane_list;

    VertexTrianglePair *insertPair(TerrainVertex* v, TerrainFace* tri, bool tri_is_from_patch) {

        VertexTrianglePair *p_pair = new VertexTrianglePair(v, tri, tri_is_from_patch);
        p_pair->insertPlane(m_plane_list);

        std::cout << "CS::insertPair(vh,fh,ck). ";
        p_pair->print();
        std::cout << std::endl;

        return p_pair;
    }
    EdgeEdgePair *insertPair(TerrainHalfEdge *he1, TerrainHalfEdge *he2, bool he1_is_from_patch) {

        EdgeEdgePair *p_pair = new EdgeEdgePair(he1, he2, he1_is_from_patch);
        p_pair->insertPlane(m_plane_list);

        std::cout << "CS::insertPair(e1,e2,ck). ";
        p_pair->print();
        std::cout << std::endl;

        return p_pair;
    }

    void removePair(CombinatorialPair *p_pair) {

        p_pair->removePlane(m_plane_list);

        std::cout << "CS::removePair(p_pair). ";
        p_pair->print();
        std::cout << std::endl;

        delete p_pair;
    }

    BasisResult solveLP(std::vector<Point> boundary_pts); // boundary_pts = counter-clockwise boundary point set

    void printPlanes();
};

#endif // COMBINATORIALSTRUCTURE_H
