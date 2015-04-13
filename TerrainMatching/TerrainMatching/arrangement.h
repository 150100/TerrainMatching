#ifndef ARRANGEMENT_H
#define ARRANGEMENT_H

#include "DCEL/Mesh.h"
#include <fstream>
#include <list>
#include <queue>
#include <set>
#include <functional>

#include "terrain.h"

///

class ArrangementEdgeData;

///

// information of vertex
class ArrangementVertexData
{
public:
    double x,y; // coordinates

    ArrangementVertexData() {x=0; y=0;}
	ArrangementVertexData(Terrain::VertexData &tvd) {x=tvd.p.x; y=tvd.p.y;}

	bool operator< (const ArrangementVertexData &vd) const {
		return x < vd.x || (x == vd.x && y < vd.y); // lexicographical ordering
	}
	bool operator== (const ArrangementVertexData &vd) const {
		return x == vd.x && y == vd.y;
	}
};

// information of halfedge

class ArrangementHalfEdgeData
{
public:
	ArrangementEdgeData *edgeData;

    ArrangementHalfEdgeData() : edgeData(NULL) {}
};

// information of face
class ArrangementFaceData
{
public:
    enum TraverseState { PASSED, NOTPASSED }; // outer face is NULL

    TraverseState state; // information for DFS

    ArrangementFaceData() : state(NOTPASSED) {}
};

// information of edge
class ArrangementEdgeData
{
public:
	typedef Terrain::TerrainMesh::Vertex TerrainVertex;
	typedef Terrain::TerrainMesh::HalfEdge TerrainHalfEdge;

	class Source {
	public:
		TerrainHalfEdge *he; // copied edge
		TerrainVertex *v; // reference point
		bool he_is_from_patch; // Terrain edge e is from patch?
	};

	typedef std::vector<Source>::iterator SourceIterator;

	ArrangementEdgeData() {}

	std::vector<Source> sources;
	HalfEdgeT<ArrangementVertexData, ArrangementHalfEdgeData, ArrangementFaceData> 
		*halfEdge_up, *halfEdge_down;
	
};

///

class Arrangement
{
public:
	typedef Terrain::TerrainMesh TerrainMesh;
    typedef TerrainMesh::Vertex TerrainVertex;
    typedef TerrainMesh::HalfEdge TerrainHalfEdge;
    typedef TerrainMesh::Face TerrainFace;

    typedef Terrain::VertexData TerrainVertexData;
    typedef Terrain::HalfEdgeData TerrainHalfEdgeData;
    typedef Terrain::EdgeData TerrainEdgeData;
    typedef Terrain::FaceData TerrainFaceData;

    typedef ArrangementVertexData VertexData;
    typedef ArrangementHalfEdgeData HalfEdgeData;
	typedef ArrangementFaceData FaceData;
	typedef ArrangementEdgeData EdgeData;

	typedef VertexT<VertexData, HalfEdgeData, FaceData> Vertex;
	typedef HalfEdgeT<VertexData, HalfEdgeData, FaceData> HalfEdge;
	typedef FaceT<VertexData, HalfEdgeData, FaceData> Face;
	typedef EdgeIteratorT<VertexData, HalfEdgeData, FaceData> EdgeIterator;

    Arrangement(Terrain *t1, Terrain *t2);

	~Arrangement() {}

	unsigned int number_of_vertices() {return vertices.size();}
	unsigned int number_of_halfedges() {return edges.size();}
	unsigned int number_of_edges() {return edgeDataContainer.size();}
    unsigned int number_of_faces() {return faces.size();}

protected:
    std::vector<Vertex> vertices;
    std::vector<HalfEdge> edges;
    std::vector<Face> faces;
	std::vector<EdgeData> edgeDataContainer;

private:

	// Sweep from left to right.
	// BBT is sorted downward.
	class SweepLine
	{
		class EventPoint
		{
		public:
			enum EventState { ENDPOINT, CROSSING };
			Vertex *vertex;
			double x, y;
			EventState state;

			EventPoint(Vertex *_vertex, EventState _state) {
				vertex = _vertex;
				x = _vertex->getData().x;
				y = _vertex->getData().y;
				state = _state;
			}
			bool operator>(const EventPoint &ep) const {
				return x > ep.x || (x == ep.x && y > ep.y); 
			}
		};

		struct EdgeDataCompare
		{
			bool operator()(const EdgeData *ed1, const EdgeData *ed2) const;
		};

		typedef std::priority_queue<EventPoint, std::vector<EventPoint>, std::greater<EventPoint>> EventQueue;
		typedef std::set<EdgeData *, EdgeDataCompare> EdgeDataSet;
		typedef std::set<EdgeData *, EdgeDataCompare>::iterator EdgeDataBBTIterator;
		
		Arrangement *parent;
		EventQueue events;
		EdgeDataSet edgeDataBBT; 

		bool handleIntersectionEventWithDCEL(EdgeData *ed1, EdgeData *ed2);

	public:
		SweepLine(Arrangement *_parent);

		void advance();
		void run() { while (!events.empty()) advance(); };
	};

};

#endif // ARRANGEMENT_H
