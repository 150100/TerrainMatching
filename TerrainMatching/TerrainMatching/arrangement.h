#ifndef ARRANGEMENT_H
#define ARRANGEMENT_H

#include "DCEL\Mesh.h"
#include <fstream>
#include <list>
#include <queue>
#include <set>
#include <map>
#include <functional>
#include <iomanip>

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

	inline bool operator< (const ArrangementVertexData &vd) const {
		return x < vd.x || (x == vd.x && y < vd.y); // lexicographical ordering
	}
	inline bool operator== (const ArrangementVertexData &vd) const {
		return x == vd.x && y == vd.y;
	}

	inline void print() {
		std::cerr << '(' << x << ',' << y << ')';
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
	
	inline void print() {
		halfEdge_up->getOrigin()->getData().print();
		std::cerr << " -> ";
		halfEdge_down->getOrigin()->getData().print();
	}
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

	static inline unsigned int number_of_vertices() {return vertices.size();}
	static inline unsigned int number_of_halfedges() { return edges.size(); }
	static inline unsigned int number_of_edges() { return edgeDataContainer.size(); }
	static inline unsigned int number_of_faces() { return faces.size(); }

protected:
    static std::vector<Vertex> vertices;
	static std::vector<HalfEdge> edges;
	static std::vector<Face> faces;
	static std::vector<EdgeData> edgeDataContainer;

private:

	// Sweep from left to right.
	// BBT is sorted downward.
	class SweepLine
	{
	private:
		class EventPoint
		{
		public:
			static enum EventState { STARTPOINT, ENDPOINT, CROSSING };
			EdgeData *ed1, *ed2, *ed1N, *ed2N;
			double x, y;
			EventState state;

			EventPoint(EdgeData *_ed1, EventState _state) {
				ed1 = _ed1;
				state = _state;

				if (state == STARTPOINT) {
					x = ed1->halfEdge_up->getOrigin()->getData().x;
					y = ed1->halfEdge_up->getOrigin()->getData().y;
				}
				else if (state == ENDPOINT) {
					x = ed1->halfEdge_down->getOrigin()->getData().x;
					y = ed1->halfEdge_down->getOrigin()->getData().y;
				}
				else throw cpp::Exception("Invalid creation of ep. (state invalid)");
			}
			EventPoint(EdgeData *_ed1, EdgeData *_ed2, EdgeData *_ed1N, EdgeData *_ed2N, Vertex *int_v, EventState _state) {
				ed1 = _ed1;
				ed2 = _ed2;
				ed1N = _ed1N;
				ed2N = _ed2N;
				state = _state;

				if (state == CROSSING) {
					x = int_v->getData().x;
					y = int_v->getData().y;
				}
				else throw cpp::Exception("Invalid creation of ep. (state invalid)");
			}

			bool operator>(const EventPoint &ep) const {
				return x > ep.x || (x == ep.x && y > ep.y); 
			}
		};

		class EdgeDataCompare
		{
		public:
			EdgeDataCompare() {}
			inline bool operator()(const EdgeData *ed1, const EdgeData *ed2) const; // ed1 < ed2
		};

	public:
		typedef std::priority_queue<EventPoint, std::vector<EventPoint>, std::greater<EventPoint>> EventQueue;
		typedef std::multiset<EdgeData *, EdgeDataCompare> EdgeDataBBT;
		typedef std::multiset<EdgeData *, EdgeDataCompare>::iterator EdgeDataBBTIterator;

	private:
		//Arrangement *parent;
		static EventQueue events;
		static EdgeDataBBT edgeDataBBT;
		static int eventCount;

		static bool handleIntersectionEventWithDCEL(EdgeData *ed1, EdgeData *ed2);

	public:
		SweepLine() {}
		
		static inline void initialize();

		static inline double getX() { return events.top().x; }

		static void advance();
		static inline void run() { while (!events.empty()) advance(); };
	};

protected:
	static SweepLine sweepLine;
};

#endif // ARRANGEMENT_H
