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
		std::cerr << std::setprecision(6) << '(' << x << ',' << y << ')';
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

// TO DO : This should be singleton.
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
	
	static inline unsigned int number_of_vertices() { return vertices.size(); }
	static inline unsigned int number_of_halfedges() { return edges.size(); }
	static inline unsigned int number_of_edges() { return edgeDataContainer.size(); }
	static inline unsigned int number_of_faces() { return faces.size(); }

protected:
	static inline Vertex* createVertex()
	{
		if (erasedVerticesIndices.empty()) {
			if (vertices.capacity() == vertices.size()) {
				throw cpp::Exception("Vertices capacity over.");
				return NULL;
			}
			else {
				vertices.push_back(Vertex());
				return &vertices.back();
			}
		}
		else {
			unsigned int id = erasedVerticesIndices.front();
			erasedVerticesIndices.pop();
			vertices[id] = Vertex(); // initialize
			return &vertices[id];
		}
	}
	static inline HalfEdge* createHalfEdge()
	{
		if (erasedEdgesIndices.empty()) {
			if (edges.capacity() == edges.size()) {
				throw cpp::Exception("Edges capacity over.");
				return NULL;
			}
			else {
				edges.push_back(HalfEdge());
				return &edges.back();
			}
		}
		else {
			unsigned int id = erasedEdgesIndices.front();
			erasedEdgesIndices.pop();
			edges[id] = HalfEdge(); // initialize
			return &edges[id];
		}
	}
	static inline Face* createFace()
	{
		if (erasedFacesIndices.empty()) {
			if (faces.capacity() == faces.size()) {
				throw cpp::Exception("Faces capacity over.");
				return NULL;
			}
			else {
				faces.push_back(Face());
				return &faces.back();
			}
		}
		else {
			unsigned int id = erasedFacesIndices.front();
			erasedFacesIndices.pop();
			faces[id] = Face(); // initialize
			return &faces[id];
		}
	}
	static inline EdgeData* createEdgeData()
	{
		if (erasedEdgeDataContainerIndices.empty()) {
			if (edgeDataContainer.capacity() == edgeDataContainer.size()) {
				throw cpp::Exception("EdgeDataContainer capacity over.");
				return NULL;
			}
			else {
				edgeDataContainer.push_back(EdgeData());
				return &edgeDataContainer.back();
			}
		}
		else {
			unsigned int id = erasedEdgeDataContainerIndices.front();
			erasedEdgeDataContainerIndices.pop();
			edgeDataContainer[id] = EdgeData(); // initialize
			return &edgeDataContainer[id];
		}
	}

	static inline void deleteVertex(unsigned int id) 
	{
		if (id >= vertices.size())
			throw cpp::Exception("Delete vertices-index exceeds the size.");
		else
			erasedVerticesIndices.push(id);
	}
	static inline void deleteHalfEdge(unsigned int id)
	{
		if (id >= edges.size())
			throw cpp::Exception("Delete edges-index exceeds the size.");
		else
			erasedEdgesIndices.push(id);
	}
	static inline void deleteFace(unsigned int id)
	{
		if (id >= faces.size())
			throw cpp::Exception("Delete faces-index exceeds the size.");
		else
			erasedFacesIndices.push(id);
	}
	static inline void deleteEdgeData(unsigned int id)
	{
		if (id >= edgeDataContainer.size())
			throw cpp::Exception("Delete edgeDataContainer-index exceeds the size.");
		else
			erasedEdgeDataContainerIndices.push(id);
	}

    static std::vector<Vertex> vertices;
	static std::vector<HalfEdge> edges;
	static std::vector<Face> faces;
	static std::vector<EdgeData> edgeDataContainer;

	static std::queue<unsigned int> erasedVerticesIndices;
	static std::queue<unsigned int> erasedEdgesIndices;
	static std::queue<unsigned int> erasedFacesIndices;
	static std::queue<unsigned int> erasedEdgeDataContainerIndices;

private:

	// Sweep from left to right.
	// BBT is sorted downward.
	class SweepLine
	{
	private:
		class EventPoint
		{
		public:
			static enum EventState { STARTPOINT, CROSSING, ENDPOINT };
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
				return x > ep.x || (x == ep.x && y > ep.y) 
					|| (x == ep.x && y == ep.y && state < ep.state); // handle endpoint first, crossing next, startpoint last (to reduce the number of BBT)
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
