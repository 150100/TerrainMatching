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

    ArrangementVertexData() { x = 0; y = 0; }
	ArrangementVertexData(Terrain::VertexData &tvd) { x = tvd.p.x; y = tvd.p.y; }

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

		bool operator==(const Source &s) const {
			return he == s.he && v == s.v && he_is_from_patch == s.he_is_from_patch;
		}
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

	Arrangement() {}
    Arrangement(Terrain *t1, Terrain *t2);

	~Arrangement() {}
	
	inline unsigned int number_of_vertices() { return vertices.size(); }
	inline unsigned int number_of_halfedges() { return edges.size(); }
	inline unsigned int number_of_edges() { return edgeDataContainer.size(); }
	inline unsigned int number_of_faces() { return faces.size(); }
	inline Vertex* createVertex()
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
	inline HalfEdge* createHalfEdge()
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
	inline Face* createFace()
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
	inline EdgeData* createEdgeData()
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

	inline void deleteVertex(unsigned int id)
	{
		if (id >= vertices.size())
			throw cpp::Exception("Delete vertices-index exceeds the size.");
		else
			erasedVerticesIndices.push(id);
	}
	inline void deleteHalfEdge(unsigned int id)
	{
		if (id >= edges.size())
			throw cpp::Exception("Delete edges-index exceeds the size.");
		else
			erasedEdgesIndices.push(id);
	}
	inline void deleteFace(unsigned int id)
	{
		if (id >= faces.size())
			throw cpp::Exception("Delete faces-index exceeds the size.");
		else
			erasedFacesIndices.push(id);
	}
	inline void deleteEdgeData(unsigned int id)
	{
		if (id >= edgeDataContainer.size())
			throw cpp::Exception("Delete edgeDataContainer-index exceeds the size.");
		else
			erasedEdgeDataContainerIndices.push(id);
	}

	void getOverlay(Arrangement *overlay) {
		sweepLine.initialize(this, overlay);
		sweepLine.getOverlay();
	}

	std::vector<Vertex> vertices;
	std::vector<HalfEdge> edges;
	std::vector<Face> faces;
	std::vector<EdgeData> edgeDataContainer;

protected:
	std::queue<unsigned int> erasedVerticesIndices;
	std::queue<unsigned int> erasedEdgesIndices;
	std::queue<unsigned int> erasedFacesIndices;
	std::queue<unsigned int> erasedEdgeDataContainerIndices;
	
private:
	static void eraseZeroLengthEdge(EdgeData *ed);

protected:
	static SweepLine sweepLine;
};

// Sweep from left to right.
// BBT is sorted downward.
class SweepLine
{
private:
	class EventPoint
	{
	public:
		static enum EventState { VERTEXPOINT, CROSSINGPOINT };
		Arrangement::Vertex *v;
		Arrangement::EdgeData *ed1, *ed2;
		double x, y;
		EventState state;

		EventPoint(Arrangement::Vertex *_v) { // vertex point
			v = _v;
			ed1 = NULL;
			ed2 = NULL;
			x = v->getData().x;
			y = v->getData().y;
			state = VERTEXPOINT;
		}
		EventPoint(Arrangement::EdgeData *_ed1, Arrangement::EdgeData *_ed2, double _x, double _y) { // crossing point
			x = NULL;
			ed1 = _ed1;
			ed2 = _ed2;
			x = _x;
			y = _y;
			state = CROSSINGPOINT;
		}

		bool operator>(const EventPoint &ep) const {
			return x > ep.x || (x == ep.x && y > ep.y); // lexicographical ordering
		}
	};

	class EdgeDataCompare
	{
	public:
		EdgeDataCompare() {}
		inline bool operator()(const Arrangement::EdgeData *ed1, const Arrangement::EdgeData *ed2) const; // ed1 < ed2
	};

public:
	typedef std::priority_queue<EventPoint, std::vector<EventPoint>, std::greater<EventPoint>> EventQueue;
	typedef std::set<Arrangement::EdgeData *, EdgeDataCompare> EdgeDataBBT;
	typedef std::set<Arrangement::EdgeData *, EdgeDataCompare>::iterator EdgeDataBBTIterator;

private:
	static Arrangement *parent;
	static Arrangement *overlay;
	static EventQueue events;
	static EdgeDataBBT edgeDataBBT;
	static int eventCount;

	static bool handleIntersectionEvent(Arrangement::EdgeData *ed1, Arrangement::EdgeData *ed2);
	static std::pair<Arrangement::EdgeData *, Arrangement::EdgeData *> updateIntersectionDCEL(Arrangement::EdgeData *ed1, Arrangement::EdgeData *ed2, double int_x, double int_y);

public:
	SweepLine() {}

	static inline void initialize(Arrangement *_parent, Arrangement *_overlay);

	static inline double getX() { return events.top().x; }

	static void advance();
	static inline void getOverlay()
	{
		while (!events.empty()) advance();
	}
};

#endif // ARRANGEMENT_H
