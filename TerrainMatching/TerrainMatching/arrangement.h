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
class Event;
class SweepLine;

///

// information of vertex
class ArrangementVertexData
{
public:
    double x,y; // coordinates
	std::multiset<Event, std::greater<Event>>::iterator it_eventQueue;

	ArrangementVertexData() { x = 0; y = 0; it_eventQueue._Ptr = NULL; }
	ArrangementVertexData(Terrain::VertexData &tvd) { x = tvd.p.x; y = tvd.p.y; }

	inline bool operator< (const ArrangementVertexData &vd) const {
		return x < vd.x || (x == vd.x && y < vd.y); // lexicographical ordering
	}
	inline bool operator== (const ArrangementVertexData &vd) const {
		return x == vd.x && y == vd.y;
	}

	inline void print(std::ostream &os) {
		os << '(' << x << ',' << y << ')';
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
	enum TraverseState { NOTPASSED, CHECKED, PASSED };

    TraverseState state; // information for DFS

    ArrangementFaceData() : state(NOTPASSED) {}
};

class ArrangementEdgeDataCompare
{
public:
	ArrangementEdgeDataCompare() {}
	inline bool operator()(const ArrangementEdgeData *ed1, const ArrangementEdgeData *ed2) const; // ed1 > ed2
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

	std::vector<Source> sources;
	HalfEdgeT<ArrangementVertexData, ArrangementHalfEdgeData, ArrangementFaceData> 
		*halfEdge_up, *halfEdge_down;
	std::multiset<ArrangementEdgeData *, ArrangementEdgeDataCompare>::iterator it_edgeDataBBT;

	ArrangementEdgeData() { halfEdge_up = NULL; halfEdge_down = NULL; it_edgeDataBBT._Ptr = NULL; }

	inline void print(std::ostream &os) {
		halfEdge_up->getOrigin()->getData().print(os);
		os << " -> ";
		halfEdge_down->getOrigin()->getData().print(os);
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
	
	inline unsigned int number_of_vertices() { return vertices.size() - erasedVerticesIndices.size(); }
	inline unsigned int number_of_halfedges() { return halfEdges.size() - erasedHalfEdgesIndices.size(); }
	inline unsigned int number_of_edges() { return edgeDataContainer.size() - erasedEdgeDataContainerIndices.size(); }
	inline unsigned int number_of_faces() { return faces.size() - erasedFacesIndices.size(); }

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
			auto idit = erasedVerticesIndices.begin();
			unsigned int id = *idit;
			erasedVerticesIndices.erase(idit);
			vertices[id] = Vertex(); // initialize
			return &vertices[id];
		}
	}
	inline HalfEdge* createHalfEdge()
	{
		if (erasedHalfEdgesIndices.empty()) {
			if (halfEdges.capacity() == halfEdges.size()) {
				throw cpp::Exception("HalfEdges capacity over.");
				return NULL;
			}
			else {
				halfEdges.push_back(HalfEdge());
				return &halfEdges.back();
			}
		}
		else {
			auto idit = erasedHalfEdgesIndices.begin();
			unsigned int id = *idit;
			erasedHalfEdgesIndices.erase(idit);
			halfEdges[id] = HalfEdge(); // initialize
			return &halfEdges[id];
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
			auto idit = erasedFacesIndices.begin();
			unsigned int id = *idit;
			erasedFacesIndices.erase(idit);
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
			auto idit = erasedEdgeDataContainerIndices.begin();
			unsigned int id = *idit;
			erasedEdgeDataContainerIndices.erase(idit);
			edgeDataContainer[id] = EdgeData(); // initialize
			return &edgeDataContainer[id];
		}
	}

	inline void deleteVertex(unsigned int id)
	{
		if (id >= vertices.size())
			throw cpp::Exception("Delete vertices-index exceeds the size.");
		else
			erasedVerticesIndices.insert(id);
	}
	inline void deleteHalfEdge(unsigned int id)
	{
		if (id >= halfEdges.size())
			throw cpp::Exception("Delete halfEdges-index exceeds the size.");
		else
			erasedHalfEdgesIndices.insert(id);
	}
	inline void deleteFace(unsigned int id)
	{
		if (id >= faces.size())
			throw cpp::Exception("Delete faces-index exceeds the size.");
		else
			erasedFacesIndices.insert(id);
	}
	inline void deleteEdgeData(unsigned int id)
	{
		if (id >= edgeDataContainer.size())
			throw cpp::Exception("Delete edgeDataContainer-index exceeds the size.");
		else
			erasedEdgeDataContainerIndices.insert(id);
	}

	// Get the outerface.
	inline Face* getOuterface() {
		return &faces.at(outerface);
	}
	// Set the outerface. (Set firstHalfEdge first.)
	inline void setOuterface(Face *f) {
		outerface = f - &faces[0];
#ifdef _DEBUG
		if (&faces[outerface] != halfEdges[firstHalfEdge].getFace())
			throw cpp::Exception("Outerface should be an adjacent face of firstHalfEdge.");
#endif
	}
	// Get a halfedge whose adjacent face is outerface.
	inline HalfEdge* getFirstHalfEdge() {
		return &halfEdges.at(firstHalfEdge);
	}
	// Set the halfedge whose adjacent face is outerface.
	inline void setFirstHalfEdge(HalfEdge *he) {
		firstHalfEdge = he - &halfEdges[0];
	}

	std::vector<Vertex> vertices;
	std::vector<HalfEdge> halfEdges;
	std::vector<Face> faces;
	std::vector<EdgeData> edgeDataContainer;

	std::set<unsigned int> erasedVerticesIndices;
	std::set<unsigned int> erasedHalfEdgesIndices;
	std::set<unsigned int> erasedFacesIndices;
	std::set<unsigned int> erasedEdgeDataContainerIndices;

	double x_min, x_max, y_min, y_max;
	double edgeLength_max;

private:
	unsigned int firstHalfEdge;
	unsigned int outerface;
	static SweepLine sweepLine;
};

// Sweepline event point.
class Event
{
public:
	Arrangement::Vertex *v;
	double x, y;

	Event() {}
	Event(Arrangement::Vertex *_v) {
		v = _v;
		x = v->getData().x;
		y = v->getData().y;
	}

	bool operator<(const Event &ep) const {
		return x < ep.x || (x == ep.x && y < ep.y); // lexicographical ordering
	}
};

// Sweep from left to right.
// BBT is sorted downward.
// Make sure that the parent arrangement has enough reserved space (at least order of nm^2).
class SweepLine
{
public:
	typedef std::multiset<Event, std::less<Event>> EventQueue;
	typedef std::multiset<Event, std::less<Event>>::iterator EventQueueIterator;
	typedef std::multiset<Arrangement::EdgeData *, ArrangementEdgeDataCompare> EdgeDataBBT;
	typedef std::multiset<Arrangement::EdgeData *, ArrangementEdgeDataCompare>::iterator EdgeDataBBTIterator;

private:
	static Arrangement *parent;
	static EventQueue events;
	static Event *currentEvent;
	static EdgeDataBBT edgeDataBBT;
	static int eventCount;
	static bool firstEvent;
	static double x_stepSize, y_stepSize;

	static EventQueueIterator events_insert(Arrangement::Vertex *v) {
		EventQueueIterator it = events.insert(Event(v));
		v->getData().it_eventQueue = it;
		return it;
	}
	static EventQueueIterator events_erase(Arrangement::Vertex *v) {
		EventQueueIterator it;
		if (v->getData().it_eventQueue._Ptr != NULL) {
			it = events.erase(v->getData().it_eventQueue);
			v->getData().it_eventQueue._Ptr = NULL;
		}
		return it;
	}
	static Event events_popfront() {
		Event e = *events.begin();
		events.erase(events.begin());
		e.v->getData().it_eventQueue._Ptr = NULL;
		return e;
	}

	static EdgeDataBBTIterator edgeDataBBT_insert(Arrangement::EdgeData *ed) {
		EdgeDataBBTIterator it = edgeDataBBT.insert(ed);
		ed->it_edgeDataBBT = it;
		return it;
	}
	static EdgeDataBBTIterator edgeDataBBT_insert(EdgeDataBBTIterator hintit, Arrangement::EdgeData *ed) {
		EdgeDataBBTIterator it = edgeDataBBT.insert(hintit, ed);
		ed->it_edgeDataBBT = it;
		return it;
	}
	static EdgeDataBBTIterator edgeDataBBT_erase(Arrangement::EdgeData *ed) {
		return edgeDataBBT.erase(ed->it_edgeDataBBT);
		ed->it_edgeDataBBT._Ptr = NULL;
	}

	static bool handleProperIntersectionEvent(Arrangement::EdgeData *ed1, Arrangement::EdgeData *ed2);
	static Arrangement::Vertex* updateDCELProperIntersection(Arrangement::EdgeData *ed1, Arrangement::EdgeData *ed2, double int_x, double int_y);
	static Arrangement::EdgeData* updateDCELVertexEdgeIntersection(Arrangement::Vertex *v, Arrangement::EdgeData *ed);
	static void updateDCEL2VertexIntersection(Arrangement::Vertex *v, Arrangement::Vertex *v_erase);
	static void updateDCELTwinEdgeWithOneSharedVertex(Arrangement::HalfEdge *he_prev, Arrangement::HalfEdge *he_next);
	static void updateDCELTwinEdgeWithTwoSharedVertex(Arrangement::HalfEdge *he_prev, Arrangement::HalfEdge *he_next);

public:
	SweepLine() {}

	static inline void initialize(Arrangement *_parent);

	static inline double getX() { return currentEvent->x; }

	static void advance();
	static inline void run()
	{
		while (!events.empty()) advance();
	}
};

#endif // ARRANGEMENT_H
