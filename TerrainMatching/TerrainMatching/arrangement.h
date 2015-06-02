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

#include "common.h"
#include "numbertype.h"
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
	typedef Terrain::TerrainMesh::Vertex TerrainVertex;

    double x,y; // coordinates
	std::multiset<Event, std::greater<Event>>::iterator it_eventQueue;
	bool insideWindow;

	ArrangementVertexData() { x = 0; y = 0; it_eventQueue._Ptr = NULL; }
	ArrangementVertexData(Terrain::VertexData &tvd) { x = tvd.p.x; y = tvd.p.y; }

	inline bool operator< (const ArrangementVertexData &vd) const { // two points should not be "very close".
		const double eps = 0.0000000000001;
		std::pair<double, double> vd_x_eps = std::minmax(vd.x * (1 - eps), vd.x * (1 + eps));
		double vd_y_epsmin = std::min(vd.y * (1 - eps), vd.y * (1 + eps));
		return x < vd_x_eps.first || (x < vd_x_eps.second && y < vd_y_epsmin); // lexicographical ordering considering floating-point error.
	}
	inline bool operator== (const ArrangementVertexData &vd) const {
		return nearlyEqual(x, vd.x) && nearlyEqual(y, vd.y);
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
	enum TraverseState { NOTPASSED, CHECKED, PASSED, OUTERFACE };

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

enum GridCellSearchState { GridCellSearch_ADVANCED, GridCellSearch_END, GridCellSearch_FAILED };

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

	// Data structure
	std::vector<Vertex> vertices;
	std::vector<HalfEdge> halfEdges;
	std::vector<Face> faces;
	std::vector<EdgeData> edgeDataContainer;

	// Erased containers
	std::set<unsigned int> erasedVerticesIndices;
	std::set<unsigned int> erasedHalfEdgesIndices;
	std::set<unsigned int> erasedFacesIndices;
	std::set<unsigned int> erasedEdgeDataContainerIndices;

	// Parent terrains
	TerrainWithGrids *t1, *t2;

	// Total arrangement information
	double x_min, x_max, y_min, y_max;
	unsigned x_gridSize, y_gridSize;
	double x_gridStepSize, y_gridStepSize;

	// Window information
	unsigned cur_x_grid, cur_y_grid;
	double cur_x_min, cur_x_max, cur_y_min, cur_y_max;

	// Constructor. Make sure the grids are made.
    Arrangement(TerrainWithGrids *_t1, TerrainWithGrids *_t2);

	// Initialize space. 
	void initializeSpace();

	// Advance to the next grid cell
	GridCellSearchState advanceGridCell();

	// Insert translated copies of t1 and -t2.
	void insertTranslatedCopies();

	// Get arrangement of existing disconnected structure
	void makeArrangement();

	// Get the number of elements
	inline unsigned int number_of_vertices() { return vertices.size() - erasedVerticesIndices.size(); }
	inline unsigned int number_of_halfedges() { return halfEdges.size() - erasedHalfEdgesIndices.size(); }
	inline unsigned int number_of_edges() { return edgeDataContainer.size() - erasedEdgeDataContainerIndices.size(); }
	inline unsigned int number_of_faces() { return faces.size() - erasedFacesIndices.size(); }

	// A point is in the window?
	inline bool isInWindow(double x, double y) {
		return nearlyInClosedRange(x, cur_x_min, cur_x_max)
			&& nearlyInClosedRange(y, cur_y_min, cur_y_max);
	}

	// Is he is the firstHalfEdge of the next grid?
	inline bool isNextFirstHalfEdge(HalfEdge *he) {
		return he == firstHalfEdge_next;
	}

	// Get the outerface.
	inline Face* getOuterface() {
		return outerface;
	}
	// Set the outerface. (Set firstHalfEdge first.)
	inline void setOuterface(Face *f) {
		outerface = f;
		f->getData().state = FaceData::OUTERFACE;
#ifdef _DEBUG
		if (outerface != firstHalfEdge->getFace())
			throw cpp::Exception("Outerface should be an adjacent face of firstHalfEdge.");
#endif
	}
	// Get a halfedge whose adjacent face is outerface.
	inline HalfEdge* getFirstHalfEdge() {
		return firstHalfEdge;
	}
	// Set the halfedge whose adjacent face is outerface.
	inline void setFirstHalfEdge(HalfEdge *he) {
		firstHalfEdge = he;
	}
	// Test emptiness
	inline bool isEmpty() { return vertices.empty(); }

	// Maintaining data structure
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

private:
	HalfEdge *firstHalfEdge;
	HalfEdge *firstHalfEdge_next; // This information is used to remember an initial CS for the next cell.
	Face *outerface;
	static SweepLine sweepLine;

};

// Sweepline event point.
class Event
{
public:
	typedef Arrangement::Vertex ArrangementVertex;

	ArrangementVertex *v;

	Event() {}
	Event(ArrangementVertex *_v) {
		v = _v;
	}

	bool operator<(const Event &ep) const {
		return v->getData() < ep.v->getData();
	}
};

// Sweep from left to right.
// BBT is sorted downward.
// Make sure that the parent arrangement has enough reserved space (at least order of nm^2).
class SweepLine
{
public:
	typedef Arrangement::Vertex ArrangementVertex;
	typedef Arrangement::HalfEdge ArrangementHalfEdge;
	typedef Arrangement::Face ArrangementFace;
	typedef Arrangement::EdgeData ArrangementEdgeData;

	typedef std::multiset<Event, std::less<Event>> EventQueue;
	typedef std::multiset<Event, std::less<Event>>::iterator EventQueueIterator;
	typedef std::multiset<ArrangementEdgeData *, ArrangementEdgeDataCompare> EdgeDataBBT;
	typedef std::multiset<ArrangementEdgeData *, ArrangementEdgeDataCompare>::iterator EdgeDataBBTIterator;

private:
	static Arrangement *parent;
	static EventQueue events;
	static Event *currentEvent;
	static EdgeDataBBT edgeDataBBT;
	static int eventCount;
	static bool firstEvent;
	static double x_stepSize, y_stepSize;

	static EventQueueIterator events_insert(ArrangementVertex *v) {
		//std::cerr << "event insert : " << v << '\n';
		EventQueueIterator it = events.emplace(v);
		v->getData().it_eventQueue = it;
		return it;
	}
	static EventQueueIterator events_erase(ArrangementVertex *v) {
		//std::cerr << "event erase : " << v << '\n';
		EventQueueIterator it;
		if (v->getData().it_eventQueue._Ptr != NULL) {
			it = events.erase(v->getData().it_eventQueue);
			v->getData().it_eventQueue._Ptr = NULL;
		}
		return it;
	}
	static Event events_popfront() {
		Event e = *events.begin();
		//std::cerr << "event popfront : " << e.v << '\n';
		events.erase(events.begin());
		e.v->getData().it_eventQueue._Ptr = NULL;
		return e;
	}

	static EdgeDataBBTIterator edgeDataBBT_insert(ArrangementEdgeData *ed) {
		EdgeDataBBTIterator it = edgeDataBBT.insert(ed);
		ed->it_edgeDataBBT = it;
		return it;
	}
	static EdgeDataBBTIterator edgeDataBBT_insert(EdgeDataBBTIterator hintit, ArrangementEdgeData *ed) {
		EdgeDataBBTIterator it = edgeDataBBT.insert(hintit, ed);
		ed->it_edgeDataBBT = it;
		return it;
	}
	static EdgeDataBBTIterator edgeDataBBT_erase(ArrangementEdgeData *ed) {
		EdgeDataBBTIterator ret = edgeDataBBT.erase(ed->it_edgeDataBBT);
		ed->it_edgeDataBBT._Ptr = NULL;
		return ret;
	}

	static ArrangementVertex* handleIntersectionEvent(ArrangementEdgeData *ed1, ArrangementEdgeData *ed2);
	static ArrangementVertex* updateDCELIntersection(ArrangementEdgeData *ed1, ArrangementEdgeData *ed2, double int_x, double int_y);
	static ArrangementVertex* updateDCELProperIntersection(ArrangementEdgeData *ed1, ArrangementEdgeData *ed2, double int_x, double int_y);
	static ArrangementEdgeData* updateDCELVertexEdgeIntersection(ArrangementVertex *v, ArrangementEdgeData *ed, bool survive_v);
	static void updateDCEL2VertexIntersection(ArrangementVertex *v, ArrangementVertex *v_erase);
	static void updateDCELTwinEdgeWithOneSharedVertex(ArrangementHalfEdge *he_prev, ArrangementHalfEdge *he_next);
	static void updateDCELTwinEdgeWithTwoSharedVertex(ArrangementHalfEdge *he_prev, ArrangementHalfEdge *he_next);

public:
	SweepLine() {}

	static inline void initialize(Arrangement *_parent);

	static inline double getX() { return currentEvent->v->getData().x; }

	static void advance();
	static inline void run()
	{
		while (!events.empty()) advance();
	}
};

#endif // ARRANGEMENT_H
