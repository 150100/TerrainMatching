#include "arrangement.h"

#include <algorithm>

#ifdef DEBUG
	std::set<Arrangement::Vertex *> debug_handledVertices;
#endif

SweepLine Arrangement::sweepLine;

Arrangement *SweepLine::parent = NULL;
SweepLine::EventQueue SweepLine::events;
Event *SweepLine::currentEvent = NULL;
SweepLine::EdgeDataBBT SweepLine::edgeDataBBT;
int SweepLine::eventCount = 0;
bool SweepLine::firstEvent = true;
double SweepLine::x_stepSize = 0;
double SweepLine::y_stepSize = 0;

template <class T>
class Vector2
{
public:
	Vector2() {}
	Vector2(double _x, double _y) : x(_x), y(_y) {}
	Vector2(const ArrangementVertexData &vd) : x(vd.x), y(vd.y) {}
	Vector2(const Arrangement::HalfEdge *he) {
		*this = Vector2<T>(he->getTwin()->getOrigin()->getData()) - Vector2<T>(he->getOrigin()->getData());
	}

	double x, y;

	// lexicographical order
	inline bool operator>(const Vector2<T> &v) const {
		return x > v.x || (x == v.x && y > v.y);
	}
	// equal position
	inline bool operator==(const Vector2<T> &v) const {
		return x == v.x && y == v.y;
	}
	// get a vector from v to *this
	inline Vector2<T> operator-(const Vector2<T> &v) const { 
		return Vector2<T>(x - v.x, y - v.y); 
	}
	// dot product
	inline double operator*(const Vector2<T> &v) const {
		return x * v.x + y * v.y; 
	}
	// determinant of *this and v (>0 if *this is left from v, <0 if *this is right from v.)
	inline double det(const Vector2<T> &v) const {
		return x * v.y - v.x * y;
	}
	// true if *this is left from v, false otherwise.
	inline bool isLeftFrom(const Vector2<T> &v) const {
		return det(v) < 0;
	}
};

typedef Vector2<double> ArrangementVector;

template <class T>
class ImplicitLine {
	T a, b, c;
public:
	ImplicitLine(T x1, T y1, T x2, T y2) {
		a = y2 - y1;
		b = x1 - x2;
		c = x2*y1 - x1*y2;
	}
	ImplicitLine(Vector2<T> &v1, Vector2<T> &v2) {
		a = v2.y - v1.y;
		b = v1.x - v2.x;
		c = v2.x*v1.y - v1.x*v2.y;
	}
	double getValue(T x, T y) {
		return a*x + b*y + c;
	}
	double getValue(Vector2<T> &v) {
		return a*v.x + b*v.y + c;
	}
};

typedef ImplicitLine<double> ArrangementImplicitLine;

bool arrVertexCompare(Arrangement::Vertex &v1, Arrangement::Vertex &v2)
{
	ArrangementVertexData &vd1 = v1.getData();
	ArrangementVertexData &vd2 = v2.getData();
	return (vd1.x < vd2.x) || (vd1.x == vd2.x && vd1.y < vd2.y);
}

///

Arrangement::Arrangement(TerrainWithGrids *_t1, TerrainWithGrids *_t2)
{
	const double WINDOW_EDGEMAX_RATIO = 1; // must be (<= 1)

	// Parent terrains
	if (!_t1->gridIsMade) throw cpp::Exception("Grid is not made on t1.");
	if (!_t2->gridIsMade) throw cpp::Exception("Grid is not made on t2.");
	t1 = _t1;
	t2 = _t2;

	// Compute total arrangement information
	x_min = t1->x_min - t2->x_max;
	x_max = t1->x_max - t2->x_min;
	y_min = t1->y_min - t2->y_max;
	y_max = t1->y_max - t2->y_min;
	const double edgeLength_max = std::max(t1->edgeLength_max, t2->edgeLength_max);
	const double windowSideLength = WINDOW_EDGEMAX_RATIO * edgeLength_max;
	x_gridStepSize = y_gridStepSize = windowSideLength;
	x_gridSize = std::ceil((x_max - x_min) / x_gridStepSize);
	y_gridSize = std::ceil((y_max - y_min) / y_gridStepSize);

	// Initialize the first window
	cur_x_grid = cur_y_grid = 0;
	cur_x_min = x_min;
	cur_x_max = x_min + x_gridStepSize;
	cur_y_min = y_min;
	cur_y_max = y_min + y_gridStepSize;

	// initialize space
	const double asymptotic_space = t1->number_of_vertices() + t2->number_of_vertices() * t2->number_of_vertices();
	vertices.reserve(10 * asymptotic_space);
	halfEdges.reserve(40 * asymptotic_space);
	edgeDataContainer.reserve(20 * asymptotic_space);

	// insert translated copies to the structure
	insertTranslatedCopies();

	// make arrangement of the inserted structure
	makeArrangement();
}

void Arrangement::initializeSpace()
{
	vertices.clear();
	halfEdges.clear();
	faces.clear();
	edgeDataContainer.clear();
	erasedVerticesIndices.clear();
	erasedHalfEdgesIndices.clear();
	erasedFacesIndices.clear();
	erasedEdgeDataContainerIndices.clear();
}

GridCellSearchState Arrangement::advanceGridCell()
{
	// initialize data structure before
	initializeSpace();

	// get grid index (zig-zag traversal)
	if (cur_y_grid % 2 == 0) {
		if (cur_x_grid == x_gridSize - 1)
			++cur_y_grid;
		else
			++cur_x_grid;
	}
	else {
		if (cur_x_grid == 0)
			++cur_y_grid;
		else
			--cur_x_grid;
	}

	if (cur_y_grid == y_gridSize)
		return GridCellSearch_END;

	// insert translated copies to the structure
	insertTranslatedCopies();

	// make arrangement of inserted structure
	makeArrangement();

	return GridCellSearch_ADVANCED;
}

void Arrangement::insertTranslatedCopies()
{
	// get grid range
	cur_x_min = cur_x_grid * x_gridStepSize + x_min;
	cur_x_max = (cur_x_grid + 1) * x_gridStepSize + x_min;
	cur_y_min = cur_y_grid * y_gridStepSize + y_min;
	cur_y_max = (cur_y_grid + 1) * y_gridStepSize + y_min;
	ArrangementVector cur_dL(cur_x_min, cur_y_min);
	ArrangementVector cur_uL(cur_x_min, cur_y_max);
	ArrangementVector cur_dR(cur_x_max, cur_y_min);
	ArrangementVector cur_uR(cur_x_max, cur_y_max);

	/* structure of the grid range */
	// collect translations of (t1) from the vertices of (-t2). ((t1 - v{t2}) ↘ grid ℅ nil ╰ v{t2} ↑ t1 - grid)
	std::vector<TerrainVertex *> transScope_t2;
	double transScope_t2_rangeX_min = t1->x_min - cur_x_max;
	double transScope_t2_rangeX_max = t1->x_max - cur_x_min;
	double transScope_t2_rangeY_min = t1->y_min - cur_y_max;
	double transScope_t2_rangeY_max = t1->y_max - cur_y_min;
	t2->appendVerticesInRange(transScope_t2_rangeX_min, transScope_t2_rangeX_max, transScope_t2_rangeY_min, transScope_t2_rangeY_max, &transScope_t2);

	// for each translation vectors, get copies inside the window.
	unsigned int verticesIdx = 0;
	unsigned int halfEdgesIdx = 0;
	unsigned int edgeDataIdx = 0;
	for (unsigned int i = 0; i < transScope_t2.size(); ++i) {
		TerrainVertex *v_t2 = transScope_t2.at(i);

		// collect structure of (t1) in such range. (t1 = grid + v_t2).
		std::vector<TerrainVertex *> vertexScope_t1;
		std::vector<TerrainHalfEdge *> halfEdgeScope_t1;
		double vertexScope_t1_rangeX_min = cur_x_min + v_t2->getData().p.x;
		double vertexScope_t1_rangeX_max = cur_x_max + v_t2->getData().p.x;
		double vertexScope_t1_rangeY_min = cur_y_min + v_t2->getData().p.y;
		double vertexScope_t1_rangeY_max = cur_y_max + v_t2->getData().p.y;
		t1->appendVerticesInRange(vertexScope_t1_rangeX_min, vertexScope_t1_rangeX_max, vertexScope_t1_rangeY_min, vertexScope_t1_rangeY_max, &vertexScope_t1);

		// vertices in the window scope
		std::vector<TerrainVertex *> vertexOutsideWindow_t1;
		unsigned int halfEdgesIdx_first = halfEdgesIdx;
		unsigned int edgeDataIdx_first = edgeDataIdx;
		for (unsigned int j = 0; j < vertexScope_t1.size(); ++j) {
			TerrainVertex *v_t1 = vertexScope_t1.at(j);
			double v_arr_x = v_t1->getData().p.x - v_t2->getData().p.x;
			double v_arr_y = v_t1->getData().p.y - v_t2->getData().p.y;

			// create vertex info only if the vertex is in the window.
			if (isInWindow(v_arr_x, v_arr_y)) {
				Vertex *v_arr = createVertex();
				v_t1->getData().arrIndex = verticesIdx++; // mark vertex index
				v_arr->getData().x = v_arr_x;
				v_arr->getData().y = v_arr_y;
				v_arr->getData().insideWindow = true;

				// mark halfedge index around the vertex
				bool first = true;
				Terrain::EdgeIterator eit(v_t1);
				while (eit.hasNext()) {
					// insert the halfedge.
					TerrainHalfEdge *he_t1 = eit.getNext();
					if (he_t1->getData().arrIndex == -1) {
						createHalfEdge(); // allocate space for the halfedge
						if (first) {
							v_arr->setIncidentEdge(&halfEdges.at(halfEdgesIdx));
							first = false;
						}
						he_t1->getData().arrIndex = halfEdgesIdx++; // mark halfEdge index
						halfEdgeScope_t1.push_back(he_t1); // collect halfedges in the window scope
					}
					else if (first) {
						v_arr->setIncidentEdge(&halfEdges.at(he_t1->getData().arrIndex));
						first = false;
					}
					// insert the twin halfedge.
					TerrainHalfEdge *he_twin_t1 = he_t1->getTwin();
					if (he_twin_t1->getData().arrIndex == -1) {
						HalfEdge *he_twin = createHalfEdge(); // allocate space for the halfedge
						he_twin_t1->getData().arrIndex = halfEdgesIdx++; // mark halfEdge index
						halfEdgeScope_t1.push_back(he_twin_t1); // collect halfedges in the window scope

						// insert a vertex just outside of the window.
						TerrainVertex *v_twin_t1 = he_twin_t1->getOrigin();
						if (v_twin_t1->getData().arrIndex == -1) {
							double v_twin_x = v_twin_t1->getData().p.x - v_t2->getData().p.x;
							double v_twin_y = v_twin_t1->getData().p.y - v_t2->getData().p.y;
							if (!isInWindow(v_twin_x, v_twin_y)) {
								Vertex *v_twin = createVertex();
								v_twin_t1->getData().arrIndex = verticesIdx++; // mark vertex index
								v_twin->getData().x = v_twin_x;
								v_twin->getData().y = v_twin_y;
								v_twin->getData().insideWindow = false;
								v_twin->setIncidentEdge(he_twin);
								vertexOutsideWindow_t1.push_back(v_twin_t1);
							}
						}
					}
				}
			}
		}

		// edges in the window scope
		halfEdgesIdx = halfEdgesIdx_first;
		for (unsigned int j = 0; j < halfEdgeScope_t1.size(); ++j) {
			TerrainHalfEdge *he_t1 = halfEdgeScope_t1.at(j);

			HalfEdge &he_arr = halfEdges.at(halfEdgesIdx++);

			// get halfEdges info indexes
			unsigned int he_nextIdx = he_t1->getNext()->getData().arrIndex;
			unsigned int he_originIdx = he_t1->getOrigin()->getData().arrIndex;
			unsigned int he_prevIdx = he_t1->getPrev()->getData().arrIndex;
			unsigned int he_twinIdx = he_t1->getTwin()->getData().arrIndex;

			if (he_nextIdx == -1) he_nextIdx = he_twinIdx;
			if (he_originIdx == -1) throw cpp::Exception("Origin should be inserted.");
			if (he_prevIdx == -1) he_prevIdx = he_twinIdx;
			if (he_twinIdx == -1) throw cpp::Exception("Twin should be inserted.");

			// copy halfEdges info
			he_arr.setFace(NULL); // originally null. To be handled.
			he_arr.setNext(&halfEdges[0] + he_nextIdx);
			he_arr.setOrigin(&vertices[0] + he_originIdx);
			he_arr.setPrev(&halfEdges[0] + he_prevIdx);
			he_arr.setTwin(&halfEdges[0] + he_twinIdx);

			// copy edgeData info
			if (he_arr.getData().edgeData == NULL) { // filtering twin edge that is already set
				// create an edgedata
				EdgeData::Source s;
				s.he = he_t1;
				s.v = v_t2;
				s.he_is_from_patch = false;
				EdgeData *ed_dest = createEdgeData();
				ed_dest->sources.push_back(s);
				he_arr.getData().edgeData = ed_dest;
				he_arr.getTwin()->getData().edgeData = ed_dest;
				++edgeDataIdx;
			}
			else { // when already set, we can get twin's origin.
				// considering halfEdge position (face is upward or downward), identify upedge and downedge.
				EdgeData *ed = he_arr.getData().edgeData;
				if (he_arr.getOrigin()->getData() < he_arr.getTwin()->getOrigin()->getData()) {
					ed->halfEdge_up = &he_arr;
					ed->halfEdge_down = he_arr.getTwin();
				}
				else {
					ed->halfEdge_up = he_arr.getTwin();
					ed->halfEdge_down = &he_arr;
				}
			}
		}
		// properly connect edges around the vertex outside the window
		for (unsigned int j = 0; j < vertexOutsideWindow_t1.size(); ++j) {
			TerrainVertex *v_t1 = vertexOutsideWindow_t1.at(j);
			Terrain::EdgeIterator eit(v_t1);
			TerrainHalfEdge *he_first_t1(NULL);
			HalfEdge *he_prev(NULL);
			while (eit.hasNext()) {
				TerrainHalfEdge *he_t1 = eit.getNext();
				if (he_t1->getData().arrIndex != -1) {
					HalfEdge *he = &halfEdges.at(he_t1->getData().arrIndex);
					if (he_first_t1 == NULL) {
						he_first_t1 = he_t1;
					}
					else {
						he_prev->getTwin()->setNext(he);
					}
					he_prev = he;
				}
			}
			he_prev->getTwin()->setNext(&halfEdges.at(he_first_t1->getData().arrIndex));
		}
		// restore markings
		for (unsigned int j = 0; j < halfEdgeScope_t1.size(); ++j) {
			TerrainHalfEdge *he_t1 = halfEdgeScope_t1.at(j);
			he_t1->getData().arrIndex = -1;
			he_t1->getOrigin()->getData().arrIndex = -1;
			he_t1->getTwin()->getOrigin()->getData().arrIndex = -1;
		}
	}

	// collect translations of (-t2) from the vertices of (t1). ((-t2 + v{t1}) ↘ grid ℅ nil ╰ v{t1} ↑ t2 + grid)
	std::vector<TerrainVertex *> transScope_t1;
	double transScope_t1_rangeX_min = t2->x_min + cur_x_min;
	double transScope_t1_rangeX_max = t2->x_max + cur_x_max;
	double transScope_t1_rangeY_min = t2->y_min + cur_y_min;
	double transScope_t1_rangeY_max = t2->y_max + cur_y_max;
	t1->appendVerticesInRange(transScope_t1_rangeX_min, transScope_t1_rangeX_max, transScope_t1_rangeY_min, transScope_t1_rangeY_max, &transScope_t1);

	// for each translation vectors, get copies inside the window.
	for (unsigned int i = 0; i < transScope_t1.size(); ++i) {
		TerrainVertex *v_t1 = transScope_t1.at(i);

		// collect structure of (-t2) in such range. (-t2 = grid range - v_t1).
		std::vector<TerrainVertex *> vertexScope_t2;
		std::vector<TerrainHalfEdge *> halfEdgeScope_t2;
		double vertexScope_t2_rangeX_min = -cur_x_max + v_t1->getData().p.x;
		double vertexScope_t2_rangeX_max = -cur_x_min + v_t1->getData().p.x;
		double vertexScope_t2_rangeY_min = -cur_y_max + v_t1->getData().p.y;
		double vertexScope_t2_rangeY_max = -cur_y_min + v_t1->getData().p.y;
		t2->appendVerticesInRange(vertexScope_t2_rangeX_min, vertexScope_t2_rangeX_max, vertexScope_t2_rangeY_min, vertexScope_t2_rangeY_max, &vertexScope_t2);

		// vertices in the window scope
		std::vector<TerrainVertex *> vertexOutsideWindow_t2;
		unsigned int halfEdgesIdx_first = halfEdgesIdx;
		unsigned int edgeDataIdx_first = edgeDataIdx;
		for (unsigned int j = 0; j < vertexScope_t2.size(); ++j) {
			TerrainVertex *v_t2 = vertexScope_t2.at(j);
			double v_arr_x = v_t1->getData().p.x - v_t2->getData().p.x;
			double v_arr_y = v_t1->getData().p.y - v_t2->getData().p.y;

			// create vertex info only if the vertex is in the window.
			if (isInWindow(v_arr_x, v_arr_y)) {
				Vertex *v_arr = createVertex();
				v_t2->getData().arrIndex = verticesIdx++; // mark vertex index
				v_arr->getData().x = v_arr_x;
				v_arr->getData().y = v_arr_y;
				v_arr->getData().insideWindow = true;

				// mark halfedge index around the vertex
				bool first = true;
				Terrain::EdgeIterator eit(v_t2);
				while (eit.hasNext()) {
					TerrainHalfEdge *he_t2 = eit.getNext();
					if (he_t2->getData().arrIndex == -1) {
						createHalfEdge(); // allocate space for the halfedge
						if (first) {
							v_arr->setIncidentEdge(&halfEdges.at(halfEdgesIdx));
							first = false;
						}
						he_t2->getData().arrIndex = halfEdgesIdx++; // mark halfEdge index
						halfEdgeScope_t2.push_back(he_t2); // collect halfedges in the window scope
					}
					else if (first) {
						v_arr->setIncidentEdge(&halfEdges.at(he_t2->getData().arrIndex));
						first = false;
					}
					// insert the twin halfedge.
					TerrainHalfEdge *he_twin_t2 = he_t2->getTwin();
					if (he_twin_t2->getData().arrIndex == -1) {
						HalfEdge *he_twin = createHalfEdge(); // allocate space for the halfedge
						he_twin_t2->getData().arrIndex = halfEdgesIdx++; // mark halfEdge index
						halfEdgeScope_t2.push_back(he_twin_t2); // collect halfedges in the window scope

						// insert a vertex just outside of the window.
						TerrainVertex *v_twin_t2 = he_twin_t2->getOrigin();
						if (v_twin_t2->getData().arrIndex == -1) {
							double v_twin_x = v_t1->getData().p.x - v_twin_t2->getData().p.x;
							double v_twin_y = v_t1->getData().p.y - v_twin_t2->getData().p.y;
							if (!isInWindow(v_twin_x, v_twin_y)) {
								Vertex *v_twin = createVertex();
								v_twin_t2->getData().arrIndex = verticesIdx++; // mark vertex index
								v_twin->getData().x = v_twin_x;
								v_twin->getData().y = v_twin_y;
								v_twin->getData().insideWindow = false;
								v_twin->setIncidentEdge(he_twin);
								vertexOutsideWindow_t2.push_back(v_twin_t2);
							}
						}
					}
				}
			}
		}

		// edges in the window scope
		halfEdgesIdx = halfEdgesIdx_first;
		for (unsigned int j = 0; j < halfEdgeScope_t2.size(); ++j) {
			TerrainHalfEdge *he_t2 = halfEdgeScope_t2.at(j);

			HalfEdge &he_arr = halfEdges.at(halfEdgesIdx++);

			// get halfEdges info indexes
			unsigned int he_nextIdx = he_t2->getNext()->getData().arrIndex;
			unsigned int he_originIdx = he_t2->getOrigin()->getData().arrIndex;
			unsigned int he_prevIdx = he_t2->getPrev()->getData().arrIndex;
			unsigned int he_twinIdx = he_t2->getTwin()->getData().arrIndex;

			if (he_nextIdx == -1) he_nextIdx = he_twinIdx;
			if (he_originIdx == -1) throw cpp::Exception("Origin should be inserted.");
			if (he_prevIdx == -1) he_prevIdx = he_twinIdx;
			if (he_twinIdx == -1) throw cpp::Exception("Twin should be inserted.");

			// copy halfEdges info
			he_arr.setFace(NULL); // originally null. To be handled.
			he_arr.setNext(&halfEdges[0] + he_nextIdx);
			he_arr.setOrigin(&vertices[0] + he_originIdx);
			he_arr.setPrev(&halfEdges[0] + he_prevIdx);
			he_arr.setTwin(&halfEdges[0] + he_twinIdx);

			// copy edgeData info
			if (he_arr.getData().edgeData == NULL) { // filtering twin edge that is already set
				// create an edgedata
				EdgeData::Source s;
				s.he = he_t2;
				s.v = v_t1;
				s.he_is_from_patch = false;
				EdgeData *ed_dest = createEdgeData();
				ed_dest->sources.push_back(s);
				he_arr.getData().edgeData = ed_dest;
				he_arr.getTwin()->getData().edgeData = ed_dest;
				++edgeDataIdx;
			}
			else { // when already set, we can get twin's origin.
				// considering halfEdge position (face is upward or downward), identify upedge and downedge.
				EdgeData *ed = he_arr.getData().edgeData;
				if (he_arr.getOrigin()->getData() < he_arr.getTwin()->getOrigin()->getData()) {
					ed->halfEdge_up = &he_arr;
					ed->halfEdge_down = he_arr.getTwin();
				}
				else {
					ed->halfEdge_up = he_arr.getTwin();
					ed->halfEdge_down = &he_arr;
				}
			}
		}

		// properly connect edges around the vertex outside the window
		for (unsigned int j = 0; j < vertexOutsideWindow_t2.size(); ++j) {
			TerrainVertex *v_t2 = vertexOutsideWindow_t2.at(j);
			Terrain::EdgeIterator eit(v_t2);
			TerrainHalfEdge *he_first_t2(NULL);
			HalfEdge *he_prev(NULL);
			while (eit.hasNext()) {
				TerrainHalfEdge *he_t2 = eit.getNext();
				if (he_t2->getData().arrIndex != -1) {
					HalfEdge *he = &halfEdges.at(he_t2->getData().arrIndex);
					if (he_first_t2 == NULL) {
						he_first_t2 = he_t2;
					}
					else {
						he_prev->getTwin()->setNext(he);
					}
					he_prev = he;
				}
			}
			he_prev->getTwin()->setNext(&halfEdges.at(he_first_t2->getData().arrIndex));
		}

		// restore markings
		for (unsigned int j = 0; j < halfEdgeScope_t2.size(); ++j) {
			TerrainHalfEdge *he_t2 = halfEdgeScope_t2.at(j);
			he_t2->getData().arrIndex = -1;
			he_t2->getOrigin()->getData().arrIndex = -1;
			he_t2->getTwin()->getOrigin()->getData().arrIndex = -1;
		}
	}

	/* insert window boundary */
	/*
	    v4   he3u   v3
	    ≒式式式式式≒
	    弛   he3d   弛
	he4u弛he4d  he2u弛he2d
	    弛   he1u   弛
	    ≒式式式式式≒
	    v1   he1d   v2
	*/

	Vertex *v1WB = createVertex();
	Vertex *v2WB = createVertex();
	Vertex *v3WB = createVertex();
	Vertex *v4WB = createVertex();
	HalfEdge *he1uWB = createHalfEdge();
	HalfEdge *he2uWB = createHalfEdge();
	HalfEdge *he3uWB = createHalfEdge();
	HalfEdge *he4uWB = createHalfEdge();
	HalfEdge *he1dWB = createHalfEdge();
	HalfEdge *he2dWB = createHalfEdge();
	HalfEdge *he3dWB = createHalfEdge();
	HalfEdge *he4dWB = createHalfEdge();
	EdgeData *ed1WB = createEdgeData();
	EdgeData *ed2WB = createEdgeData();
	EdgeData *ed3WB = createEdgeData();
	EdgeData *ed4WB = createEdgeData();
	v1WB->setIncidentEdge(he1uWB);
	v1WB->getData().x = cur_x_min;
	v1WB->getData().y = cur_y_min;
	v2WB->setIncidentEdge(he2uWB);
	v2WB->getData().x = cur_x_max;
	v2WB->getData().y = cur_y_min;
	v3WB->setIncidentEdge(he3dWB);
	v3WB->getData().x = cur_x_max;
	v3WB->getData().y = cur_y_max;
	v4WB->setIncidentEdge(he4dWB);
	v4WB->getData().x = cur_x_min;
	v4WB->getData().y = cur_y_max;
	he1uWB->setNext(he2uWB);
	he1uWB->setOrigin(v1WB);
	he1uWB->setPrev(he4dWB);
	he1uWB->setTwin(he1dWB);
	he1uWB->getData().edgeData = ed1WB;
	he1dWB->setNext(he4uWB);
	he1dWB->setOrigin(v2WB);
	he1dWB->setPrev(he2dWB);
	//he1dWB->setTwin(he1uWB);
	he1dWB->getData().edgeData = ed1WB;
	he2uWB->setNext(he3dWB);
	he2uWB->setOrigin(v2WB);
	//he2uWB->setPrev(he1uWB);
	he2uWB->setTwin(he2dWB);
	he2uWB->getData().edgeData = ed2WB;
	//he2dWB->setNext(he1dWB);
	he2dWB->setOrigin(v3WB);
	he2dWB->setPrev(he3uWB);
	//he2dWB->setTwin(he2uWB);
	he2dWB->getData().edgeData = ed2WB;
	//he3uWB->setNext(he2dWB);
	he3uWB->setOrigin(v4WB);
	he3uWB->setPrev(he4uWB);
	he3uWB->setTwin(he3dWB);
	he3uWB->getData().edgeData = ed3WB;
	he3dWB->setNext(he4dWB);
	he3dWB->setOrigin(v3WB);
	//he3dWB->setPrev(he2uWB);
	//he3dWB->setTwin(he3uWB);
	he3dWB->getData().edgeData = ed3WB;
	//he4uWB->setNext(he1dWB);
	he4uWB->setOrigin(v1WB);
	//he4uWB->setPrev(he1dWB);
	he4uWB->setTwin(he4dWB);
	he4uWB->getData().edgeData = ed4WB;
	//he4dWB->setNext(he1uWB);
	he4dWB->setOrigin(v4WB);
	//he4dWB->setPrev(he3dWB);
	//he4dWB->setTwin(he4uWB);
	he4dWB->getData().edgeData = ed4WB;
	ed1WB->halfEdge_up = he1uWB;
	ed1WB->halfEdge_down = he1dWB;
	ed2WB->halfEdge_up = he2uWB;
	ed2WB->halfEdge_down = he2dWB;
	ed3WB->halfEdge_up = he3uWB;
	ed3WB->halfEdge_down = he3dWB;
	ed4WB->halfEdge_up = he4uWB;
	ed4WB->halfEdge_down = he4dWB;

	// Take the firstHalfEdge of the current grid and that of the next grid.
	if (cur_y_grid % 2 == 0) {
		if (cur_x_grid == 0) { // first cell of the row
			firstHalfEdge = he1dWB;
			firstHalfEdge_next = he2uWB;
		}
		else if (cur_x_grid == x_gridSize - 1) { // last cell of the row
			firstHalfEdge = he4uWB;
			firstHalfEdge_next = he3dWB;
		}
		else {
			firstHalfEdge = he4uWB;
			firstHalfEdge_next = he2uWB;
		}
	}
	else {
		if (cur_x_grid == x_gridSize - 1) { // first cell of the row
			firstHalfEdge = he1dWB;
			firstHalfEdge_next = he4dWB;
		}
		else if (cur_x_grid == 0) { // last cell of the row
			firstHalfEdge = he2dWB;
			firstHalfEdge_next = he3dWB;
		}
		else {
			firstHalfEdge = he2dWB;
			firstHalfEdge_next = he4dWB;
		}
	}
}

void Arrangement::makeArrangement()
{
	// Run sweepline algorithm.
	sweepLine.initialize(this);
	sweepLine.run();

	// Erase outside-window structure.
	auto erased_it = erasedVerticesIndices.begin();
	for (unsigned int i = 0; i < vertices.size(); ++i) {
		if (erased_it != erasedVerticesIndices.end() && i == *erased_it) { // discard the erased vertex.
			++erased_it;
			continue;
		}
		Vertex &v = vertices.at(i);
		if (!isInWindow(v.getData().x, v.getData().y)) { // if the vertex is outside the window, erase adjacent structure.
			EdgeIterator eit(&v);
			while (eit.hasNext()) {
				HalfEdge *he = eit.getNext();
				HalfEdge *he_twin = he->getTwin();
				VertexData &vd_twin = he_twin->getOrigin()->getData();

				if (((nearlyEqual(vd_twin.x, cur_x_min) || nearlyEqual(vd_twin.x, cur_x_max)) && (cur_y_min <= vd_twin.y && vd_twin.y <= cur_y_max))
					|| ((nearlyEqual(vd_twin.y, cur_y_min) || nearlyEqual(vd_twin.y, cur_y_max)) && (cur_x_min <= vd_twin.x && vd_twin.x <= cur_x_max))) {
					// he_twin->origin is a boundary point, link the boundary structure.
#ifdef DEBUG
					if (he->getNext()->getOrigin() != he_twin->getOrigin())
						std::cerr << "ERROR vertex : " << i << '\n';
#endif
					he->getNext()->setPrev(he_twin->getPrev());
				}
				deleteEdgeData(he->getData().edgeData - &edgeDataContainer[0]);
				deleteHalfEdge(he->getTwin() - &halfEdges[0]);
				deleteHalfEdge(he - &halfEdges[0]);
			}
			deleteVertex(i);
		}
	}

	// Construct faces.
	faces.reserve(number_of_edges() - number_of_vertices() + 2); // Euler's law. e - v + d(=2) = f
	auto erased_heit = erasedHalfEdgesIndices.begin();
	for (unsigned int i = 0; i < halfEdges.size(); ++i) {
		if (erased_heit != erasedHalfEdgesIndices.end() && *erased_heit == i) { // discard the erased halfEdge.
			++erased_heit;
		}
		else {
			HalfEdge &he = halfEdges.at(i);
			if (he.getFace() == NULL) { // if face is not set yet, set face of all boundary halfEdges of the face.
				Face *f = createFace();
				f->setBoundary(&he);
				EdgeIterator eit(f);
				while (eit.hasNext()) {
					HalfEdge *eit_he = eit.getNext();
#ifdef DEBUG
					if (eit_he->getFace() != NULL)
						std::cout << "ERROR FACE : " << eit_he->getFace() - &faces[0] << '\n';
#endif
					eit_he->setFace(f);
					if (firstHalfEdge == eit_he) { // if he is firstHalfEdge, set this face to the outerface.
						setOuterface(f);
					}
				}
			}
		}
	}
#ifdef DEBUG
	if (faces.size() != number_of_edges() - number_of_vertices() + 2) {
		std::cerr << "# of edges = " << number_of_edges() << '\n';
		std::cerr << "# of vertices = " << number_of_vertices() << '\n';
		std::cerr << "# of faces = " << number_of_faces() << ' ' << faces.size() << '\n';
		throw cpp::Exception("Face num not match.");
	}
#endif
}

bool ArrangementEdgeDataCompare::operator()(const ArrangementEdgeData *ed1, const ArrangementEdgeData *ed2) const // return ed1 > ed2
{
	Arrangement::VertexData &vd1L = ed1->halfEdge_up->getOrigin()->getData();
	Arrangement::VertexData &vd1R = ed1->halfEdge_down->getOrigin()->getData();
	Arrangement::VertexData &vd2L = ed2->halfEdge_up->getOrigin()->getData();
	Arrangement::VertexData &vd2R = ed2->halfEdge_down->getOrigin()->getData();

	double sx = SweepLine::getX();
	
	// Compare y-coordinate of the intersection with the sweepLine.
	// (vd1R.y - vd1L.y) * (sx - vd1L.x) / (vd1R.x - vd1L.x) + vd1L.y = y1 > y2 = (vd2R.y - vd2L.y) * (sx - vd2L.x) / (vd2R.x - vd2L.x) + vd2L.y;
	if (vd1L.x == vd1R.x) { // if ed1 is vertical, consider as range.
		if (vd2L.x == vd2R.x) { // if ed2 is also vertical, consider as range.
			return vd1L.y > vd2R.y; // vd1L.y < vd1R.y and vd2L.y < vd2R.y, so y1 strictly larger than y2 is equivalent to ...
		}
		else {
			double divisor = vd2R.x - vd2L.x;
			return vd1L.y * divisor > (vd2R.y - vd2L.y) * (sx - vd2L.x) + vd2L.y * divisor; // vd1L.y < vd1R.y, so y1 strictly larger than y2 is equivalent to ...
		}
	}
	else { 
		if (vd2L.x == vd2R.x) {
			double divisor = vd1R.x - vd1L.x;
			return (vd1R.y - vd1L.y) * (sx - vd1L.x) + vd1L.y * divisor > vd2R.y * divisor; // vd2L.y < vd2R.y, so y1 strictly larger than y2 is equivalent to ...
		}
		else {
			double divisor1 = vd1R.x - vd1L.x;
			double divisor2 = vd2R.x - vd2L.x;
			double divisor12 = divisor1 * divisor2;
			return ((vd1R.y - vd1L.y) * (sx - vd1L.x) * divisor2 + vd1L.y * divisor12)
				 > ((vd2R.y - vd2L.y) * (sx - vd2L.x) * divisor1 + vd2L.y * divisor12);
		}
	}
}

// ed2 was below ed1. ed2 will go up, and ed1 will go down relatively.
Arrangement::Vertex* SweepLine::handleIntersectionEvent(ArrangementEdgeData *ed1, ArrangementEdgeData *ed2)
{
	//// check source
	//for (unsigned int i = 0; i < ed1->sources.size(); ++i) {
	//	for (unsigned int j = 0; j < ed2->sources.size(); ++j) {
	//		if (ed1->sources[i] == ed2->sources[j]) { // if the edges are from the same source, don't handle intersection event (already proved not to be intersect).
	//			return false;
	//		}
	//	}
	//}

	ArrangementVertexData &v1L(ed1->halfEdge_up->getOrigin()->getData());
	ArrangementVertexData &v1R(ed1->halfEdge_down->getOrigin()->getData());
	ArrangementVertexData &v2L(ed2->halfEdge_up->getOrigin()->getData());
	ArrangementVertexData &v2R(ed2->halfEdge_down->getOrigin()->getData());

	// filter boundary-intersecting set
	if (v1L == v2L || v1R == v2R) {
		return NULL;
	}

	// Ax + By = C
	double A_1 = v1R.y - v1L.y;
	double B_1 = v1L.x - v1R.x;
	double C_1 = A_1 * v1L.x + B_1 * v1L.y;
	double A_2 = v2R.y - v2L.y;
	double B_2 = v2L.x - v2R.x;
	double C_2 = A_2 * v2L.x + B_2 * v2L.y;

	// determinant
	double det = A_1 * B_2 - A_2 * B_1;
	if (det == 0) { // Lines are parallel
		return NULL;
	}

	// Intersection point of two "lines"
	double x_det = (B_2 * C_1 - B_1 * C_2);
	double y_det = (A_1 * C_2 - A_2 * C_1);
	std::pair<double, double> y1_range = std::minmax(v1L.y, v1R.y);
	std::pair<double, double> y2_range = std::minmax(v2L.y, v2R.y);
	std::pair<double, double> x_range(std::max(v1L.x, v2L.x), std::min(v1R.x, v2R.x));
	std::pair<double, double> y_range(std::max(y1_range.first, y2_range.first), std::min(y1_range.second, y2_range.second));

	// (x,y) in the segment range?
	if ((det > 0
		&& (nearlyInClosedRange(x_det, x_range.first * det, x_range.second * det) &&
			nearlyInClosedRange(y_det, y_range.first * det, y_range.second * det)))
	 || (det < 0
		&& (nearlyInClosedRange(x_det, x_range.second * det, x_range.first * det) &&
			nearlyInClosedRange(y_det, y_range.second * det, y_range.first * det))))
	{ // if two edges intersect, create new intersection event.
		ArrangementVertex *v = updateDCELIntersection(ed1, ed2, x_det / det, y_det / det);
#ifdef DEBUG
		std::cerr << "-- Intersection :" << v << " (" << x_det / det << ',' << y_det / det << ") " << ed1 << ' ' << ed2 << '\n';
		debug_handledVertices.insert(v);
#endif
		v->getData().it_eventQueue = events_insert(v);
		return v;
	}
	return NULL;
}

// ed2 was below ed1. ed2 will go up, and ed1 will go down relatively.
// returns the intersection vertex
Arrangement::Vertex *
SweepLine::updateDCELIntersection(ArrangementEdgeData *ed1, ArrangementEdgeData *ed2, double int_x, double int_y)
{
	ArrangementVertexData &v1L(ed1->halfEdge_up->getOrigin()->getData());
	ArrangementVertexData &v1R(ed1->halfEdge_down->getOrigin()->getData());
	ArrangementVertexData &v2L(ed2->halfEdge_up->getOrigin()->getData());
	ArrangementVertexData &v2R(ed2->halfEdge_down->getOrigin()->getData());

	ArrangementVertexData v_int;
	v_int.x = int_x;
	v_int.y = int_y;

	if (v1L == v_int) {
		if (v2L == v_int)
			throw cpp::Exception("twin-vertex should be handled before.");
		else if (v2R == v_int)
			throw cpp::Exception("twin-vertex should be handled before.");
		else {
			ArrangementVertex *v = ed1->halfEdge_up->getOrigin();
			updateDCELVertexEdgeIntersection(v, ed2, true);
			return v;
		}
	}
	else if (v1R == v_int) {
		if (v2L == v_int)
			throw cpp::Exception("twin-vertex should be handled before.");
		else if (v2R == v_int)
			throw cpp::Exception("twin-vertex should be handled before.");
		else {
			ArrangementVertex *v = ed1->halfEdge_down->getOrigin();
			updateDCELVertexEdgeIntersection(v, ed2, true);
			return v;
		}
	}
	else if (v2L == v_int) {
		if (v1L == v_int)
			throw cpp::Exception("twin-vertex should be handled before.");
		else if (v1R == v_int)
			throw cpp::Exception("twin-vertex should be handled before.");
		else {
			ArrangementVertex *v = ed2->halfEdge_up->getOrigin();
			updateDCELVertexEdgeIntersection(v, ed1, true);
			return v;
		}
	}
	else if (v2R == v_int) {
		if (v1L == v_int)
			throw cpp::Exception("twin-vertex should be handled before.");
		else if (v1R == v_int)
			throw cpp::Exception("twin-vertex should be handled before.");
		else {
			ArrangementVertex *v = ed2->halfEdge_down->getOrigin();
			updateDCELVertexEdgeIntersection(v, ed1, true);
			return v;
		}
	}
	else
		return updateDCELProperIntersection(ed1, ed2, int_x, int_y);
}

// ed2 was below ed1. ed2 will go up, and ed1 will go down relatively.
// returns the intersection vertex
Arrangement::Vertex * 
SweepLine::updateDCELProperIntersection(ArrangementEdgeData *ed1, ArrangementEdgeData *ed2, double int_x, double int_y)
{
	// existing sturcture of ed1
	ArrangementHalfEdge *ed1u = ed1->halfEdge_up;
	ArrangementHalfEdge *ed1d = ed1->halfEdge_down;

	// create new structure extended from ed1
	ArrangementEdgeData *ed1N = parent->createEdgeData();
	ArrangementHalfEdge *ed1Nu = parent->createHalfEdge();
	ArrangementHalfEdge *ed1Nd = parent->createHalfEdge();
	ArrangementVertex *v_int = parent->createVertex();
	v_int->getData().x = int_x;
	v_int->getData().y = int_y;

#ifdef DEBUG
	if (!(ed1u->getOrigin()->getData() < v_int->getData() && v_int->getData() < ed1d->getOrigin()->getData()))
		throw cpp::Exception("ed1 -> v_int -> ed1N");
#endif

	// link edges and halfedges
	ed1->halfEdge_down = ed1Nd;
	ed1N->halfEdge_up = ed1Nu;
	ed1N->halfEdge_down = ed1d;
	ed1N->sources = ed1->sources;
	ed1u->getData().edgeData = ed1;
	ed1Nd->getData().edgeData = ed1;
	ed1Nu->getData().edgeData = ed1N;
	ed1d->getData().edgeData = ed1N;

	// link halfedges and v_int
	ed1Nu->setNext(ed1u->getNext());
	ed1Nu->setOrigin(v_int);
	ed1Nu->setPrev(ed1u);
	ed1Nu->setTwin(ed1d);
	ed1Nd->setNext(ed1d->getNext());
	ed1Nd->setOrigin(v_int);
	ed1Nd->setPrev(ed1d);
	ed1Nd->setTwin(ed1u);
	v_int->setIncidentEdge(ed1Nu);

	// merge v_int to the edge ed2
	updateDCELVertexEdgeIntersection(v_int, ed2, true);

	// returns v_int
	return v_int;
}

// returns edN
ArrangementEdgeData *
SweepLine::updateDCELVertexEdgeIntersection(ArrangementVertex *v, ArrangementEdgeData *ed, bool survive_v)
{
	/*
			  v_dummy			   v弛				 edu  v弛  edNu
		(ed)-----≒-----(edN)	+	≒	  =		(ed)-------≒-------(edN)
									弛				 edNd  弛  edd
	*/
	
	// existing structure of ed
	ArrangementHalfEdge *edu = ed->halfEdge_up;
	ArrangementHalfEdge *edd = ed->halfEdge_down;

	// create new structure
	ArrangementEdgeData *edN = parent->createEdgeData();
	ArrangementHalfEdge *edNu = parent->createHalfEdge();
	ArrangementHalfEdge *edNd = parent->createHalfEdge();
	ArrangementVertex *v_dummy = parent->createVertex();
	v_dummy->getData().x = v->getData().x;
	v_dummy->getData().y = v->getData().y;

#ifdef DEBUG
	/*if (!(edu->getOrigin()->getData() < v_dummy->getData() && v_dummy->getData() < edd->getOrigin()->getData()))
		throw cpp::Exception("ed -> v -> edN");*/

	for (unsigned int i = 0; i < parent->halfEdges.size(); ++i) {
		if (parent->erasedHalfEdgesIndices.find(i) != parent->erasedHalfEdgesIndices.end())
			continue;
		if (parent->halfEdges[i].getOrigin() == v_dummy)
			std::cerr << "ERROR : " << i << '\n';
	}
#endif
	// link edges and halfedges
	ed->halfEdge_down = edNd;
	edN->halfEdge_up = edNu;
	edN->halfEdge_down = edd;
	edN->sources = ed->sources;
	edu->getData().edgeData = ed;
	edNd->getData().edgeData = ed;
	edNu->getData().edgeData = edN;
	edd->getData().edgeData = edN;

	// link halfedges and v_dummy
	edNu->setNext(edu->getNext());
	edNu->setOrigin(v_dummy);
	edNu->setPrev(edu);
	edNu->setTwin(edd);
	edNd->setNext(edd->getNext());
	edNd->setOrigin(v_dummy);
	edNd->setPrev(edd);
	edNd->setTwin(edu);
	v_dummy->setIncidentEdge(edNu);

	// merge v_dummy to v
	if (survive_v)
		updateDCEL2VertexIntersection(v, v_dummy);
	else {
		updateDCEL2VertexIntersection(v_dummy, v); // event of v will be deleted
		events_insert(v_dummy); // newly insert event
	}

	return edN;
}

void SweepLine::updateDCEL2VertexIntersection(ArrangementVertex *v, ArrangementVertex *v_del)
{
	if (v == v_del) return; // if two are the same point, discard the operation.

	ArrangementVector vec_v(v->getData());
	ArrangementVector vec_v_del(v_del->getData());
	ArrangementVector vec_zero(0, 0);

	// check the same point
	if (!(v->getData() == v_del->getData())) 
		throw cpp::Exception("v and v_del are not the same point.");
	
	// merge incident edges of v_del into v
	/*
		忙------- he_prev (from v)
		弛
		戍------- he (from v_del)
		弛
		戌------- he_next (from v)
	*/

	// initialize edges from v
	Arrangement::EdgeIterator eit(v);
	ArrangementHalfEdge *he_prev = eit.getNext();
	ArrangementHalfEdge *he_next = eit.getNext();

	// copy incident edges of v_del
	std::vector<ArrangementHalfEdge *> incidentHalfEdges_v_del;
	Arrangement::EdgeIterator eit_del(v_del);
	while (eit_del.hasNext()) {
		ArrangementHalfEdge *he = eit_del.getNext();
		incidentHalfEdges_v_del.push_back(he);
	}

	// if v has only one edge...
	while (he_next == NULL) { 
		ArrangementHalfEdge *he = incidentHalfEdges_v_del.back();
		incidentHalfEdges_v_del.pop_back();
		he->setPrev(he_prev->getTwin());
		he->getTwin()->setNext(he_prev);
		he->setOrigin(v);
		he_next = he;

		// twin-edge handling (prev and he)
		ArrangementVector vec_he_prev = ArrangementVector(he_prev->getTwin()->getOrigin()->getData()) - vec_v;
		ArrangementVector vec_he = ArrangementVector(he->getTwin()->getOrigin()->getData()) - vec_v;
		if (vec_he_prev.det(vec_he) == 0) { // if prev and he have the same bi-direction,
			if (vec_zero > vec_he && vec_zero > vec_he_prev) { // he and he_prev are decreasing edges.
				throw cpp::Exception("Decreasing edges...1");
			}
			// cut and merge as a twin-edge.
			if ((vec_he > vec_zero && vec_he_prev > vec_he) || (vec_zero > vec_he && vec_he > vec_he_prev)) { // he is shorter than he_prev
				// link he to v. he will be the only edge of v.
				he->setOrigin(v);
				he->setPrev(he->getTwin());
				// cut he_prev. (actually shrink he_prev to he->he_prev)
				v->setIncidentEdge(he);
				ArrangementVertex *v_new = parent->createVertex();
				ArrangementVertex *v_hetwin = he->getTwin()->getOrigin();
				v_new->getData().insideWindow = v_hetwin->getData().insideWindow;
				v_new->getData().x = v_hetwin->getData().x;
				v_new->getData().y = v_hetwin->getData().y;
				v_new->setIncidentEdge(he_prev);
				he_prev->setOrigin(v_new);
				he_prev->setPrev(he_prev->getTwin());
				events_insert(v_new);
				// BBT handling
				if (he_prev->getData().edgeData->inEdgeDataBBT) {
					if (he->getData().edgeData->inEdgeDataBBT)
						throw cpp::Exception("two twin-edges are in BBT.");
					auto it = edgeDataBBT_erase(he_prev->getData().edgeData);
					edgeDataBBT_insert(it, he->getData().edgeData);
				}
				// copy source data.
				he->getData().edgeData->sources.insert(he->getData().edgeData->sources.end(), he_prev->getData().edgeData->sources.begin(), he_prev->getData().edgeData->sources.end());
			}
			else if ((vec_he_prev > vec_zero && vec_he > vec_he_prev) || (vec_zero > vec_he_prev && vec_he_prev > vec_he)) { // he_prev is shorter than he
				// cut he. (actually shrink he_prev to he_prev->he)
				ArrangementVertex *v_new = parent->createVertex();
				ArrangementVertex *v_heprevtwin = he_prev->getTwin()->getOrigin();
				v_new->getData().insideWindow = v_heprevtwin->getData().insideWindow;
				v_new->getData().x = v_heprevtwin->getData().x;
				v_new->getData().y = v_heprevtwin->getData().y;
				v_new->setIncidentEdge(he);
				he->setOrigin(v_new);
				he->setPrev(he->getTwin());
				events_insert(v_new);
				// BBT handling
				if (he->getData().edgeData->inEdgeDataBBT) {
					if (he_prev->getData().edgeData->inEdgeDataBBT)
						throw cpp::Exception("two twin-edges are in BBT.");
					auto it = edgeDataBBT_erase(he->getData().edgeData);
					edgeDataBBT_insert(it, he_prev->getData().edgeData);
				}
				// copy source data.
				he_prev->getData().edgeData->sources.insert(he_prev->getData().edgeData->sources.end(), he->getData().edgeData->sources.begin(), he->getData().edgeData->sources.end());
			}
			else if (vec_he_prev == vec_he) { // length of he_prev and he are the same
				// remove the edge of he. (lazy) (handle events also)
				std::pair<bool, bool> vDelRes = parent->deleteEdge(he->getData().edgeData);
				if (vDelRes.first) events_erase(he->getData().edgeData->halfEdge_up->getOrigin());
				if (vDelRes.second) events_erase(he->getData().edgeData->halfEdge_down->getOrigin());
				// BBT handling
				if (he->getData().edgeData->inEdgeDataBBT) {
					if (he_prev->getData().edgeData->inEdgeDataBBT)
						throw cpp::Exception("two twin-edges are in BBT.");
					auto it = edgeDataBBT_erase(he->getData().edgeData);
					edgeDataBBT_insert(it, he_prev->getData().edgeData);
				}
				// copy source data.
				he_prev->getData().edgeData->sources.insert(he_prev->getData().edgeData->sources.end(), he->getData().edgeData->sources.begin(), he->getData().edgeData->sources.end());
			}
			else
				throw cpp::Exception("Length should be shorter, longer, or the same.");

			if (incidentHalfEdges_v_del.empty()) // if there is no other edge, escape.
				break;
			else // if there is another edge, continue to find.
				he_next = NULL;	
		}
	}

	// traverse incident edges of v_del and attach them to v
	ArrangementVector vec_he_prev = ArrangementVector(he_prev->getTwin()->getOrigin()->getData()) - vec_v;
	ArrangementVector vec_he_next = ArrangementVector(he_next->getTwin()->getOrigin()->getData()) - vec_v;
	for (unsigned int i = 0; i < incidentHalfEdges_v_del.size(); ++i) {
		ArrangementHalfEdge *he = incidentHalfEdges_v_del[i];
		ArrangementVector vec_he = ArrangementVector(he->getTwin()->getOrigin()->getData()) - vec_v_del;
		if (vec_he == vec_zero)
			throw cpp::Exception("An edge should not be zero-length (vec_he).");

		// traverse incident edges of v while he becomes between he_prev and he_next.
		bool insert = false;
		while (!insert) {
			if (vec_he_prev.isLeftFrom(vec_he_next)) { // not reflex angle
				if (!vec_he.isLeftFrom(vec_he_prev) && vec_he.isLeftFrom(vec_he_next)) insert = true; // if he is in [prev,next)
			}
			else if (vec_he_next.isLeftFrom(vec_he_prev)) { // reflex angle
				if (!vec_he.isLeftFrom(vec_he_prev) || vec_he.isLeftFrom(vec_he_next)) insert = true; // if he is in [prev,next)
			}
			else if (vec_he_prev * vec_he_next < 0) { // 180-degree angle
				if (vec_he_prev.x > 0) {
					if (vec_he.x > 0 && !vec_he.isLeftFrom(vec_he_prev)) insert = true; // if he is in [prev,(0,-1))
					else if (vec_he.x <= 0 && vec_he.isLeftFrom(vec_he_next)) insert = true; // if he is in [(0,-1),next)
				}
				else if (vec_he_prev.x < 0) {
					if (vec_he.x < 0 && !vec_he.isLeftFrom(vec_he_prev)) insert = true; // if he is in [prev,(0,1))
					else if (vec_he.x >= 0 && vec_he.isLeftFrom(vec_he_next)) insert = true; // if he is in [(0,1),next)
				}
				else { // if vec_he_prev.x == 0
					if (vec_he_prev.y * vec_he.y > 0 && !vec_he.isLeftFrom(vec_he_prev)) insert = true; // if he is in [prev,(1,0) or (-1,0))
					else if (vec_he_prev.y * vec_he.y <= 0 && vec_he.isLeftFrom(vec_he_next)) insert = true; // if he is in [(1,0) or (-1,0)),next)
				}
			}
			else if (vec_he_prev * vec_he_next > 0) { // 0-degree angle
				// twin-edge handling (prev and next)
				if (vec_zero > vec_he_next && vec_zero > vec_he_prev) { // he and he_prev are decreasing edges.
					throw cpp::Exception("Decreasing edges...2");
				}
				if (vec_he_prev > vec_zero && vec_he_next > vec_zero) // if prev and next have the forward direction, merge them as a twin-edge.
					// cut and merge as a twin-edge.
					if ((vec_he_next > vec_zero && vec_he_prev > vec_he_next) || (vec_zero > vec_he_next && vec_he_next > vec_he_prev)) { // he_next is shorter than he_prev
						// cut he_prev. (actually shrink he_prev to he_next->he_prev)
						v->setIncidentEdge(he_next);
						ArrangementVertex *v_new = parent->createVertex();
						ArrangementVertex *v_henexttwin = he_next->getTwin()->getOrigin();
						v_new->getData().insideWindow = v_henexttwin->getData().insideWindow;
						v_new->getData().x = v_henexttwin->getData().x;
						v_new->getData().y = v_henexttwin->getData().y;
						v_new->setIncidentEdge(he_prev);
						he_prev->setOrigin(v_new);
						he_prev->setPrev(he_prev->getTwin());
						events_insert(v_new);
						// BBT handling
						if (he_prev->getData().edgeData->inEdgeDataBBT) {
							if (he_next->getData().edgeData->inEdgeDataBBT)
								throw cpp::Exception("two twin-edges are in BBT.");
							auto it = edgeDataBBT_erase(he_prev->getData().edgeData);
							edgeDataBBT_insert(it, he_next->getData().edgeData);
						}
						// copy source data.
						he_next->getData().edgeData->sources.insert(he_next->getData().edgeData->sources.end(), he_prev->getData().edgeData->sources.begin(), he_prev->getData().edgeData->sources.end());
					}
					else if ((vec_he_prev > vec_zero && vec_he_next > vec_he_prev) || (vec_zero > vec_he_prev && vec_he_prev > vec_he_next)) { // he_prev is shorter than he_next
						// cut he_next. (actually shrink he_prev to he_prev->he_next)
						ArrangementVertex *v_new = parent->createVertex();
						ArrangementVertex *v_heprevtwin = he_prev->getTwin()->getOrigin();
						v_new->getData().insideWindow = v_heprevtwin->getData().insideWindow;
						v_new->getData().x = v_heprevtwin->getData().x;
						v_new->getData().y = v_heprevtwin->getData().y;
						v_new->setIncidentEdge(he_next);
						he_next->setOrigin(v_new);
						he_next->setPrev(he_next->getTwin());
						events_insert(v_new);
						// BBT handling
						if (he_next->getData().edgeData->inEdgeDataBBT) {
							if (he_prev->getData().edgeData->inEdgeDataBBT)
								throw cpp::Exception("two twin-edges are in BBT.");
							auto it = edgeDataBBT_erase(he_next->getData().edgeData);
							edgeDataBBT_insert(it, he_prev->getData().edgeData);
						}
						// copy source data.
						he_prev->getData().edgeData->sources.insert(he_prev->getData().edgeData->sources.end(), he_next->getData().edgeData->sources.begin(), he_next->getData().edgeData->sources.end());
					}
					else if (vec_he_prev == vec_he_next) { // length of he_prev and he_next are the same
						// remove the edge of he_next. (lazy) (handle events also)
						std::pair<bool, bool> vDelRes = parent->deleteEdge(he_next->getData().edgeData);
						if (vDelRes.first) events_erase(he_next->getData().edgeData->halfEdge_up->getOrigin());
						if (vDelRes.second) events_erase(he_next->getData().edgeData->halfEdge_down->getOrigin());
						// BBT handling
						if (he_next->getData().edgeData->inEdgeDataBBT) {
							if (he_prev->getData().edgeData->inEdgeDataBBT)
								throw cpp::Exception("two twin-edges are in BBT.");
							auto it = edgeDataBBT_erase(he_next->getData().edgeData);
							edgeDataBBT_insert(it, he_prev->getData().edgeData);
						}
						// copy source data.
						he_prev->getData().edgeData->sources.insert(he_prev->getData().edgeData->sources.end(), he_next->getData().edgeData->sources.begin(), he_next->getData().edgeData->sources.end());
					}
					else
						throw cpp::Exception("Length should be shorter, longer, or the same.");
			}
			else {
				throw cpp::Exception("An edge should not be zero-length (vec_he_prev or vec_he_next).");
			}

			if (!insert) { // if not insert, get next range of incident edges of v
				if (!eit.hasNext()) 
					eit.reset(); // rewind iterator
				he_prev = he_next;
				vec_he_prev = vec_he_next;
				he_next = eit.getNext();
				vec_he_next = ArrangementVector(he_next->getTwin()->getOrigin()->getData()) - vec_v;
			}
			else { // if insert, 
				if (he_prev == he) {
					continue; // filter same-source edges
				}

				// twin-edge handling (prev and he)
				if (vec_he_prev.det(vec_he) == 0) { // if prev and he have the same bidirection.
					if (vec_he_prev * vec_he < 0) { // 180-degree.
						// link he to v between he_prev and he_next
						he->setPrev(he_prev->getTwin());
						he->getTwin()->setNext(he_next);
						he->setOrigin(v);
						// set insert edge as the "prev"
						he_prev = he;
						vec_he_prev = vec_he;
					}
					else if (vec_he_prev * vec_he > 0) { // same direction.
						if (vec_zero > vec_he && vec_zero > vec_he_prev) { // he and he_prev are decreasing edges.
							throw cpp::Exception("Decreasing edges...3");
						}
						// cut and merge as a twin-edge.
						if ((vec_he > vec_zero && vec_he_prev > vec_he) || (vec_zero > vec_he && vec_he > vec_he_prev)) { // he is shorter than he_prev
							// link he to v between he_prev and he_next
							he->setOrigin(v);
							he->setPrev(he_prev->getPrev());
							he->getTwin()->setNext(he_next);
							// cut he_prev. (actually shrink he_prev to he->he_prev)
							v->setIncidentEdge(he);
							ArrangementVertex *v_new = parent->createVertex();
							ArrangementVertex *v_hetwin = he->getTwin()->getOrigin();
							v_new->getData().insideWindow = v_hetwin->getData().insideWindow;
							v_new->getData().x = v_hetwin->getData().x;
							v_new->getData().y = v_hetwin->getData().y;
							v_new->setIncidentEdge(he_prev);
							he_prev->setOrigin(v_new);
							he_prev->setPrev(he_prev->getTwin());
							events_insert(v_new);
							// BBT handling
							if (he_prev->getData().edgeData->inEdgeDataBBT) {
								if (he->getData().edgeData->inEdgeDataBBT)
									throw cpp::Exception("two twin-edges are in BBT.");
								auto it = edgeDataBBT_erase(he_prev->getData().edgeData);
								edgeDataBBT_insert(it, he->getData().edgeData);
							}
							// copy source data.
							he->getData().edgeData->sources.insert(he->getData().edgeData->sources.end(), he_prev->getData().edgeData->sources.begin(), he_prev->getData().edgeData->sources.end());
						}
						else if ((vec_he_prev > vec_zero && vec_he > vec_he_prev) || (vec_zero > vec_he_prev && vec_he_prev > vec_he)) { // he_prev is shorter than he
							// cut he. (actually shrink he_prev to he_prev->he)
							ArrangementVertex *v_new = parent->createVertex();
							ArrangementVertex *v_heprevtwin = he_prev->getTwin()->getOrigin();
							v_new->getData().insideWindow = v_heprevtwin->getData().insideWindow;
							v_new->getData().x = v_heprevtwin->getData().x;
							v_new->getData().y = v_heprevtwin->getData().y;
							v_new->setIncidentEdge(he);
							he->setOrigin(v_new);
							he->setPrev(he->getTwin());
							events_insert(v_new);
							// BBT handling
							if (he->getData().edgeData->inEdgeDataBBT) {
								if (he_prev->getData().edgeData->inEdgeDataBBT)
									throw cpp::Exception("two twin-edges are in BBT.");
								auto it = edgeDataBBT_erase(he->getData().edgeData);
								edgeDataBBT_insert(it, he_prev->getData().edgeData);
							}
							// copy source data.
							he_prev->getData().edgeData->sources.insert(he_prev->getData().edgeData->sources.end(), he->getData().edgeData->sources.begin(), he->getData().edgeData->sources.end());
						}
						else if (vec_he_prev == vec_he) { // length of he_prev and he are the same
							// remove the edge of he. (lazy) (handle events also)
							std::pair<bool, bool> vDelRes = parent->deleteEdge(he->getData().edgeData);
							if (vDelRes.first) events_erase(he->getData().edgeData->halfEdge_up->getOrigin());
							if (vDelRes.second) events_erase(he->getData().edgeData->halfEdge_down->getOrigin());
							// BBT handling
							if (he->getData().edgeData->inEdgeDataBBT) {
								if (he_prev->getData().edgeData->inEdgeDataBBT)
									throw cpp::Exception("two twin-edges are in BBT.");
								auto it = edgeDataBBT_erase(he->getData().edgeData);
								edgeDataBBT_insert(it, he_prev->getData().edgeData);
							}
							// copy source data.
							he_prev->getData().edgeData->sources.insert(he_prev->getData().edgeData->sources.end(), he->getData().edgeData->sources.begin(), he->getData().edgeData->sources.end());
						}
						else
							throw cpp::Exception("Length should be shorter, longer, or the same.");
					}
				}
				else {
					// link he to v between he_prev and he_next
					he->setPrev(he_prev->getTwin());
					he->getTwin()->setNext(he_next);
					he->setOrigin(v);
					// set insert edge as the "prev"
					he_prev = he;
					vec_he_prev = vec_he;
				}
			}
		}
	}

	// erase v_del (lazy)
	unsigned int id = v_del - &(parent->vertices[0]);
	parent->deleteVertex(id);
	events_erase(v_del);

#ifdef DEBUG
	//std::cerr << "--- TwoVertexIntersection Result ---\n";
	//std::cerr << "vertex : " << v << " (" << v->getData().x << ',' << v->getData().y << ")\n";
	//Arrangement::EdgeIterator eit_debug(v);
	//while (eit_debug.hasNext()) {
	//	ArrangementHalfEdge *he = eit_debug.getNext();
	//	std::cerr << he->getData().edgeData << '(' << he << ',' << he->getTwin() << ')' << " - ";
	//	he->getData().edgeData->print(std::cerr);
	//	std::cerr << '\n';
	//}
	//std::cerr << '\n';

	for (unsigned int i = 0; i < parent->halfEdges.size(); ++i) {
		if (parent->erasedHalfEdgesIndices.find(i) != parent->erasedHalfEdgesIndices.end())
			continue;
		if (parent->halfEdges[i].getOrigin() == v_del)
			std::cerr << "ERROR : " << i << '\n';
	}
#endif
}

void SweepLine::initialize(Arrangement *_parent)
{
	parent = _parent;

	// Initialize event queue
	events = EventQueue();
	firstEvent = true;

	// Insert event points related to the vertices
	auto erased_it = parent->erasedVerticesIndices.begin();
	for (unsigned int i = 0; i < parent->vertices.size(); ++i)
	{
		if (erased_it != parent->erasedVerticesIndices.end() && *erased_it == i) {
			++erased_it;
			continue;
		}
		ArrangementVertex *v = &parent->vertices.at(i);
		v->getData().it_eventQueue = events_insert(v);
	}
}

void SweepLine::advance()
{
	// Take the first event
	Event ep = events_popfront();
	++eventCount;
	currentEvent = &ep;

#ifdef DEBUG
	std::cerr << "================== Event " << eventCount << " ======================\n";
	std::cerr << "ep = (" << ep.v->getData().x << ',' << ep.v->getData().y << ')' << ep.v << "\n";
#endif

	// while two event points have the same position, merge the point to ep.
	while (!events.empty() 
		&& events.begin()->v->getData().x == ep.v->getData().x 
		&& events.begin()->v->getData().y == ep.v->getData().y) {


		Event ep_next = events_popfront();
#ifdef DEBUG
		std::cerr << "ep_next = (" << ep_next.v->getData().x << ',' << ep_next.v->getData().y << ')' << ep_next.v << "\n";
#endif
		updateDCEL2VertexIntersection(ep.v, ep_next.v);
	}

	// merge edges that passes ep.
	ArrangementEdgeData sample_zero_ep; // an edge of ep.v -> ep.v
	sample_zero_ep.halfEdge_up = ep.v->getIncidentEdge();
	sample_zero_ep.halfEdge_down = ep.v->getIncidentEdge();
	auto eqRange = edgeDataBBT.equal_range(&sample_zero_ep);
	for (auto eqRange_it = eqRange.first; eqRange_it != eqRange.second; ++eqRange_it) {
		ArrangementEdgeData *ed = *eqRange_it;
		if (ed->halfEdge_up->getOrigin() != ep.v && ed->halfEdge_down->getOrigin() != ep.v) { // not a ep.v edge
			updateDCELVertexEdgeIntersection(ep.v, ed, true);
		}
	}

#ifdef DEBUG
	std::cerr << "Merging same-position events done. Arrangement structure.\n";
	Arrangement::EdgeIterator eit_event(ep.v);
	while (eit_event.hasNext()) {
		ArrangementEdgeData *ed = eit_event.getNext()->getData().edgeData;
		std::cerr << ed << " = ";
		ed->print(std::cerr);
		std::cerr << '\n';
	}
	std::cerr << "Merging same-position events done.\n";
	EdgeDataBBTIterator it_debug = edgeDataBBT.begin();
	while (it_debug != edgeDataBBT.end()) {
		std::cerr << *it_debug << " = ";
		(*it_debug)->print(std::cerr);
		std::cerr << '\n';
		++it_debug;
	}
	std::cerr << '\n';
#endif

	// Insert the edges before current sweepline.
	Arrangement::EdgeIterator eit(ep.v);

	// Find the first before-event edge and the first after-event edge.
	ArrangementVector vec_v(ep.v->getData());
	ArrangementHalfEdge *BEedge(NULL), *AEedge(NULL);
	ArrangementVector vec_edge_before_event, vec_edge_after_event;
	while (eit.hasNext()) {
		ArrangementHalfEdge *he = eit.getNext();
		const Arrangement::VertexData vd = he->getTwin()->getOrigin()->getData();
		ArrangementVector vec_edge = ArrangementVector(vd) - vec_v;

		if (vd < ep.v->getData())
		{ // if he is the first before-event edge, set it to the first edge of before-edges.
			if (BEedge == NULL || vec_edge.isLeftFrom(vec_edge_before_event)) {
				BEedge = he;
				vec_edge_before_event = vec_edge;
			}
		}
		else if (ep.v->getData() < vd) { // if he is the first after-event edge, set it to the first edge of after-edges.
			if (AEedge == NULL || vec_edge.isLeftFrom(vec_edge_after_event)) {
				AEedge = he;
				vec_edge_after_event = vec_edge;
			}
		}
		else {
			throw cpp::Exception("zero-edge cannot be exist.");
		}
	}

	// Erase mode. 〢
	EdgeDataBBTIterator bbt_entry; // remember the position that edges erased.
	bbt_entry._Ptr = NULL;
	if (BEedge != NULL) {

		ep.v->setIncidentEdge(BEedge);
		Arrangement::EdgeIterator BEeit(ep.v);
		ArrangementHalfEdge *he = BEeit.getNext();
		bbt_entry = he->getData().edgeData->it_edgeDataBBT;
		//bbt_entry = edgeDataBBT.upper_bound(BEedge->getData().edgeData);
		do { // erase all the existing edges containing ep.v.
			bbt_entry = edgeDataBBT_erase(he->getData().edgeData);
			he = BEeit.getNext();
		} while (he != AEedge && he != BEedge);

		if (AEedge == NULL && bbt_entry != edgeDataBBT.end() && bbt_entry != edgeDataBBT.begin()) {
			// if erased edge is not the highest or the lowest and there is no edges to be continuously inserted, check intersection event.
			ArrangementEdgeData *ed1 = *bbt_entry;
			ArrangementEdgeData *ed2 = *--bbt_entry;
			handleIntersectionEvent(ed2, ed1);
		}
	}

#ifdef DEBUG
	std::cerr << "Erase mode done.\n";
	it_debug = edgeDataBBT.begin();
	while (it_debug != edgeDataBBT.end()) {
		std::cerr << *it_debug << " = ";
		(*it_debug)->print(std::cerr);
		std::cerr << '\n';
		++it_debug;
	}
	std::cerr << '\n';
#endif

	// Insert mode. 〣
	if (AEedge != NULL) {

		ep.v->setIncidentEdge(AEedge);
		Arrangement::EdgeIterator AEeit(ep.v);
		ArrangementHalfEdge *he = AEeit.getNext();
		//EdgeDataBBTIterator eit_same_y = edgeDataBBT.find(he->getData().edgeData);
		//if (eit_same_y != edgeDataBBT.end()) { // if there is already an edge passing ep.v (so cannot be inserted), error.
		//	throw cpp::Exception("An edge passing ep.v should be handled before.");
		//}
		EdgeDataBBTIterator inserted_lowerbound, inserted_upperbound;
		bool after_first = false;
		EdgeDataBBTIterator hintit;
		if (bbt_entry._Ptr == NULL)
			hintit = edgeDataBBT.upper_bound(he->getData().edgeData);
		else
			hintit = bbt_entry;
		do { // insert only when he does not reach to BEedge or AEedge yet.
			EdgeDataBBTIterator hintit_prev = hintit;
			EdgeDataBBTIterator inserted_it = edgeDataBBT_insert(hintit, he->getData().edgeData);
			// update lowerbound and upperbound.
			if (!after_first) {
				after_first = true;
				inserted_lowerbound = inserted_it;
			}
			inserted_upperbound = inserted_it;

			hintit = ++inserted_it; // update hint.
			he = AEeit.getNext(); // next step.
		} while (he != BEedge && he != AEedge);

		// handle intersection with lower neighbor and higher neighbor
		if (inserted_lowerbound != edgeDataBBT.begin()) { // if it is not the lowest, check intersection with the lower neighbor.
			EdgeDataBBTIterator it_above = inserted_lowerbound;
			if (ep.v == handleIntersectionEvent(*--it_above, *inserted_lowerbound)) { // if it_above passes ep.v,
				throw cpp::Exception("Cannot be done. (it_above passes ep.v)");
			}
		}
		if (inserted_upperbound != --edgeDataBBT.end()) { // if it is not the highest, check intersection with the higher neighbor.
			EdgeDataBBTIterator it_below = inserted_upperbound;
			if (ep.v == handleIntersectionEvent(*inserted_upperbound, *++it_below)) { // if it_below passes ep.v,
				throw cpp::Exception("Cannot be done. (it_below passes ep.v)");
			}
		}
	}
	
#ifdef DEBUG
	std::cerr << "Insert mode done.\n";
	it_debug = edgeDataBBT.begin();
	while (it_debug != edgeDataBBT.end()) {
		std::cerr << *it_debug <<  " = ";
		(*it_debug)->print(std::cerr);
		std::cerr << '\n';
		++it_debug;
	}
	std::cerr << '\n';
#endif
}

