#include "arrangement.h"

#include <algorithm>

SweepLine Arrangement::sweepLine;

SweepLine::EventQueue SweepLine::events;
SweepLine::EdgeDataBBT SweepLine::edgeDataBBT;
int SweepLine::eventCount = 0;

template <class T>
class Vector2
{
private:
	double x, y;

public:
	Vector2() {}
	Vector2(double _x, double _y) : x(_x), y(_y) {}
	Vector2(ArrangementVertexData &vd) : x(vd.x), y(vd.y) {}

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
		return det(v) > 0;
	}
};

bool arrVertexCompare(Arrangement::Vertex &v1, Arrangement::Vertex &v2)
{
	ArrangementVertexData &vd1 = v1.getData();
	ArrangementVertexData &vd2 = v2.getData();
	return (vd1.x < vd2.x) || (vd1.x == vd2.x && vd1.y < vd2.y);
}

///

Arrangement::Arrangement(Terrain *t1, Terrain *t2)
{
	// Datastructure of t1 & t2
    TerrainMesh &mesh_t1 = t1->getMesh();
    TerrainMesh &mesh_t2 = t2->getMesh();
	std::vector<TerrainVertex> &vertices_t1 = mesh_t1.getVertices();
	std::vector<TerrainVertex> &vertices_t2 = mesh_t2.getVertices();
	std::vector<TerrainHalfEdge> &edges_t1 = mesh_t1.getHalfEdges();
	std::vector<TerrainHalfEdge> &edges_t2 = mesh_t2.getHalfEdges();
	std::vector<TerrainFace> &faces_t1 = mesh_t1.getFaces();
	std::vector<TerrainFace> &faces_t2 = mesh_t2.getFaces();

	// Allocate space for vectors
	vertices.resize(2 * vertices_t1.size() * vertices_t2.size());
	edges.resize(vertices_t2.size() * edges_t1.size() + vertices_t1.size() * edges_t2.size());
	faces.resize(vertices_t2.size() * faces_t1.size() + vertices_t1.size() * faces_t2.size());
	unsigned int vertices_first_idx = 0;
	unsigned int edges_first_idx = 0;
	unsigned int faces_first_idx = 0;

	// nm^2 estimation of the overlay size.
	vertices.reserve(10 * vertices_t1.size() * vertices_t2.size() * vertices_t2.size());
	edges.reserve(80 * vertices_t1.size() * vertices_t2.size() * vertices_t2.size());
	faces.reserve(20 * vertices_t1.size() * vertices_t2.size() * vertices_t2.size());
	edgeDataContainer.reserve(40 * vertices_t1.size() * vertices_t2.size() * vertices_t2.size());

	// Insert all copied (t1) with translation of -(vertex of t2).
	for (unsigned int i=0; i < vertices_t2.size(); ++i)
	{
		TerrainVertex &v_t2 = vertices_t2.at(i);

		// copy vertex info
		for (unsigned int j=0; j < vertices_t1.size(); ++j)
		{
			Vertex &v = vertices.at(vertices_first_idx + j);
			TerrainVertex &v_t1 = vertices_t1.at(j);
			v.getData().x = v_t1.getData().p.x - v_t2.getData().p.x;
			v.getData().y = v_t1.getData().p.y - v_t2.getData().p.y;

			if (v_t1.getIncidentEdge() == NULL) v.setIncidentEdge(NULL);
			else v.setIncidentEdge(&edges[0] + edges_first_idx
				+ mesh_t1.getHalfEdgeId(v_t1.getIncidentEdge()));
		}
		// copy edge info
		for (unsigned int j=0; j < edges_t1.size(); ++j)
		{
			HalfEdge &he = edges.at(edges_first_idx + j);
			TerrainHalfEdge &he_t1 = edges_t1.at(j);

			if (he_t1.getFace() == NULL) he.setFace(NULL);
			else he.setFace(&faces[0] + faces_first_idx 
				+ mesh_t1.getFaceId(he_t1.getFace()));

			he.setNext(&edges[0] + edges_first_idx
				+ mesh_t1.getHalfEdgeId(he_t1.getNext()));
			he.setOrigin(&vertices[0] + vertices_first_idx 
				+ mesh_t1.getVertexId(he_t1.getOrigin()));
			he.setPrev(&edges[0] + edges_first_idx
				+ mesh_t1.getHalfEdgeId(he_t1.getPrev()));
			he.setTwin(&edges[0] + edges_first_idx
				+ mesh_t1.getHalfEdgeId(he_t1.getTwin()));

			if (he.getData().edgeData == NULL) { // filtering twin edge that is already set

				// create an edgedata
				EdgeData ed;
				EdgeData::Source s;
				s.he = &he_t1;
				s.v = &v_t2;
				s.he_is_from_patch = false;
				ed.sources.push_back(s);
				edgeDataContainer.push_back(ed);
				he.getData().edgeData = &edgeDataContainer.back();
				he.getTwin()->getData().edgeData = &edgeDataContainer.back();
			}
			else { // when already set, we can get twin's origin.

				// considering halfEdge position (face is upward or downward), identify upedge and downedge.
				EdgeData *ed = he.getData().edgeData;
				if (he.getOrigin()->getData() < he.getTwin()->getOrigin()->getData()) {
					ed->halfEdge_up = &he;
					ed->halfEdge_down = he.getTwin();
				}
				else {
					ed->halfEdge_up = he.getTwin();
					ed->halfEdge_down = &he;
				}
			}
		}

		// faces will be computed later

		vertices_first_idx += vertices_t1.size();
		edges_first_idx += edges_t1.size();
	}
		
	// Insert all copied (-t2) with translation of (vertex of t1).
	for (unsigned int i=0; i < vertices_t1.size(); ++i)
	{
		TerrainVertex &v_t1 = vertices_t1.at(i);

		// copy vertex info
		for (unsigned int j=0; j < vertices_t2.size(); ++j)
		{
			Vertex &v = vertices[vertices_first_idx + j];
			TerrainVertex &v_t2 = vertices_t2.at(j);
			v.getData().x = v_t1.getData().p.x - v_t2.getData().p.x;
			v.getData().y = v_t1.getData().p.y - v_t2.getData().p.y;

			if (v_t2.getIncidentEdge() == NULL) v.setIncidentEdge(NULL);
			else v.setIncidentEdge(&edges[0] + edges_first_idx
				+ mesh_t1.getHalfEdgeId(v_t2.getIncidentEdge()));
		}
		// copy edge info
		for (unsigned int j=0; j < edges_t2.size(); ++j)
		{
			HalfEdge &he = edges[edges_first_idx + j];
			TerrainHalfEdge &he_t2 = edges_t2.at(j);

			if (he_t2.getFace() == NULL) he.setFace(NULL);
			else he.setFace(&faces[0] + faces_first_idx 
				+ mesh_t2.getFaceId(he_t2.getFace()));

			he.setNext(&edges[0] + edges_first_idx 
				+ mesh_t2.getHalfEdgeId(he_t2.getNext()));
			he.setOrigin(&vertices[0] + vertices_first_idx 
				+ mesh_t2.getVertexId(he_t2.getOrigin()));
			he.setPrev(&edges[0] + edges_first_idx
				+ mesh_t2.getHalfEdgeId(he_t2.getPrev()));
			he.setTwin(&edges[0] + edges_first_idx
				+ mesh_t2.getHalfEdgeId(he_t2.getTwin()));

			if (he.getData().edgeData == NULL) { // when edge data is not ready, we should create data.

				// create an edgedata
				EdgeData ed;
				EdgeData::Source s;
				s.he = &he_t2;
				s.v = &v_t1;
				s.he_is_from_patch = true;
				ed.sources.push_back(s);
				edgeDataContainer.push_back(ed);
				he.getData().edgeData = &edgeDataContainer.back();
				he.getTwin()->getData().edgeData = &edgeDataContainer.back();
			}
			else { // when already set, we can get twin's origin.

				// considering halfEdge position (face is upward or downward), identify upedge and downedge.
				EdgeData *ed = he.getData().edgeData;
				if (he.getOrigin()->getData() < he.getTwin()->getOrigin()->getData()) {
					ed->halfEdge_up = &he;
					ed->halfEdge_down = he.getTwin();
				}
				else {
					ed->halfEdge_up = he.getTwin();
					ed->halfEdge_down = &he;
				}
			}
		}
		
		// face will be computed later
		
		vertices_first_idx += vertices_t2.size();
		edges_first_idx += edges_t2.size();
		faces_first_idx += faces_t2.size();
	}

	// Test step and size.
	if (vertices_first_idx != vertices.size() &&
		edges_first_idx != edges.size() &&
		faces_first_idx != faces.size())
		throw cpp::Exception("Sum of steps is not equal to the total size.");

	// Test edgeDataContainer.
	if (edgeDataContainer.size() != number_of_halfedges() / 2)
		throw cpp::Exception("# of Edges is not (# of Halfedges)/2.");

	// Run sweepline algorithm.
	sweepLine.initialize(this);
	sweepLine.run();
}


bool SweepLine::EdgeDataCompare::operator()(const Arrangement::EdgeData *ed1, const Arrangement::EdgeData *ed2) const // return ed1 < ed2
{
	Arrangement::VertexData &vd1L = ed1->halfEdge_up->getOrigin()->getData();
	Arrangement::VertexData &vd1R = ed1->halfEdge_down->getOrigin()->getData();
	Arrangement::VertexData &vd2L = ed2->halfEdge_up->getOrigin()->getData();
	Arrangement::VertexData &vd2R = ed2->halfEdge_down->getOrigin()->getData();

	double sx = getX();
	
	// Compare y-coordinate of the intersection with the sweepLine.
	// (vd1R.y - vd1L.y) * (sx - vd1L.x) / (vd1R.x - vd1L.x) = y1 < y2 = (vd2R.y - vd2L.y) * (sx - vd2L.x) / (vd2R.x - vd2L.x);
	double comp_y = (vd1R.y - vd1L.y) * (sx - vd1L.x) * (vd2R.x - vd2L.x) - (vd2R.y - vd2L.y) * (sx - vd2L.x) * (vd1R.x - vd1L.x);

	if (comp_y == 0) {
		// Compare slope
		// (vd1R.y - vd1L.y) / (vd1R.x - vd1L.x) = s1 < s2 = (vd2R.y - vd2L.y) / (vd2R.x - vd2L.x)
		double comp_s = (vd1R.y - vd1L.y) * (vd2R.x - vd2L.x) - (vd2R.y - vd2L.y) * (vd1R.x - vd1L.x);

		//if (comp_s == 0) {
		//	return ed1 < ed2; // not comparable. just distinguish.
		//}
		//else
			return comp_s < 0;
	}
	else
		return comp_y < 0;
}

// ed2 was below ed1. ed2 will go up, and ed1 will go down relatively.
bool SweepLine::handleProperIntersectionEvent(Arrangement::EdgeData *ed1, Arrangement::EdgeData *ed2)
{
	// check source
	for (unsigned int i = 0; i < ed1->sources.size(); ++i) {
		for (unsigned int j = 0; j < ed2->sources.size(); ++j) {
			if (ed1->sources[i] == ed2->sources[j]) { // if the edges are from the same source, don't handle intersection event (already proved not to be intersect).
				return false;
			}
		}
	}

	Arrangement::VertexData &vd1L = ed1->halfEdge_up->getOrigin()->getData();
	Arrangement::VertexData &vd1R = ed1->halfEdge_down->getOrigin()->getData();
	Arrangement::VertexData &vd2L = ed2->halfEdge_up->getOrigin()->getData();
	Arrangement::VertexData &vd2R = ed2->halfEdge_down->getOrigin()->getData();

	// Ax + By = C
	double A_1 = vd1R.y - vd1L.y;
	double B_1 = vd1L.x - vd1R.x;
	double C_1 = A_1 * vd1L.x + B_1 * vd1L.y;
	double A_2 = vd2R.y - vd2L.y;
	double B_2 = vd2L.x - vd2R.x;
	double C_2 = A_2 * vd2L.x + B_2 * vd2L.y;

	// determinant
	double det = A_1 * B_2 - A_2 * B_1;
	//bool sameLine = false;
	if (det == 0) { // Lines are parallel
		//if (B_1 * C_2 == B_2 * C_1)
		//	sameLine = true;
		//else
			return false;
	}

	// Intersection point of two "lines"
	double x_det = (B_2 * C_1 - B_1 * C_2);
	double y_det = (A_1 * C_2 - A_2 * C_1);
	double max_y = vd1L.y > vd1R.y ? vd1L.y : vd1R.y;
	double min_y = vd1L.y > vd1R.y ? vd1R.y : vd1L.y;

	// (x,y) in the segment range?
	if ((det > 0 && vd1L.x * det < x_det && x_det < vd1R.x * det && min_y * det < y_det && y_det < max_y * det)
	 || (det < 0 && vd1L.x * det > x_det && x_det > vd1R.x * det && min_y * det > y_det && y_det > max_y * det)) { // if two edges properly intersect, create new intersection event.
		Arrangement::Vertex *v = updateDCELProperIntersection(ed1, ed2, x_det / det, y_det / det);
		EventPoint ep(v, EventPoint::VERTEXPOINT);
		events.push(ep);
		return true;
	}
	else { // if segments not intersect, return false.
		return false;
	}
}

// ed2 was below ed1. ed2 will go up, and ed1 will go down relatively.
// returns (ed1N, ed2N)
Arrangement::Vertex * 
SweepLine::updateDCELProperIntersection(Arrangement::EdgeData *ed1, Arrangement::EdgeData *ed2, double int_x, double int_y)
{
	/*
			he1u		 [f1u]		    he2Nu
		ed1----------¦¤			¦£------------ed2N (N for new)
			he1d	 ¦¦¦¤v_int¦£¦¥	    he2Nd
			  [f1d]	   ¦§-----¦©	 [f2u]
			he2u	 ¦£¦¥	  ¦¦¦¤	    he1Nu
		ed2----------¦¥			¦¦------------ed1N
			he2d		 [f2d]		    he1Nd
	*/

	Arrangement::HalfEdge *he1u = ed1->halfEdge_up;
	Arrangement::HalfEdge *he1d = ed1->halfEdge_down;
	Arrangement::HalfEdge *he2u = ed2->halfEdge_up;
	Arrangement::HalfEdge *he2d = ed2->halfEdge_down;

	// create a new intersection vertex
	Arrangement::Vertex *v_int = parent->createVertex();
	v_int->getData().x = int_x;
	v_int->getData().y = int_y;

	// create four halfedges for ed1 and ed2 after v
	Arrangement::HalfEdge *he1Nu = parent->createHalfEdge();
	Arrangement::HalfEdge *he1Nd = parent->createHalfEdge();
	Arrangement::HalfEdge *he2Nu = parent->createHalfEdge();
	Arrangement::HalfEdge *he2Nd = parent->createHalfEdge();

	// create two edgedata for ed1 and ed2 after v
	Arrangement::EdgeData *ed1N = parent->createEdgeData();
	Arrangement::EdgeData *ed2N = parent->createEdgeData();
	ed1N->sources = ed1->sources;
	ed1N->halfEdge_up = he1Nu;
	ed1N->halfEdge_down = he1Nd;
	ed2N->sources = ed2->sources;
	ed2N->halfEdge_up = he2Nu;
	ed2N->halfEdge_down = he2Nd;
	he1Nu->getData().edgeData = ed1N;
	he1Nd->getData().edgeData = ed1N;
	he2Nu->getData().edgeData = ed2N;
	he2Nd->getData().edgeData = ed2N;

	// store required old structures (because updating new structures change the old structures automatically)
	Arrangement::HalfEdge *he1u_next = he1u->getNext();
	Arrangement::HalfEdge *he1d_prev = he1d->getPrev();
	Arrangement::HalfEdge *he2u_next = he2u->getNext();
	Arrangement::HalfEdge *he2d_prev = he2d->getPrev();

	// update new structures
	v_int->setIncidentEdge(he2Nu);
	he2Nu->setNext(he2u_next);
	he2Nu->setOrigin(v_int);
	he2Nu->setPrev(he1u);
	he2Nu->setTwin(he2Nd);
	he2Nd->setNext(he1Nu);
	he2Nd->setOrigin(he2d->getOrigin());
	he2Nd->setPrev(he2d_prev);
	//he2Nd->setTwin(he2Nu);
	he1Nu->setNext(he1u_next);
	he1Nu->setOrigin(v_int);
	//he1Nu->setPrev(he2Nd);
	he1Nu->setTwin(he1Nd);
	he1Nd->setNext(he2d);
	he1Nd->setOrigin(he1d->getOrigin());
	he1Nd->setPrev(he1d_prev);
	//he1Nd->setTwin(he1Nu);

	// update old structures
	//he1u->setNext(he2Nu);
	he1d->setOrigin(v_int);
	he1d->setPrev(he2u);
	//he2u->setNext(he1d);
	he2d->setOrigin(v_int);
	//he2d->setPrev(he1Nd);

	// face will be handled later

	// returns (ed1N, ed2N)
	return v_int;
}

void SweepLine::updateDCELVertexEdgeIntersection(Arrangement::Vertex *v, Arrangement::EdgeData *ed)
{
	/*
			  v_dummy			   v¦¢				 edu  v¦¢  edNu
		(ed)-----¡Ü-----(edN)	+	¡Ü	  =		(ed)-------¡Ü-------(edN)
									¦¢				 edNd  ¦¢  edd
	*/

	// existing structure of ed
	Arrangement::HalfEdge *edu = ed->halfEdge_up;
	Arrangement::HalfEdge *edd = ed->halfEdge_down;

	// create new structure
	Arrangement::EdgeData *edN = parent->createEdgeData();
	Arrangement::HalfEdge *edNu = parent->createHalfEdge();
	Arrangement::HalfEdge *edNd = parent->createHalfEdge();
	Arrangement::Vertex *v_dummy = parent->createVertex();
	v_dummy->getData() = v->getData(); // copy constructor

	// link edges, halfedges, v_dummy
	edN->halfEdge_up = edNu;
	edN->halfEdge_down = edd;
	ed->halfEdge_down = edNd;
	edNu->setNext(edu->getNext());
	edNu->setOrigin(v_dummy);
	edNu->setPrev(edu);
	edNu->setTwin(edd);
	edNd->setNext(edd->getNext());
	edNd->setOrigin(v_dummy);
	edNd->setPrev(edd);
	edNd->setTwin(edu);

	// merge v_dummy to v
	updateDCEL2VertexIntersection(v, v_dummy);
}

void SweepLine::updateDCEL2VertexIntersection(Arrangement::Vertex *v, Arrangement::Vertex *v_del)
{
	Vector2<double> vec_v(v->getData());
	Vector2<double> vec_v_del(v_del->getData());

	// check the same point
	if (!(v->getData() == v_del->getData())) throw cpp::Exception("v and v_del are not the same point.");
	
	// merge incident edges of v_del into v
	/*
		¦£------- he_prev (from v)
		¦¢
		¦§------- he (from v_del)
		¦¢
		¦¦------- he_next (from v)
	*/

	// initialize edges from v
	Arrangement::EdgeIterator eit(v);
	Arrangement::HalfEdge *he_prev = eit.getNext();
	Vector2<double> vec_he_prev = Vector2<double>(he_prev->getTwin()->getOrigin()->getData()) - vec_v;
	Arrangement::HalfEdge *he_next = eit.getNext();
	Vector2<double> vec_he_next = Vector2<double>(he_next->getTwin()->getOrigin()->getData()) - vec_v;

	// traverse incident edges of v_del and attach them to v
	Arrangement::EdgeIterator eit_del(v_del);
	while (eit_del.hasNext()) {
		Arrangement::HalfEdge *he = eit_del.getNext();
		Vector2<double> vec_he = Vector2<double>(he->getTwin()->getOrigin()->getData()) - vec_v_del;

		// traverse incident edges of v while he becomes between he_prev and he_next.
		bool inserted = false;
		while (!inserted) {
			if (vec_he_prev * vec_he_next < 0) { // reflex angle
				if (!vec_he.isLeftFrom(vec_he_prev) || vec_he.isLeftFrom(vec_he_next)) { // if he is in prev~next
					// link edges of v_del to v
					he_prev->getTwin()->setNext(he);
					he_next->setPrev(he);
					inserted = true;
				}
			}
			else { // not reflex angle
				if (!vec_he.isLeftFrom(vec_he_prev) && vec_he.isLeftFrom(vec_he_next)) { // if he is in prev~next
					// link edges of v_del to v
					he_prev->getTwin()->setNext(he);
					he_next->setPrev(he);
					inserted = true;
				}
			}
			// get next one
			if (!inserted) {
				if (!eit.hasNext()) eit.reset();
				he_prev = he_next;
				vec_he_prev = vec_he_next;
				he_next = eit.getNext();
				vec_he_next = Vector2<double>(he_next->getTwin()->getOrigin()->getData()) - vec_v;
			}
		}

		// attach he to v
		he->setOrigin(v);
	}

	// erase v_del (lazy)
	unsigned int id = v_del - &(parent->vertices[0]);
	parent->deleteVertex(id);
}

void SweepLine::initialize(Arrangement *_parent)
{
	parent = _parent;

	// Initialize event queue
	events = EventQueue();

	// Insert event points related to the vertices
	for (unsigned int i = 0; i < parent->vertices.size(); ++i)
	{
		Arrangement::Vertex *v = &parent->vertices.at(i);
		events.push(EventPoint(v, EventPoint::VERTEXPOINT));
	}
}

void SweepLine::advance()
{
	EventPoint ep = events.top();
	events.pop();

	if (ep.state == EventPoint::VERTEXPOINT)
	{
		// Insert the edges before current sweepline.
		Arrangement::EdgeIterator eit(ep.v);

		// Find the first before-event edge and the first after-event edge.
		Arrangement::HalfEdge *BEedge, *AEedge;
		bool he_is_before_event = true;
		bool he_is_after_event = true;
		while (eit.hasNext()) {
			Arrangement::HalfEdge *he = eit.getNext();
			bool cur_he_is_before_event = he->getTwin()->getOrigin()->getData() < ep.v->getData();
			if (!he_is_before_event && he_is_after_event && cur_he_is_before_event) { // if he is the first before-event edge, set it to the first edge of eit.
				BEedge = he;
				break; // when this break does not accur, original ep.v was the first before-event edge.
			}
			else { // else, update before-event flags.
				he_is_before_event = cur_he_is_before_event;
				he_is_after_event = !cur_he_is_before_event;
			}
		}
		ep.v->setIncidentEdge(BEedge);

		// Iterator from before-edges to after-edges
		Arrangement::EdgeIterator BE2AEeit(ep.v);

		// Erase mode.
		while (BE2AEeit.hasNext()) {
			Arrangement::HalfEdge *he = BE2AEeit.getNext();
			if (he == AEedge) break; // if it becomes after-event edge, break.

			EdgeDataBBTIterator bbt_entry = edgeDataBBT.find(he->getData().edgeData);
			bbt_entry = edgeDataBBT.erase(bbt_entry);
			if (bbt_entry != edgeDataBBT.end() || bbt_entry != edgeDataBBT.begin()) { // if erased edge is not the highest of lowest, check intersection event.
				Arrangement::EdgeData *ed1 = *bbt_entry;
				Arrangement::EdgeData *ed2 = *--bbt_entry;
				handleProperIntersectionEvent(ed1, ed2);
			}
		}

		// Insert mode.
		EdgeDataBBTIterator hintit = edgeDataBBT.begin();
		while (BE2AEeit.hasNext()) {
			Arrangement::HalfEdge *he = BE2AEeit.getNext();

			EdgeDataBBTIterator hintit_next = edgeDataBBT.insert(hintit, he->getData().edgeData);
			if (hintit_next == hintit) { // if there is already an edge of the same slope (so cannot be inserted), make it twin-edge. -----¡Ü======¡Ü-----

			}
			else { // else, check intersections with neighbors of BBT
				EdgeDataBBTIterator it = hintit_next;
				if (it != edgeDataBBT.begin()) { // if it is not the lowest, check intersection with the lower neighbor.
					EdgeDataBBTIterator it_low = it;
					handleProperIntersectionEvent(*--it_low, *it);
				}
				if (it != --edgeDataBBT.end()) { // if it is not the highest, check intersection with the higher neighbor.
					EdgeDataBBTIterator it_high = it;
					handleProperIntersectionEvent(*++it_high, *it);
				}
				hintit = hintit_next; // update hint.
			}
		}
	}
	else if (ep.state == EventPoint::CROSSINGPOINT)
	{
		// Swap contents of BBT
		ep.v
	}
	else
	{
		throw cpp::Exception("Error event.");
	}

	EdgeDataBBTIterator it_debug = edgeDataBBT.begin();
	while (it_debug != edgeDataBBT.end()) {
		std::cerr << " = ";
		(*it_debug)->print();
		std::cerr << '\n';
		++it_debug;
	}

	++eventCount;
}
