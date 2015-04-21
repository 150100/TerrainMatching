#include "arrangement.h"

#include <algorithm>

std::vector<Arrangement::Vertex> Arrangement::vertices;
std::vector<Arrangement::HalfEdge> Arrangement::edges;
std::vector<Arrangement::Face> Arrangement::faces;
std::vector<Arrangement::EdgeData> Arrangement::edgeDataContainer;

std::queue<unsigned int> Arrangement::erasedVerticesIndices;
std::queue<unsigned int> Arrangement::erasedEdgesIndices;
std::queue<unsigned int> Arrangement::erasedFacesIndices;
std::queue<unsigned int> Arrangement::erasedEdgeDataContainerIndices;

Arrangement::SweepLine::EventQueue Arrangement::SweepLine::events;
Arrangement::SweepLine::EdgeDataBBT Arrangement::SweepLine::edgeDataBBT;
int Arrangement::SweepLine::eventCount = 0;
Arrangement::SweepLine Arrangement::sweepLine;

template <class T>
class Vector2
{
private:
	double x, y;

public:
	Vector2() {}
	Vector2(double _x, double _y) : x(_x), y(_y) {}
	Vector2(ArrangementVertexData &vd) : x(vd.x), y(vd.y) {}

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

	vertices.reserve(vertices_t1.size() * vertices_t2.size() * vertices_t2.size() * 10); // nm^2 estimation
	edges.reserve(vertices_t1.size() * vertices_t2.size() * vertices_t2.size() * 60); // nm^2 estimation
	faces.reserve(vertices_t1.size() * vertices_t2.size() * vertices_t2.size() * 20); // nm^2 estimation
	edgeDataContainer.reserve(edges.capacity() / 2);

	unsigned int vertices_first_idx = 0;
	unsigned int edges_first_idx = 0;
	unsigned int faces_first_idx = 0;

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
		// copy face info
		for (unsigned int j=0; j < faces_t1.size(); ++j)
		{
			Face &f = faces[faces_first_idx + j];
			TerrainFace &f_t1 = faces_t1.at(j);

			if (f_t1.getBoundary() == NULL) f.setBoundary(NULL);
			else f.setBoundary(&edges[0] + edges_first_idx
				+ mesh_t1.getHalfEdgeId(f_t1.getBoundary()));
		}

		vertices_first_idx += vertices_t1.size();
		edges_first_idx += edges_t1.size();
		faces_first_idx += faces_t1.size();
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
		// copy face info
		for (unsigned int j=0; j < faces_t2.size(); ++j)
		{
			Face &f = faces[faces_first_idx + j];
			TerrainFace &f_t2 = faces_t2.at(j);

			if (f_t2.getBoundary() == NULL) f.setBoundary(NULL);
			else f.setBoundary(&edges[0] + edges_first_idx
				+ mesh_t2.getHalfEdgeId(f_t2.getBoundary()));
		}
		
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

	/// Sweepline Algorithm...
	sweepLine.initialize();
	sweepLine.run();
}

// erase edge(with halfedges) itself and right vertex (= halfedge_down->origin).
// merge them into left vertex (= halfedge_up->origin).
// make sure the edge length is 0.
void Arrangement::eraseZeroLengthEdge(EdgeData *ed)
{
	Vertex *v = ed->halfEdge_up->getOrigin();
	Vector2<double> vec_v(v->getData());
	Vertex *v_del = ed->halfEdge_down->getOrigin();
	Vector2<double> vec_v_del(v_del->getData());

	// link edges
	ed->halfEdge_up->getPrev()->setNext(ed->halfEdge_up->getNext());
	ed->halfEdge_down->getPrev()->setNext(ed->halfEdge_down->getNext());
	v->setIncidentEdge(ed->halfEdge_up->getNext());
	
	// merge incident edges of v_del into v
	/*
		   ¦£-¡æhe_prev (from v)
		¦£-¦¥(face from v)
		¦«----¡æhe (from v_del)
		¦¦-¦¤(face from v_del)
		   ¦¦-¡æhe_next (from v)
	*/

	EdgeIterator eit(v);

	std::vector<HalfEdge *> hes_before;
	while (eit.hasNext()) {
		HalfEdge *he = eit.getNext();
		hes_before.push_back(he);
	}
	eit.reset();

	HalfEdge *he_prev = eit.getNext();
	Vector2<double> vec_he_prev = Vector2<double>(he_prev->getTwin()->getOrigin()->getData()) - vec_v;
	HalfEdge *he_next;
	// get next one
	if (eit.hasNext()) {
		he_next = eit.getNext();
	}
	else {
		eit.reset();
		he_next = eit.getNext();
	}
	EdgeIterator eit_del(v_del);
	std::vector<HalfEdge *> hes_del;
	while (eit_del.hasNext()) {
		HalfEdge *he = eit_del.getNext();
		hes_del.push_back(he);
	}
	for (unsigned int i = 0; i < hes_del.size(); ++i) {
		HalfEdge *he = hes_del[i];
		
		Vector2<double> vec_he = Vector2<double>(he->getTwin()->getOrigin()->getData()) - vec_v_del;
		Vector2<double> vec_he_next = Vector2<double>(he_next->getTwin()->getOrigin()->getData()) - vec_v;
		if (vec_he == Vector2<double>(0,0))
			continue; // bypass
		else if (vec_he_prev * vec_he_next < 0) { // reflex angle
			while (!(!vec_he.isLeftFrom(vec_he_prev) || vec_he.isLeftFrom(vec_he_next))) { // if he is not in prev~next
				// proceed eit
				// store previous one
				he_prev = he_next;
				vec_he_prev = vec_he_next;
				// get next one
				if (eit.hasNext()) {
					he_next = eit.getNext();
				}
				else {
					eit.reset();
					he_next = eit.getNext();
				}
			}
			// link edges of v_del to v
			he_prev->getTwin()->setNext(he);
			he_next->setPrev(he);
			// link faces of v_del to v
			he->setFace(he_prev->getTwin()->getFace());
			he_next->setFace(he->getTwin()->getFace());
		}
		else { // not reflex angle
			while (! (!vec_he.isLeftFrom(vec_he_prev) && vec_he.isLeftFrom(vec_he_next))) { // if he is not in prev~next
				// proceed eit
				// store previous one
				he_prev = he_next;
				vec_he_prev = vec_he_next;
				// get next one
				if (eit.hasNext()) {
					he_next = eit.getNext();
				}
				else {
					eit.reset();
					he_next = eit.getNext();
				}
			}
			// link edges of v_del to v
			he_prev->getTwin()->setNext(he);
			he_next->setPrev(he);
			// link faces of v_del to v
			he->setFace(he_prev->getTwin()->getFace());
			he_next->setFace(he->getTwin()->getFace());
		}

	}

	std::vector<HalfEdge *> hes_after;
	eit.reset();
	while (eit.hasNext()) {
		HalfEdge *he = eit.getNext();
		hes_after.push_back(he);
	}

	// erase v (lazy)
	unsigned int id = v_del - &(vertices[0]);
	deleteVertex(id);
}

bool Arrangement::SweepLine::EdgeDataCompare::operator()(const EdgeData *ed1, const EdgeData *ed2) const // return ed1 < ed2
{
	VertexData &vd1L = ed1->halfEdge_up->getOrigin()->getData();
	VertexData &vd1R = ed1->halfEdge_down->getOrigin()->getData();
	VertexData &vd2L = ed2->halfEdge_up->getOrigin()->getData();
	VertexData &vd2R = ed2->halfEdge_down->getOrigin()->getData();

	double sx = sweepLine.getX();
	
	// Compare y-coordinate of the intersection with the sweepLine.
	// (vd1R.y - vd1L.y) * (sx - vd1L.x) / (vd1R.x - vd1L.x) = y1 < y2 = (vd2R.y - vd2L.y) * (sx - vd2L.x) / (vd2R.x - vd2L.x);
	double comp_y = (vd1R.y - vd1L.y) * (sx - vd1L.x) * (vd2R.x - vd2L.x) - (vd2R.y - vd2L.y) * (sx - vd2L.x) * (vd1R.x - vd1L.x);

	if (comp_y == 0) {
		// Compare slope
		// (vd1R.y - vd1L.y) / (vd1R.x - vd1L.x) = s1 < s2 = (vd2R.y - vd2L.y) / (vd2R.x - vd2L.x)
		double comp_s = (vd1R.y - vd1L.y) * (vd2R.x - vd2L.x) - (vd2R.y - vd2L.y) * (vd1R.x - vd1L.x);

		if (comp_s == 0) {
			return ed1 < ed2; // not comparable. just distinguish.
		}
		else
			return comp_s < 0;
	}
	else
		return comp_y < 0;
}

// ed2 was below ed1. ed2 will go up, and ed1 will go down relatively.
bool Arrangement::SweepLine::handleIntersectionEventWithDCEL(EdgeData *ed1, EdgeData *ed2)
{
	VertexData &vd1L = ed1->halfEdge_up->getOrigin()->getData();
	VertexData &vd1R = ed1->halfEdge_down->getOrigin()->getData();
	VertexData &vd2L = ed2->halfEdge_up->getOrigin()->getData();
	VertexData &vd2R = ed2->halfEdge_down->getOrigin()->getData();

	// Ax + By = C
	double A_1 = vd1R.y - vd1L.y;
	double B_1 = vd1L.x - vd1R.x;
	double C_1 = A_1 * vd1L.x + B_1 * vd1L.y;
	double A_2 = vd2R.y - vd2L.y;
	double B_2 = vd2L.x - vd2R.x;
	double C_2 = A_2 * vd2L.x + B_2 * vd2L.y;

	// line intersection
	double det = A_1 * B_2 - A_2 * B_1;
	bool sameLine = false;
	if (det == 0) { // Lines are parallel
		if (B_1 * C_2 == B_2 * C_1)
			sameLine = true;
		else
			return false;
	}

	// Intersection point of two "lines"
	double x_det = (B_2 * C_1 - B_1 * C_2);
	double y_det = (A_1 * C_2 - A_2 * C_1);
	double max_y = vd1L.y > vd1R.y ? vd1L.y : vd1R.y;
	double min_y = vd1L.y > vd1R.y ? vd1R.y : vd1L.y;

	// (x,y) in the segment range?
	if ((det > 0 && vd1L.x * det <= x_det && x_det <= vd1R.x * det && min_y * det <= y_det && y_det <= max_y * det)
	 || (det < 0 && vd1L.x * det >= x_det && x_det >= vd1R.x * det && min_y * det >= y_det && y_det >= max_y * det)) 
	{ // Segments intersect!
		updateIntersectionDCEL(ed1, ed2, x_det / det, y_det / det);
		return true;
	}
	else if (sameLine)
	{ // two edges are intersecting but parallel. Handle as two intersections, and remove 2-edge-face.
		
		/*
								   edR
			edL		------------------
			------------------     

						v
						v
				(first intersection)
						v
						v

					   edRN
			edL		------------------
			-------¡Ü f2edge
					----------
					   edLN

						v
						v
				(second intersection)
						v
						v

					   edRN
			edL		----------	 edRNN	
			-------¡Ü f2edge ¡Ü-------	     
					----------
					   edLN

					    v
					    v
				(merge edRN into edLN)
						v
						v
				
				       edLN
			-------¡Ü========¡Ü-------


			(edR and edLNN become 0-length edges)

		*/
		EdgeData *edL, *edR;
		if (vd1L.x < vd2L.x) { edL = ed1; edR = ed2; }
		else { edL = ed2; edR = ed1; }

		// first intersection
		VertexData &int1 = edR->halfEdge_up->getOrigin()->getData();
		auto edN = updateIntersectionDCEL(ed1, ed2, int1.x, int1.y);
		EdgeData *edLN = edN.first;
		EdgeData *edRN = edN.second;

		// second intersection
		VertexData &int2 = edLN->halfEdge_down->getOrigin()->getData();
		auto edNN = updateIntersectionDCEL(edRN, edLN, int2.x, int2.y);
		Face *f2edge = edLN->halfEdge_up->getFace();
		if (f2edge != edRN->halfEdge_down->getFace())
			throw cpp::Exception("f2edge error : not a same face of edRN and edLN");
		if (f2edge->getBoundary()->getNext()->getNext() != f2edge->getBoundary())
			throw cpp::Exception("f2edge error : not a 2-edge-face.");

		// merge edRN into edLN
		edRN->halfEdge_up->getFace()->setBoundary(edLN->halfEdge_up);
		edLN->halfEdge_up->setFace(edRN->halfEdge_up->getFace());
		edLN->sources.insert(edLN->sources.end(), edRN->sources.begin(), edRN->sources.end());
		unsigned int heRNu_id = edRN->halfEdge_up - &edges[0];
		unsigned int heRNd_id = edRN->halfEdge_down - &edges[0];
		unsigned int f2edge_id = f2edge - &faces[0];
		unsigned int edRN_id = edRN - &edgeDataContainer[0];
		deleteHalfEdge(heRNu_id);
		deleteHalfEdge(heRNd_id);
		deleteFace(f2edge_id);
		deleteEdgeData(edRN_id);
	}
	else
	{ // segments not intersect.
		return false;
	}
}

// ed2 was below ed1. ed2 will go up, and ed1 will go down relatively.
// returns (ed1N, ed2N)
std::pair<Arrangement::EdgeData *, Arrangement::EdgeData *> 
Arrangement::SweepLine::updateIntersectionDCEL(EdgeData *ed1, EdgeData *ed2, double int_x, double int_y)
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

	HalfEdge *he1u = ed1->halfEdge_up;
	HalfEdge *he1d = ed1->halfEdge_down;
	HalfEdge *he2u = ed2->halfEdge_up;
	HalfEdge *he2d = ed2->halfEdge_down;

	// lazy handling of intersections before related to faces
	he1u->setFace(he1u->getPrev()->getFace());
	he1d->setFace(he1d->getNext()->getFace());
	he2u->setFace(he2u->getPrev()->getFace());
	he2d->setFace(he2d->getNext()->getFace());

	// create a new intersection vertex
	Vertex *v_int = createVertex();
	v_int->getData().x = int_x;
	v_int->getData().y = int_y;

	// create four halfedges for ed1 and ed2 after v
	HalfEdge *he1Nu = createHalfEdge();
	HalfEdge *he1Nd = createHalfEdge();
	HalfEdge *he2Nu = createHalfEdge();
	HalfEdge *he2Nd = createHalfEdge();

	// create two edgedata for ed1 and ed2 after v
	EdgeData *ed1N = createEdgeData();
	EdgeData *ed2N = createEdgeData();
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

	// faces to be modified
	Face *f1u = he1u->getFace();
	Face *f1d = he1d->getFace();
	Face *f2u = he2u->getFace();
	Face *f2d = he2d->getFace();

	// store required old structures (because updating new structures change the old structures automatically)
	HalfEdge *he1u_next = he1u->getNext();
	HalfEdge *he1d_prev = he1d->getPrev();
	HalfEdge *he2u_next = he2u->getNext();
	HalfEdge *he2d_prev = he2d->getPrev();

	// update new structures
	v_int->setIncidentEdge(he2Nu);
	he2Nu->setFace(f1u);
	he2Nu->setNext(he2u_next);
	he2Nu->setOrigin(v_int);
	he2Nu->setPrev(he1u);
	he2Nu->setTwin(he2Nd);
	he2Nd->setFace(f2u);
	he2Nd->setNext(he1Nu);
	he2Nd->setOrigin(he2d->getOrigin());
	he2Nd->setPrev(he2d_prev);
	//he2Nd->setTwin(he2Nu);
	he1Nu->setFace(f2u);
	he1Nu->setNext(he1u_next);
	he1Nu->setOrigin(v_int);
	//he1Nu->setPrev(he2Nd);
	he1Nu->setTwin(he1Nd);
	he1Nd->setFace(f2d);
	he1Nd->setNext(he2d);
	he1Nd->setOrigin(he1d->getOrigin());
	he1Nd->setPrev(he1d_prev);
	//he1Nd->setTwin(he1Nu);

	// update old structures
	//he1u->setNext(he2Nu);
	he1d->setOrigin(v_int);
	he1d->setPrev(he2u);
	he2u->setFace(f1d);
	//he2u->setNext(he1d);
	he2d->setOrigin(v_int);
	//he2d->setPrev(he1Nd);

	// delete zero-distance edges (possibly ed1, ed2) from BBT, and merge DCEL structure.
	if (he1u->getOrigin()->getData() == he1d->getOrigin()->getData()) {
		eraseZeroLengthEdge(ed1);
	}
	if (he2u->getOrigin()->getData() == he2d->getOrigin()->getData()) {
		eraseZeroLengthEdge(ed2);
	}
	if (he1Nu->getOrigin()->getData() == he1Nd->getOrigin()->getData()) {
		eraseZeroLengthEdge(ed1N);
	}
	if (he2Nu->getOrigin()->getData() == he2Nd->getOrigin()->getData()) {
		eraseZeroLengthEdge(ed2N);
	}

	// update face structures
	if (f1u != NULL) {
		f1u->setBoundary(he1u);
		//EdgeIterator eit(f1u);
		//while (eit.hasNext()) {
		//	HalfEdge *he = eit.getNext();
		//	he->setFace(f1u);
		//}
	}
	if (f1d != NULL) {
		f1d->setBoundary(he1d);
		//EdgeIterator eit(f1d);
		//while (eit.hasNext()) {
		//	HalfEdge *he = eit.getNext();
		//	he->setFace(f1d);
		//}
	}
	if (f2u != NULL) {
		f2u->setBoundary(he2Nd);
		//EdgeIterator eit(f2u);
		//while (eit.hasNext()) {
		//	HalfEdge *he = eit.getNext();
		//	he->setFace(f2u);
		//}
	}
	if (f2d != NULL) {
		f2d->setBoundary(he2d);
		//EdgeIterator eit(f2d);
		//while (eit.hasNext()) {
		//	HalfEdge *he = eit.getNext();
		//	he->setFace(f2d);
		//}
	}

	// create new intersection event
	EventPoint ep(ed1, ed2, ed1N, ed2N, v_int, EventPoint::CROSSING);
	events.push(ep);

	// returns (ed1N, ed2N)
	return std::pair<EdgeData *, EdgeData *>(ed1N, ed2N);
}

void Arrangement::SweepLine::initialize()
{
	// Insert event points related to the end points of edges
	for (unsigned int i = 0; i < edgeDataContainer.size(); ++i)
	{
		EdgeData *ed = &edgeDataContainer.at(i);
		events.push(EventPoint(ed, EventPoint::STARTPOINT));
	}

	edgeDataBBT.insert(events.top().ed1);
	events.pop();
}

void Arrangement::SweepLine::advance()
{
	// Split all the edges, and merge them. (inefficient)
	EventPoint ep = events.top();
	events.pop();
	++eventCount;

	if (eventCount % 1 == 100) {
		std::cerr << "Event " << eventCount << "\n";
		std::cerr << "# of edges in BBT : " << edgeDataBBT.size() << "\n";
	}

	if (ep.state == EventPoint::STARTPOINT)
	{
		// insert the edge into BBT
		EdgeDataBBTIterator it = edgeDataBBT.insert(ep.ed1);
		std::cerr << "Insert : ";
		ep.ed1->print();
		std::cerr << '\n';
		
		// check future intersections (with adjacent edges) & insert them into event
		EdgeDataBBTIterator itPrev, itNext;

		itPrev = it;
		if (itPrev != edgeDataBBT.begin()) {
			--itPrev;
			handleIntersectionEventWithDCEL(*itPrev, ep.ed1);
		}

		itNext = it;
		++itNext;
		if (itNext != edgeDataBBT.end()) {
			handleIntersectionEventWithDCEL(ep.ed1, *itNext);
		}

		// insert future deletion
		events.push(EventPoint(ep.ed1, EventPoint::ENDPOINT));
	}
	else if (ep.state == EventPoint::ENDPOINT) 
	{
		unsigned int erase_num = edgeDataBBT.erase(ep.ed1);
		std::cerr << "Erase : ";
		ep.ed1->print();
		std::cerr << '\n';
	}
	else if (ep.state == EventPoint::CROSSING)
	{
		edgeDataBBT.erase(ep.ed1);
		edgeDataBBT.erase(ep.ed2);
		edgeDataBBT.insert(ep.ed1N);
		edgeDataBBT.insert(ep.ed2N);
		std::cerr << "Erase : ";
		ep.ed1->print();
		std::cerr << '\n';
		std::cerr << "Erase : ";
		ep.ed2->print();
		std::cerr << '\n';
		std::cerr << "Insert : ";
		ep.ed1N->print();
		std::cerr << '\n';
		std::cerr << "Insert : ";
		ep.ed2N->print();
		std::cerr << '\n';

		// insert future deletion
		events.push(EventPoint(ep.ed1N, EventPoint::ENDPOINT));
		events.push(EventPoint(ep.ed2N, EventPoint::ENDPOINT));
	}
	else
		throw cpp::Exception("Invalid event state.");

	EdgeDataBBTIterator it_debug = edgeDataBBT.begin();
	while (it_debug != edgeDataBBT.end()) {
		std::cerr << " = ";
		(*it_debug)->print();
		std::cerr << '\n';
		++it_debug;
	}
}
