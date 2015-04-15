#include "arrangement.h"

#include <algorithm>

std::vector<Arrangement::Vertex> Arrangement::vertices;
std::vector<Arrangement::HalfEdge> Arrangement::edges;
std::vector<Arrangement::Face> Arrangement::faces;
std::vector<Arrangement::EdgeData> Arrangement::edgeDataContainer;

Arrangement::SweepLine::EventQueue Arrangement::SweepLine::events;
Arrangement::SweepLine::EdgeDataBBT Arrangement::SweepLine::edgeDataBBT;
int Arrangement::SweepLine::eventCount = 0;
Arrangement::SweepLine Arrangement::sweepLine;


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
		// Compare slopes
		// (vd1R.y - vd1L.y) / (vd1R.x - vd1L.x) = s1 < s2 = (vd2R.y - vd2L.y) / (vd2R.x - vd2L.x)
		double comp_s = (vd1R.y - vd1L.y) * (vd2R.x - vd2L.x) - (vd2R.y - vd2L.y) * (vd1R.x - vd1L.x);

		if (comp_s == 0) {

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

	//// Same endpoint
	//if (vd1L == vd2L) {
	//	Vertex *v1 = ed1->halfEdge_up->getOrigin();
	//	Vertex *v2 = ed2->halfEdge_up->getOrigin();
	//	if (v1 == v2) return false; // two edges are from the same terrain, so already connected.
	//	else ed2->halfEdge_up->setOrigin(v1); // merge vertex to v1

	//	Face *f1u = ed1->halfEdge_up->getFace();
	//	Face *f1d = ed1->halfEdge_down->getFace();
	//	Face *f2u = ed2->halfEdge_up->getFace();
	//	Face *f2d = ed2->halfEdge_down->getFace();
	//	f2u->setBoundary()

	//	ed2->halfEdge_up->setPrev(ed1->halfEdge_down);
	//	ed2->halfEdge_up->setFace
	//}

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

	// (x,y) in the segment range? (excluding endpoints)
	if ((det > 0 && vd1L.x * det <= x_det && x_det <= vd1R.x * det && min_y * det <= y_det && y_det <= max_y * det)
	 || (det < 0 && vd1L.x * det >= x_det && x_det >= vd1R.x * det && min_y * det >= y_det && y_det >= max_y * det)) { // Segments intersect!


		/*
				he1u		 [f1u]		    he2Nu
			ed1----------��			��------------ed2N (N for new)
				he1d	 ����v_int����	    he2Nd
				  [f1d]	   ��-----��	 [f2u]
				he2u	 ����	  ����	    he1Nu
			ed2----------��			��------------ed1N
				he2d		 [f2d]		    he1Nd
		*/

		HalfEdge *he1u = ed1->halfEdge_up;
		HalfEdge *he1d = ed1->halfEdge_down;
		HalfEdge *he2u = ed2->halfEdge_up;
		HalfEdge *he2d = ed2->halfEdge_down;

		// create a new intersection vertex
		vertices.push_back( Vertex() );
		Vertex *v_int = &vertices.back();
		v_int->getData().x = x_det / det;
		v_int->getData().y = y_det / det;

		// create four halfedges for ed1 and ed2 after v
		edges.insert(edges.end(), 4, HalfEdge());
		auto heit = edges.end();
		HalfEdge *he1Nu = &*(--heit);
		HalfEdge *he1Nd = &*(--heit);
		HalfEdge *he2Nu = &*(--heit);
		HalfEdge *he2Nd = &*(--heit);

		// create two edgedata for ed1 and ed2 after v
		edgeDataContainer.insert(edgeDataContainer.end(), 2, EdgeData());
		auto edit = edgeDataContainer.end();
		EdgeData *ed1N = &*(--edit);
		EdgeData *ed2N = &*(--edit);
		ed1N->sources = ed1->sources;
		ed1N->halfEdge_up = he1Nu;
		ed1N->halfEdge_down = he1Nd;
		ed2N->sources = ed2->sources;
		ed2N->halfEdge_up = he2Nu;
		ed2N->halfEdge_down = he2Nd;
		he1Nu->getData().edgeData = ed1N;
		he2Nu->getData().edgeData = ed2N;

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
		if (f1u != NULL) f1u->setBoundary(he1u);
		if (f1d != NULL) f1d->setBoundary(he1d);
		if (f2u != NULL) f2u->setBoundary(he2Nd);
		if (f2d != NULL) f2d->setBoundary(he2d);

		// create new intersection event
		EventPoint ep(ed1, ed2, ed1N, ed2N, v_int, EventPoint::CROSSING);
		events.push(ep);
		return true;
	}
	else { // segments not intersect.
		return false;
	}
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
		edgeDataBBT.erase(ep.ed1);
		std::cerr << "Erase : ";
		ep.ed1->print();
		std::cerr << '\n';
	}
	else if (ep.state == EventPoint::CROSSING)
	{
		//edgeDataBBT.erase(ep.ed1);
		//edgeDataBBT.erase(ep.ed2);
		edgeDataBBT.insert(ep.ed1N);
		edgeDataBBT.insert(ep.ed2N);
		//std::cerr << "Erase : ";
		//ep.ed1->print();
		//std::cerr << '\n';
		//std::cerr << "Erase : ";
		//ep.ed2->print();
		//std::cerr << '\n';
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

	++eventCount;
}
