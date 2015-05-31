#include "translationspacesubdivision.h"

#include "DCEL/DCELStream.h"

TranslationSpaceSubdivision::TranslationSpaceSubdivision(TerrainWithGrids *_t1, TerrainWithGrids *_t2)
	: t1(_t1), t2(_t2), arr(_t1, _t2)
{
    // initialize search data and CS
    init();
}

void TranslationSpaceSubdivision::init()
{
    // find a boundary halfedge incident to the unbounded face.
	Arrangement::HalfEdge* he = arr.getFirstHalfEdge();

	// test of he
	if (he == NULL) throw cpp::Exception("There is no outer-face.");
	if (he->getTwin()->getFace() == NULL) throw cpp::Exception("First face is NULL. (TSS::init())");
	if (he->getFace() != arr.getOuterface()) throw cpp::Exception("Outer face is not outer. (TSS::init())");

    // initialize first cell
    update_CS(he);
    m_he_stack.push(he);
    m_he_path.push(he);
}

DFSState TranslationSpaceSubdivision::advanceDFS()
{
    if (m_he_stack.empty()) throw cpp::Exception("TSS stack should not be empty.");
	
    // previous target
	Arrangement::HalfEdge* cur_he = m_he_stack.top(); // enterence to the target face
	Arrangement::Face* cur_f = cur_he->getTwin()->getFace(); // target face

    // update the previous one
	cur_f->getData().state = Arrangement::FaceData::PASSED;
    m_he_stack.pop();

    // push adjacent faces to stack
	Arrangement::EdgeIterator eit(cur_he->getTwin()->getFace());
    while (eit.hasNext())
    {
		Arrangement::HalfEdge* next_he = eit.getNext();
		Arrangement::Face* adj_f = next_he->getTwin()->getFace();
		if (adj_f->getData().state == Arrangement::FaceData::NOTPASSED) {
			m_he_stack.push(next_he);
			adj_f->getData().state = Arrangement::FaceData::CHECKED;
		}
    }

	// ending condition of DFS
	if (m_he_stack.empty()) return DFS_END;

    // setting the next target
    cur_he = m_he_stack.top();
    cur_f = cur_he->getTwin()->getFace();

    // rewind path to the new target
    // backtrack until reaching to the parent of current target
    while (cur_he->getFace() != m_he_path.top()->getTwin()->getFace())
    {
        // backtracking with CS updates
		Arrangement::HalfEdge* back_he = m_he_path.top();
        update_CS(back_he->getTwin()); // backward

        m_he_path.pop();
    }
    // connect path to the target
    m_he_path.push(cur_he);

    // CS update
    update_CS(cur_he);

    return DFS_ADVANCED;
}

bool TranslationSpaceSubdivision::solveLP(BasisResult &b)
{
    bool all_matched = true;
    for (unsigned int i=0; i < t2->getMesh().getNumVertices(); ++i) {
        if (t2->getMesh().getVertex(i)->getData().VTpair == NULL) { // if unmatched one exists,
            all_matched = false; // patch is not in domain.
            break;
        }
    }
    // Execute LP-type problem solver only if the patch points are all on the valid domain.
    if (all_matched) {
		Arrangement::HalfEdge* he = m_he_stack.top();
		Arrangement::Face* f = he->getTwin()->getFace();
        std::vector<Point> boundary_pts;

		Arrangement::EdgeIterator eit(f);
		bool hasConnectingCorner = false;
        while (eit.hasNext())
        {
			Arrangement::HalfEdge* he = eit.getNext();
			Arrangement::VertexData& vd = he->getOrigin()->getData();
            Point p(vd.x, vd.y, 0);
            boundary_pts.push_back(p);
			if (arr.isNextFirstHalfEdge(he)) hasConnectingCorner = true;
        }
		nextGrid_CS = m_CS; // remember CS for the next grid.
        b = m_CS.solveLP(boundary_pts);
        return true;
    }
    else {
        //Translation is in the invalid domain.
        return false;
    }
}

void TranslationSpaceSubdivision::switch_VTpair(TerrainHalfEdge *he, TerrainVertex *v, bool he_is_from_patch)
{
    TerrainFace* f1 = he->getFace();
    TerrainFace* f2 = he->getTwin()->getFace();

    // previous vertex-triangle pair
    VertexTrianglePair *&ref_VTpair = v->getData().VTpair;

    // update vertex-triangle pair
    // remove previous pair
    if (ref_VTpair != NULL) {
        TerrainFace* fh = ref_VTpair->pTriangle;
        //bool fh_ck = ref_VTpair->tri_is_from_patch;

        m_CS.removePair(ref_VTpair);

        if (f1 == fh) {
            if (f2 == NULL)
                ref_VTpair = NULL;
            else
                ref_VTpair = m_CS.insertPair(v, f2, he_is_from_patch);
        }
        else if (f2 == fh) {
            if (f1 == NULL)
                ref_VTpair = NULL;
            else
                ref_VTpair = m_CS.insertPair(v, f1, he_is_from_patch);
        }
        else
            throw cpp::Exception("The triangle should contain e as boundary.");
    }
    else {
        // e should be boundary (one side is unbounded face)
        if (f1 != NULL && f2 != NULL) throw cpp::Exception("e should be boundary when p_VTpair==NULL.");

        if (f1 == NULL)
            ref_VTpair = m_CS.insertPair(v, f2, he_is_from_patch);
        else
            ref_VTpair = m_CS.insertPair(v, f1, he_is_from_patch);
    }
}

void TranslationSpaceSubdivision::switch_EEpairs(TerrainHalfEdge *he, TerrainVertex *v, bool he_is_from_patch)
{
	std::list<EdgeEdgePair *> &EEpair_list = he->getData().edgeData->EEpair_list;

    std::list<std::list<EdgeEdgePair *>::iterator> EEpair_it_v; // subset of EEpair_list related to v.

    // filter EEpair_list
    std::list<EdgeEdgePair *>::iterator eepit = EEpair_list.begin();
    while (eepit != EEpair_list.end())
    {
		if ((*eepit)->ed1->halfEdge_up->getOrigin() == v 
			|| (*eepit)->ed1->halfEdge_down->getOrigin() == v
			|| (*eepit)->ed2->halfEdge_up->getOrigin() == v
			|| (*eepit)->ed2->halfEdge_down->getOrigin() == v) 
		{
            EEpair_it_v.push_back(eepit);
        }
		++eepit;
    }

    // k^2 algorithm...
    Terrain::EdgeIterator eit_v(v);
    while (eit_v.hasNext())
    {
        TerrainHalfEdge *cur_he = eit_v.getNext();
		std::list<EdgeEdgePair *> &cur_EEpair_list = cur_he->getData().edgeData->EEpair_list;

        bool exist = false; // true if (cur_he, he) pair is already in EEpair_list.
        std::list<std::list<EdgeEdgePair *>::iterator>::iterator eepit_v = EEpair_it_v.begin();
        while (eepit_v != EEpair_it_v.end())
        {
			if ((**eepit_v)->ed1 == cur_he->getData().edgeData
				|| (**eepit_v)->ed2 == cur_he->getData().edgeData) 
			{
                exist = true;
                break;
            }
            else
                ++eepit_v;
        }

        if (exist) {
			EdgeEdgePair *toBeRemoved = **eepit_v;
            EEpair_list.erase(*eepit_v);
            EEpair_it_v.erase(eepit_v);
			std::list<EdgeEdgePair *>::iterator it_cur_EEpair_list = cur_EEpair_list.begin();
			while (it_cur_EEpair_list != cur_EEpair_list.end())
			{
				if ((*it_cur_EEpair_list)->ed1 == he->getData().edgeData
					|| (*it_cur_EEpair_list)->ed2 == he->getData().edgeData) 
				{
					cur_EEpair_list.erase(it_cur_EEpair_list);
					break;
				}
				else
					++it_cur_EEpair_list;
			}
			m_CS.removePair(toBeRemoved);
        }
        else {
            EdgeEdgePair *eepair = m_CS.insertPair(he->getData().edgeData, cur_he->getData().edgeData, he_is_from_patch);
            EEpair_list.push_back(eepair);
			cur_EEpair_list.push_back(eepair);
        }
    }
}

void TranslationSpaceSubdivision::update_CS(Arrangement::HalfEdge *he)
{
	Arrangement::EdgeData::SourceIterator sit;

    int multiple_edge_number = 0;

	// for all the data of hh, update CS
	std::vector<Arrangement::EdgeData::Source> &sources = he->getData().edgeData->sources;
	for (sit = sources.begin(); sit != sources.end(); ++sit) {
		// For vertex-triangle pairs...
		switch_VTpair(sit->he, sit->v, sit->he_is_from_patch);

		// For edge-edge pairs...
		switch_EEpairs(sit->he, sit->v, sit->he_is_from_patch);

		// count multiple edge
        ++multiple_edge_number;
    }

    if (multiple_edge_number > 2) {
        std::cout << "Warning: multi-edge has >2 edges";
    }
}