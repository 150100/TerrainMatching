#include "terrain.h"

#include <iostream>
#include <map>

Terrain::Terrain()
{
    /* initialize members */
	x_min = y_min = z_min = edgeLength_min = std::numeric_limits<double>::infinity();
	x_max = y_max = z_max = edgeLength_max = -std::numeric_limits<double>::infinity();
}

void Terrain::loadData(const char *filename)
{
    WavefrontObjImporter<WavefrontMesh> importer;
    WavefrontMesh wMesh;

    /* import the file */
    importer.import(filename, wMesh);

    /* copy only DCEL structure of wMash */
    DCELStream<WavefrontMesh>::copyDcelData(wMesh, mesh);

    /* initialize */
	x_min = y_min = z_min = std::numeric_limits<double>::infinity();
	x_max = y_max = z_max = -std::numeric_limits<double>::infinity();

    /* update vertexData */
	bool isolated_start = false;
	unsigned int isolated_start_idx = -1;
    for (unsigned int i=0; i < wMesh.getNumVertices(); ++i)
    {
        WavefrontMesh::Vertex* originVertex = wMesh.getVertex(i);		
        TerrainMesh::Vertex* targetVertex = mesh.getVertex(i);

		// check isolated vertex.
		if (originVertex->getIncidentEdge() == NULL) {
			if (isolated_start == false) {
				isolated_start_idx = i;
				isolated_start = true;
			}
		}
		else if (isolated_start) {
			throw cpp::Exception("isolated vertex should be only on the end of the structure.");
		}

		// copy vertex data.
        TerrainMesh::VertexData &targetVertexData = targetVertex->getData();
        targetVertexData.setCoordinates(originVertex->getData().position.x,
                                        originVertex->getData().position.y,
                                        originVertex->getData().position.z);

		// update range.
        if (targetVertexData.p.x < x_min) x_min = targetVertexData.p.x;
        if (targetVertexData.p.y < y_min) y_min = targetVertexData.p.y;
        if (targetVertexData.p.z < z_min) z_min = targetVertexData.p.z;
        if (targetVertexData.p.x > x_max) x_max = targetVertexData.p.x;
        if (targetVertexData.p.y > y_max) y_max = targetVertexData.p.y;
        if (targetVertexData.p.z > z_max) z_max = targetVertexData.p.z;
    }

	/* remove isolated vertices. */
	if (isolated_start) {
		auto &meshVertices = mesh.getVertices();
		meshVertices.erase(meshVertices.begin() + isolated_start_idx, meshVertices.end());
	}
	
	/* update edgeDataContainer */
	edgeDataContainer.reserve(number_of_halfedges() / 2);

	double edgeLength2_min(std::numeric_limits<double>::infinity()), edgeLength2_max(-std::numeric_limits<double>::infinity());
	for (unsigned int i=0; i < number_of_halfedges(); ++i)
	{
		TerrainMesh::HalfEdge *he = mesh.getHalfEdge(i);
		if (he->getData().edgeData == NULL) {
			// check zero-length edge.
			if (he->getOrigin()->getData().p.x == he->getTwin()->getOrigin()->getData().p.x &&
				he->getOrigin()->getData().p.y == he->getTwin()->getOrigin()->getData().p.y) {
				throw cpp::Exception("There should not be a zero-length edge.");
			}

			EdgeData ed;
			edgeDataContainer.push_back(ed);
			he->getData().edgeData = &edgeDataContainer.back();
			he->getTwin()->getData().edgeData = &edgeDataContainer.back();
			
			// update edgeLength range.
			Point &p1 = he->getOrigin()->getData().p;
			Point &p2 = he->getTwin()->getOrigin()->getData().p;
			double length2 = pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2);
			if (length2 < edgeLength2_min) edgeLength2_min = length2;
			if (length2 > edgeLength2_max) edgeLength2_max = length2;
		}
	}
	edgeLength_min = sqrt(edgeLength2_min);
	edgeLength_max = sqrt(edgeLength2_max);

	if (edgeDataContainer.size() != number_of_halfedges() / 2)
		throw cpp::Exception("# of Edges is not (# of Halfedges)/2.");
}

void Terrain::loadData(std::vector<double> coordinates, std::vector<unsigned int> triangles)
{
    std::vector<TerrainMesh::Vertex>& vertices = mesh.getVertices();
    std::vector<TerrainMesh::HalfEdge>& halfEdges = mesh.getHalfEdges();
    std::vector<TerrainMesh::Face>& faces = mesh.getFaces();

    /* coords and tris should be triples */
    if (coordinates.size() % 3 != 0)
        throw cpp::Exception("Coordinates are not triples.");
    if (triangles.size() % 3 != 0)
        throw cpp::Exception("Triangle indexes are not triples.");

    /* allocate memory for mesh */
    unsigned int vSize = coordinates.size() / 3;
    unsigned int tSize = triangles.size() / 3;
    unsigned int heSize = 2 * (vSize + tSize - 1); // F(excluding unbounded face) + V - E = 1
    vertices.reserve(vSize);
    halfEdges.reserve(heSize);
    faces.reserve(tSize);

	/* initialize */
	x_min = y_min = z_min = std::numeric_limits<double>::infinity();
	x_max = y_max = z_max = -std::numeric_limits<double>::infinity();

    /* create vertices */
    for (unsigned int i=0; i < coordinates.size(); i=i+3)
    {
        mesh.createGetVertex()->getData().setCoordinates(coordinates[i], coordinates[i+1], coordinates[i+2]);

		// update range.
		if (coordinates[i] < x_min) x_min = coordinates[i];
		if (coordinates[i+1] < y_min) y_min = coordinates[i+1];
		if (coordinates[i+2] < z_min) z_min = coordinates[i+2];
		if (coordinates[i] > x_max) x_max = coordinates[i];
		if (coordinates[i+1] > y_max) y_max = coordinates[i+1];
		if (coordinates[i+2] > z_max) z_max = coordinates[i+2];
    }

    /* create triangles */
    for (unsigned int i=0; i < triangles.size(); i=i+3)
    {
        mesh.createTriangularFace(triangles[i], triangles[i+1], triangles[i+2]);
    }

    /* check triangles validity */
    mesh.checkAllFaces();

	/* update edgeDataContainer */
	edgeDataContainer.reserve(number_of_halfedges() / 2);

	double edgeLength2_min(std::numeric_limits<double>::infinity()), edgeLength2_max(-std::numeric_limits<double>::infinity());
	for (unsigned int i = 0; i < number_of_halfedges(); ++i)
	{
		TerrainMesh::HalfEdge *he = mesh.getHalfEdge(i);
		if (he->getData().edgeData == NULL) {
			// check zero-length edge.
			if (he->getOrigin()->getData().p.x == he->getTwin()->getOrigin()->getData().p.x &&
				he->getOrigin()->getData().p.y == he->getTwin()->getOrigin()->getData().p.y) {
				throw cpp::Exception("There should not be a zero-length edge.");
			}

			EdgeData ed;
			edgeDataContainer.push_back(ed);
			he->getData().edgeData = &edgeDataContainer.back();
			he->getTwin()->getData().edgeData = &edgeDataContainer.back();

			// update edgeLength range.
			Point &p1 = he->getOrigin()->getData().p;
			Point &p2 = he->getTwin()->getOrigin()->getData().p;
			double length2 = pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2);
			if (length2 < edgeLength2_min) edgeLength2_min = length2;
			if (length2 > edgeLength2_max) edgeLength2_max = length2;
		}
	}
	edgeLength_min = sqrt(edgeLength2_min);
	edgeLength_max = sqrt(edgeLength2_max);

	if (edgeDataContainer.size() != number_of_halfedges() / 2)
		throw cpp::Exception("# of Edges is not (# of Halfedges)/2.");
}

//Arrangement Terrain::createGetXYArrangement()
//{
//}
void Terrain::createGetTrimTerrain(double x_min, double x_max, double y_min, double y_max, Terrain *trimTerrain)
{
	std::vector<TerrainMesh::Vertex> &vertices = mesh.getVertices();
	std::vector<TerrainMesh::HalfEdge> &halfEdges = mesh.getHalfEdges();
	std::vector<TerrainMesh::Face> &faces = mesh.getFaces();

	double x_min_coord = (this->x_max - this->x_min) * x_min + this->x_min;
	double x_max_coord = (this->x_max - this->x_min) * x_max + this->x_min;
	double y_min_coord = (this->y_max - this->y_min) * y_min + this->y_min;
	double y_max_coord = (this->y_max - this->y_min) * y_max + this->y_min;

	std::vector<double> coordinates;
	std::vector<unsigned int> triangles;

	std::map<unsigned int, unsigned int> map_inserted_vertices; // map from *this to t
	
	unsigned int t_idx = 0;
	for (unsigned int i = 0; i < vertices.size(); ++i) { // search all the vertices
		TerrainMesh::Vertex &v = vertices.at(i);
		Point &v_p = v.getData().p;
		if (x_min_coord <= v_p.x && v_p.x <= x_max_coord &&
			y_min_coord <= v_p.y && v_p.y <= y_max_coord) { // if a vertex is in the range, insert the vertex into the trim-terrain.
			coordinates.push_back(v_p.x);
			coordinates.push_back(v_p.y);
			coordinates.push_back(v_p.z);
			map_inserted_vertices[i] = t_idx++;
		}
	}

	for (unsigned int i = 0; i < faces.size(); ++i) { // search all the faces
		TerrainMesh::Face &f = faces.at(i);
		unsigned int v1, v2, v3;
		TerrainMesh::HalfEdge *he = f.getBoundary();
#ifdef _DEBUG
		TerrainMesh::HalfEdge *he_first = he;
#endif
		v1 = mesh.getVertexId(he->getOrigin());
		he = he->getNext();
		v2 = mesh.getVertexId(he->getOrigin());
		he = he->getNext();
		v3 = mesh.getVertexId(he->getOrigin());
#ifdef _DEBUG
		if (he->getNext() != he_first)
			throw cpp::Exception("f is not a triangle.");
#endif
		std::map<unsigned int, unsigned int>::iterator it_v1 = map_inserted_vertices.find(v1);
		std::map<unsigned int, unsigned int>::iterator it_v2 = map_inserted_vertices.find(v2);
		std::map<unsigned int, unsigned int>::iterator it_v3 = map_inserted_vertices.find(v3);
		if (it_v1 != map_inserted_vertices.end() && it_v2 != map_inserted_vertices.end() && it_v3 != map_inserted_vertices.end()) { // if all the vertices of the triangle is inserted, insert that triangle.
			triangles.push_back(it_v1->second);
			triangles.push_back(it_v2->second);
			triangles.push_back(it_v3->second);
		}
	}

	trimTerrain->loadData(coordinates, triangles);
}
