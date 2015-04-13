#include "terrain.h"

#include <iostream>

Terrain::Terrain()
{
    /* initialize members */
    x_min = y_min = z_min = std::numeric_limits<double>::infinity();
    x_max = y_max = z_max = -std::numeric_limits<double>::infinity();
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

    /* update min-max (data part) */
    for (unsigned int i=0; i < wMesh.getNumVertices(); ++i)
    {
        WavefrontMesh::Vertex* originVertex = wMesh.getVertex(i);
        TerrainMesh::Vertex* targetVertex = mesh.getVertex(i);
        TerrainMesh::VertexData &targetVertexData = targetVertex->getData();
        targetVertexData.setCoordinates(originVertex->getData().position.x,
                                        originVertex->getData().position.y,
                                        originVertex->getData().position.z);

        if (targetVertexData.p.x < x_min) x_min = targetVertexData.p.x;
        if (targetVertexData.p.y < y_min) y_min = targetVertexData.p.y;
        if (targetVertexData.p.z < z_min) z_min = targetVertexData.p.z;

        if (targetVertexData.p.x > x_max) x_max = targetVertexData.p.x;
        if (targetVertexData.p.y > y_max) y_max = targetVertexData.p.y;
        if (targetVertexData.p.z > z_max) z_max = targetVertexData.p.z;
    }
	
	/* update edgeDataContainer */
	edgeDataContainer.reserve(number_of_halfedges() / 2);

	for (unsigned int i=0; i < number_of_halfedges(); ++i)
	{
		if (mesh.getHalfEdge(i)->getData().edgeData == NULL) {
			EdgeData ed;
			edgeDataContainer.push_back(ed);
			mesh.getHalfEdge(i)->getData().edgeData = &edgeDataContainer.back();
			mesh.getHalfEdge(i)->getTwin()->getData().edgeData = &edgeDataContainer.back();
		}
	}

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

    /* create vertices */
    for (unsigned int i=0; i < coordinates.size(); i=i+3)
    {
        mesh.createGetVertex()->getData().setCoordinates(coordinates[i], coordinates[i+1], coordinates[i+2]);
    }

    /* create triangles */
    for (unsigned int i=0; i < triangles.size(); i=i+3)
    {
        mesh.createTriangularFace(triangles[i], triangles[i+1], triangles[i+2]);
    }

    /* check triangles validity */
    mesh.checkAllFaces();
}

//Arrangement Terrain::createGetXYArrangement()
//{
//}

//Terrain Terrain::createGetTrimTerrain(double x_min, double x_max, double y_min, double y_max)
//{
//}
