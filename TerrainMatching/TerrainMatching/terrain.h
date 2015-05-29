#ifndef TERRAIN_H
#define TERRAIN_H

#include <cstddef>
#include <fstream>
#include <iostream>
#include <list>

#include "DCEL/DCELStream.h"
#include "DCEL/Mesh.h"
#include "DCEL/WavefrontObjImporter.h"
#include "DCEL/Vector3.h"
#include <fstream>
#include <limits>

#include "plane.h"

///

// auxiliary definition of types of combinatorial pairs
class VertexTrianglePair;
class EdgeEdgePair;

// auxiliary definition of arrangement datastructure
class Arrangement;

// auxiliary definition of edge datastructure
class TerrainEdgeData;

///

// information of vertex
class TerrainVertexData
{
public:
    Point p; // coordinates
    VertexTrianglePair *VTpair; // paired one while translating.
	unsigned int arrIndex; // index for arrangement structure.

    TerrainVertexData()
        : VTpair(NULL), arrIndex(-1) {}

    void setCoordinates(double _x, double _y, double _z) {p.x = _x; p.y = _y; p.z = _z;}
};

// auxiliary datastructure for WavefrontObjImporter
class VertexData3f
{
public:
    Vector3f position;
};

// information of halfedge.
class TerrainHalfEdgeData
{
public:
    TerrainEdgeData* edgeData; // link to edge data.
	unsigned int arrIndex; // index for arrangement structure.

	TerrainHalfEdgeData() : edgeData(NULL), arrIndex(-1) {}
};

// information of face
class TerrainFaceData
{
public:
    // The paired one while translating will be stored on a vertex.

    TerrainFaceData() {}
};

// information of edge
class TerrainEdgeData
{
public:
    bool check_removed;
    std::list<EdgeEdgePair *> EEpair_list; // paired one while translating

    TerrainEdgeData() : check_removed(false) {}
};


typedef Mesh<VertexData3f, TerrainHalfEdgeData, TerrainFaceData> WavefrontMesh;

///

class Terrain
{
public:
    typedef TerrainVertexData VertexData;
    typedef TerrainHalfEdgeData HalfEdgeData;
    typedef TerrainFaceData FaceData;

    typedef Mesh<TerrainVertexData, TerrainHalfEdgeData, TerrainFaceData> TerrainMesh;

	typedef TerrainMesh::Vertex Vertex;
	typedef TerrainMesh::HalfEdge HalfEdge;
	typedef TerrainMesh::Face Face;
    typedef TerrainMesh::EdgeIterator EdgeIterator;
	
    typedef TerrainEdgeData EdgeData;

    /**
     * @brief Terrain Constructor.
     */
    Terrain();

    /**
     * @brief loadData Load a terrain data from Wavefront OBJ file.
     * @param filename OBJ file name.
     */
    void loadData(const char* filename);

    /**
     * @brief loadData Load a terrain data from the built-in code.
     * @param coordinates Vector of coordinates. The size should be multiple of 3.
     * @param triangles Vector of triangle indexes. The size should be multiple of 3.
     */
    void loadData(std::vector<double> coordinates, std::vector<unsigned int> triangles);

    /**
     * @brief createGetTrimTerrain Make a trimmed terrain and return it.
     * @param x_min Minimum x-axis-ratio of the trim area
     * @param x_max Maximum x-axis-ratio of the trim area
     * @param y_min Minimum y-axis-ratio of the trim area
     * @param y_max Maximum y-axis-ratio of the trim area
     * @return The terrain trimmed by the specified rectangular region.
     */
	void createGetTrimTerrain(double x_min, double x_max, double y_min, double y_max, Terrain *trimTerrain);

	// Sort the vertices and map data of edges to the corresponding vertices.
	void sortVerticesAndUpdate();

    unsigned int number_of_vertices() {return mesh.getNumVertices();}
    unsigned int number_of_halfedges() {return mesh.getNumHalfEdges();}
	unsigned int number_of_edges() {return edgeDataContainer.size();}
    unsigned int number_of_faces() {return mesh.getNumFaces();}

    TerrainMesh& getMesh() { return mesh; }
	std::vector<EdgeData>& getEdgeDataContainer() { return edgeDataContainer; }

	void print(std::ostream &os) {
		std::cerr << "V = " << number_of_vertices() << std::endl;
		std::cerr << "E = " << number_of_edges() << std::endl;
		std::cerr << "F = " << number_of_faces() << std::endl;
		std::cerr << "range = ([" << x_min << ',' << x_max << "], [" << y_min << ',' << y_max << "], [" << z_min << ',' << z_max << "])\n";
		std::cerr << "edge lengths = [" << edgeLength_min << ',' << edgeLength_max << "]\n";
	}

	double x_min, x_max, y_min, y_max, z_min, z_max;
	double edgeLength_min, edgeLength_max;

protected:
    TerrainMesh mesh;
	std::vector<EdgeData> edgeDataContainer;
};

class TerrainWithGrids : public Terrain
{
public:
	bool gridIsMade;

	TerrainWithGrids() : Terrain(), gridIsMade(false) {}

	// Make grids and distribute vertices into them. Make sure that the mesh data is loaded.
	void makeGrids(double _gridStepSize);

	// Get vertices in the range.
	void appendVerticesInRange(double rangeX_min, double rangeX_max, double rangeY_min, double rangeY_max, std::vector<Vertex *> *verticesInRange);
	
protected:
	unsigned gridSizeX, gridSizeY;
	double gridStepSize;
	std::vector<Vertex *> **grids;
};

#endif // TERRAIN_H
