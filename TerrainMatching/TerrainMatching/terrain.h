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
	bool isolated;

    TerrainVertexData()
        : VTpair(NULL), isolated(false) {}

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
    TerrainEdgeData* edgeData; // Link to edge data.

	TerrainHalfEdgeData() : edgeData(NULL) {}
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

//    /**
//     * @brief createGetXYArrangement Make XY-plane projection and return it.
//     * @return XY-plane projected arrangement of the terrain.
//     */
//    Arrangement createGetXYArrangement();

    /**
     * @brief createGetTrimTerrain Make a trimmed terrain and return it.
     * @param x_min Minimum x-axis-ratio of the trim area
     * @param x_max Maximum x-axis-ratio of the trim area
     * @param y_min Minimum y-axis-ratio of the trim area
     * @param y_max Maximum y-axis-ratio of the trim area
     * @return The terrain trimmed by the specified rectangular region.
     */
    Terrain createGetTrimTerrain(double x_min, double x_max, double y_min, double y_max);

    unsigned int number_of_vertices() {return mesh.getNumVertices();}
    unsigned int number_of_halfedges() {return mesh.getNumHalfEdges();}
	unsigned int number_of_edges() {return edgeDataContainer.size();}
    unsigned int number_of_faces() {return mesh.getNumFaces();}

    TerrainMesh& getMesh() {return mesh;}
	std::vector<EdgeData>& getEdgeDataContainer() { return edgeDataContainer; }

private:
    double x_min, x_max, y_min, y_max, z_min, z_max;
    TerrainMesh mesh;
	std::vector<EdgeData> edgeDataContainer;
};


#endif // TERRAIN_H
