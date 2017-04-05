#ifndef UTIL_H
#define UTIL_H

#include "mesh.h"
#include <vector>
#include <unordered_map>

using namespace std;

#define PI 3.14159265359

num_t getArea (TriMesh & mesh);
num_t getVolume ( TriMesh & mesh );

num_t getAverageEdgeLength(TriMesh &mesh);
void getAllFaceAreas(TriMesh &mesh, vector<num_t> &areas);
void getAllFaceCentroids(TriMesh &mesh, vector<TriMesh::Point> &centroids);
void getAllFaceNormals(TriMesh &mesh, vector<TriMesh::Normal> &normals);
void getFaceVertexAngle(TriMesh &mesh,TriMesh::FaceHandle fh,TriMesh::VertexHandle vh, num_t & angle);
num_t getVertexArea(TriMesh & mesh,TriMesh::VertexHandle vh,vector<num_t> & areas);
void getAllVertexAreas(TriMesh & mesh, vector<num_t> & areas);
void getAllPoints(TriMesh &mesh, vector<TriMesh::Point> & points);
void getAllConnectedComponentVertices(TriMesh & mesh, vector<vector<int> > & connectedComponentVertices);
void getAllConnectedComponentAreas(TriMesh & mesh, vector<vector<int> > & connectedComponentVertices, vector<num_t> & connectedComponentAreas);
void removeSmallConnectedComponents(TriMesh & mesh, num_t threshold);

inline num_t GaussianWeight(num_t distance, num_t sigma)
{
    return exp( -0.5 * distance * distance / (sigma * sigma));
}

inline num_t NormalDistance(const TriMesh::Normal &n1, const TriMesh::Normal &n2)
{
    return (n1 - n2).length();
}

void uniformFaceSampling(TriMesh & mesh,vector<int> & sampled_faces);

#endif // UTIL_H
