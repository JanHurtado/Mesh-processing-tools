#ifndef CURVATURE_H
#define CURVATURE_H

//#include "mesh.h"
#include "util.h"

void getVertexGaussianCurvature(TriMesh &mesh, vector<num_t> & curvature);

void getVertexMeanCurvature(TriMesh & mesh,vector<num_t> & curvature);

void getVertexPrincipalCurvaturesSum(TriMesh & mesh,vector<num_t> & mean_curvature,vector<num_t> & gaussian_curvature,vector<num_t> & principal_curvature);

void getAllCurvatures(TriMesh & mesh,vector<num_t> & mean_curvature,vector<num_t> & gaussian_curvature,vector<num_t> & principal_curvature);


#endif // CURVATURE_H
