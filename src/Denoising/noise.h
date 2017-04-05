#ifndef NOISE_H
#define NOISE_H

#include "mesh.h"

num_t generateRandomGaussian(num_t mean, num_t StandardDeviation);

void randomGaussianNumbers(num_t mean, num_t StandardDeviation, int number, vector<num_t> &RandomNumbers);

TriMesh addNoise(TriMesh & _mesh);

TriMesh::Normal generateRandomDirection();

void randomDirections(int number, vector<TriMesh::Normal> &RandomDirections);

TriMesh addNoiseR(TriMesh & _mesh);

#endif // NOISE_H
