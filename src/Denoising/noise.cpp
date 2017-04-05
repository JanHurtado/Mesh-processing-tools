#include "noise.h"

const num_t default_noise_level = 0.3;

num_t generateRandomGaussian(num_t mean, num_t StandardDeviation)
{
    static num_t v1, v2, s;
    static int phase = 0;
    num_t x;

    if(phase == 0)
    {
        do
        {
            v1 = -1 + 2 * (num_t)rand() / (num_t) RAND_MAX;
            v2 = -1 + 2 * (num_t)rand() / (num_t) RAND_MAX;
            s = v1 * v1 + v2 * v2;
        }while (s >= 1 || s == 0);

        x = v1 * sqrt(-2 * log(s) / s);
    }
    else
        x = v2 * sqrt(-2 * log(s) / s);

    phase = 1 - phase;

    return x * StandardDeviation + mean;
}

void randomGaussianNumbers(num_t mean, num_t StandardDeviation, int number, vector<num_t> &RandomNumbers)
{
    RandomNumbers.resize(number, 0.0);

    srand((unsigned int)time(NULL));
    for(int i = 0; i < number; i++){
        RandomNumbers[i] = generateRandomGaussian(mean, StandardDeviation);
    }
}

TriMesh addNoise(TriMesh & _mesh)
{
    TriMesh mesh = _mesh;
    // compute average length of mesh
    num_t average_length = 0.0;
    for(TriMesh::EdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end(); e_it++)
        average_length += mesh.calc_edge_length(*e_it);
    num_t edge_numbers = (num_t)mesh.n_edges();
    average_length /= edge_numbers;
    num_t noise_level = default_noise_level;
    num_t standard_deviation = average_length * noise_level;

    mesh.request_face_normals();
    mesh.request_vertex_normals();
    mesh.update_normals();

    vector<num_t> GaussianNumbers;

    randomGaussianNumbers(0, standard_deviation, (int)mesh.n_vertices(), GaussianNumbers);

    for(TriMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++){

        TriMesh::Point p = mesh.point(*v_it) + mesh.normal(*v_it) * GaussianNumbers[v_it->idx()];
        mesh.set_point(*v_it, p);
    }
    return mesh;
}

TriMesh::Normal generateRandomDirection()
{
    double x, y, z, length;
    do{
        x = -1 + 2 * (double)rand() / (double) RAND_MAX;
        y = -1 + 2 * (double)rand() / (double) RAND_MAX;
        length = x * x + y * y;
    }while(length>1);

    const double r = 2 * std::sqrt(1 - length);

    x *= r;
    y *= r;
    z = 1 - 2 * length;

    return TriMesh::Normal(x, y, z);
}

void randomDirections(int number, std::vector<TriMesh::Normal> &RandomDirections)
{
    RandomDirections.resize(number, TriMesh::Normal(0.0, 0.0, 0.0));

    srand((unsigned int)time(NULL));
    for(int i = 0; i < number; i++){
        RandomDirections[i] = generateRandomDirection();
    }
}

TriMesh addNoiseR(TriMesh & _mesh)
{
    TriMesh mesh = _mesh;
    // compute average length of mesh
    num_t average_length = 0.0;
    for(TriMesh::EdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end(); e_it++)
        average_length += mesh.calc_edge_length(*e_it);
    num_t edge_numbers = (num_t)mesh.n_edges();
    average_length /= edge_numbers;
    num_t noise_level = default_noise_level;
    num_t standard_deviation = average_length * noise_level;

    mesh.request_face_normals();
    mesh.request_vertex_normals();
    mesh.update_normals();

    vector<num_t> GaussianNumbers;
    vector<TriMesh::Normal> RandomDirections;

    randomGaussianNumbers(0, standard_deviation, (int)mesh.n_vertices(), GaussianNumbers);
    randomDirections((int)mesh.n_vertices(),RandomDirections);

    for(TriMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++){

        TriMesh::Point p = mesh.point(*v_it) + RandomDirections[v_it->idx()] * GaussianNumbers[v_it->idx()];
        mesh.set_point(*v_it, p);
    }
    return mesh;
}
