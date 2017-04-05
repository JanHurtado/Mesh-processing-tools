#ifndef DENOISING_H
#define DENOISING_H

#include "neighborhood.h"
#include "util.h"
#include "custom.h"

/**
 * @brief The FaceNeighborType enum
 * kVertexBased : vertex based face neighborhood
 * kEdgeBased : edge based face neighborhood
 * kRadiusBased : radius based face neighborhood
 */
enum FaceNeighborType{kVertexBased, kEdgeBased, kRadiusBased};

/** @defgroup groupAuxiliarFunctions Auxiliar functions group
 *  This group contains auxiliar functions for mesh denoising algorithms
 *  @{
 */

/**
 * @brief updateVertexPositions: this function compute new vertex positions adapted to an input normal field.
 * @param mesh: input mesh.
 * @param filtered_normals: input normal field.
 * @param iteration_number: number of iterations of the optimization algorithm (gradient descent approach).
 * @param fixed_boundary: determines if a boundary vertex will be updated
 */
void updateVertexPositions(TriMesh &mesh, vector<TriMesh::Normal> &filtered_normals, int iteration_number, bool fixed_boundary);

/**
 * @brief getSigmaC: sigma_c computation based on the average of face centroid distances (only between adjacent faces), and multiplied by a scalar
 * @param mesh: input mesh
 * @param face_centroids: precomputed face centroids
 * @param sigma_c_scalar: scalar
 * @return value of sigma_c
 */
num_t getSigmaC(TriMesh &mesh, vector<TriMesh::Point> &face_centroids, num_t sigma_c_scalar);

/**
 * @brief getRadius: computation of radius for radius-based face neighborhood, based on the average of face centroid distances (only between adjacent faces) and multiplied by a scalar
 * @param mesh: input mesh
 * @param scalar
 * @return radius value
 */
num_t getRadius(TriMesh &mesh, num_t scalar);

/** @} */ // end of group1



/** @defgroup isotropicMeshSmoothingAlgorithms Isotropic mesh denoising algorithms
 *  This group contains isotropic mesh denoising algorithms
 *  @{
 */

/////////////////////////////////////////////////
/// Uniform Laplacian smoothing
/////////////////////////////////////////////////

TriMesh uniformLaplacian(TriMesh & _mesh,int iteration_number,num_t scale);

/////////////////////////////////////////////////
/// Taubin smoothing (G. Taubin)
/////////////////////////////////////////////////

TriMesh taubin(TriMesh & _mesh,int iteration_number,num_t lambda,num_t mu);

/////////////////////////////////////////////////
/// HC Laplacian smoothing (Vollmer et al.)
/////////////////////////////////////////////////

TriMesh HCLaplacian(TriMesh & _mesh,int iteration_number,num_t alpha,num_t beta);

/** @} */


/** @defgroup anisotropicMeshSmoothingAlgorithms Anisotropic mesh denoising algorithms
 *  This group contains anisotropic mesh denoising algorithms
 *  @{
 */

/////////////////////////////////////////////////
/// Hildebrandt and Polthier smoohing
/////////////////////////////////////////////////


num_t HildebrandtAndPolthier_weight_function(num_t lambda,num_t transition_width,num_t value);

TriMesh HildebrandtAndPolthier(TriMesh & _mesh,int iteration_number,num_t lambda,num_t transition_width,num_t step_width);

/////////////////////////////////////////////////
/// Bilateral mesh denoising (Fleishman et al.)
/////////////////////////////////////////////////

TriMesh bilateral(TriMesh & _mesh,int iteration_number);

/////////////////////////////////////////////////
/// Fast and effective feature-presearving mesh denoising (Sun et al.)
/////////////////////////////////////////////////

TriMesh FastAndEffectiveFeaturePreserving(TriMesh & _mesh, int normal_iteration_number, int vertex_iteration_number, num_t threshold_T, FaceNeighborType fnt);

void updateFilteredNormalsFast(TriMesh &mesh, int normal_iteration_number, num_t threshold_T,FaceNeighborType fnt,vector<TriMesh::Normal> &filtered_normals);


/////////////////////////////////////////////////
/// Bilateral normal filtering for mesh denoising (Zheng et al.)
/////////////////////////////////////////////////

void updateFilteredNormals(TriMesh &mesh, int normal_iteration_number, num_t sigma_c_scalar, num_t sigma_s, vector<TriMesh::Normal> &filtered_normals);

TriMesh bilateralNormal(TriMesh &_mesh, int normal_iteration_number, int vertex_iteration_number, num_t sigma_c_scalar,  num_t sigma_s );


/////////////////////////////////////////////////
/// Guided mesh normal filtering (Zhang et al.)
/////////////////////////////////////////////////

void getVertexBasedFaceNeighbor(TriMesh &mesh, TriMesh::FaceHandle fh, vector<TriMesh::FaceHandle> &face_neighbor);

void getAllFaceNeighborGMNF(TriMesh &mesh, FaceNeighborType face_neighbor_type, num_t radius, bool include_central_face,
                                                       vector<vector<TriMesh::FaceHandle> > &all_face_neighbor);

void getAllGuidedNeighborGMNF(TriMesh &mesh, vector<vector<TriMesh::FaceHandle> > &all_guided_neighbor);



void getFaceNeighborInnerEdge(TriMesh &mesh, vector<TriMesh::FaceHandle> &face_neighbor, vector<TriMesh::EdgeHandle> &inner_edge);

void getRangeAndMeanNormal(TriMesh &mesh, vector<vector<TriMesh::FaceHandle> > &all_guided_neighbor,
                                                       vector<num_t> &face_areas, vector<TriMesh::Normal> &normals,
                                                       vector<pair<num_t, TriMesh::Normal> > &range_and_mean_normal);

void getGuidedNormals(TriMesh &mesh, vector<vector<TriMesh::FaceHandle> > &all_guided_neighbor,
                                                 vector<num_t> &face_areas, vector<TriMesh::Normal> &normals,
                                                 vector<pair<num_t, TriMesh::Normal> > range_and_mean_normal,
                                                 vector<TriMesh::Normal> &guided_normals);


void updateFilteredNormalsGuided(TriMesh &mesh, vector<TriMesh::Normal> &filtered_normals,num_t radius_scalar,
                                 num_t sigma_c_scalar, int normal_iteration_number, num_t sigma_s, int vertex_iteration_number);

TriMesh guided(TriMesh _mesh,int normal_iteration_number, int vertex_iteration_number, num_t sigma_c_scalar, num_t sigma_s, num_t radius_scalar);


/** @} */


class DenoisingAlgorithm
{
protected:
    TriMesh inputMesh;
    TriMesh outputMesh;
    bool executed;
public:
    DenoisingAlgorithm()
    {
        executed = 0;
    }
    DenoisingAlgorithm(TriMesh & _inputMesh)
    {
        inputMesh = _inputMesh;
        executed = 0;
    }
    virtual ~DenoisingAlgorithm(){}
    virtual void run()=0;
    virtual void setParameters(vector<num_t> & parameters)=0;
    virtual void setDefaultParameters()=0;

    void setInput(TriMesh & _inputMesh){inputMesh = _inputMesh; executed=0;}
    TriMesh getInput(){return inputMesh;}
    TriMesh getOutput()
    {
        if (executed)
            return outputMesh;
        else
        {
            run();
            return outputMesh;
        }
    }
};


class AlgorithmUniformLaplacian : public DenoisingAlgorithm
{
    int num_iterations;
public:
    AlgorithmUniformLaplacian():DenoisingAlgorithm(){}
    AlgorithmUniformLaplacian(TriMesh & _inputMesh):DenoisingAlgorithm(_inputMesh){}
    void run()
    {
        printf("Running Uniform Laplacian ... %i \n",num_iterations);
        printf("Parameter order: %s \n",uniform_parameter_order.c_str());
        outputMesh = uniformLaplacian(inputMesh,num_iterations,0.5);
        executed = 1;
    }
    void setParameters(vector<num_t> & _parameters)
    {
        vector<num_t> parameters = uniform_parameters_default;
        for(size_t i = 0; i<_parameters.size() && i<parameters.size() ; i++ )
            parameters[i] = _parameters[i];
        num_iterations = (int)parameters[0];
    }
    void setDefaultParameters()
    {
        num_iterations = (int)uniform_num_iterations_default;
    }
};

class AlgorithmTaubin : public DenoisingAlgorithm
{
    int num_iterations;
    num_t lambda;
    num_t mu;
public:
    AlgorithmTaubin():DenoisingAlgorithm(){}
    AlgorithmTaubin(TriMesh & _inputMesh):DenoisingAlgorithm(_inputMesh){}
    void run()
    {
        printf("Running Taubin ... %i %f %f \n",num_iterations,lambda,mu);
        printf("Parameter order: %s \n",taubin_parameter_order.c_str());
        outputMesh = taubin(inputMesh,num_iterations,lambda,mu);
        executed = 1;
    }
    void setParameters(vector<num_t> & _parameters)
    {
        vector<num_t> parameters = taubin_parameters_default;
        for(size_t i = 0; i<_parameters.size() && i<parameters.size() ; i++ )
            parameters[i] = _parameters[i];
        num_iterations = (int)parameters[0];
        lambda = parameters[1];
        mu = parameters[2];
    }
    void setDefaultParameters()
    {
        num_iterations = (int)taubin_num_iterations_default;
        lambda = taubin_lambda_default;
        mu = taubin_mu_default;
    }
};

class AlgorithmHCLaplacian : public DenoisingAlgorithm
{
    int num_iterations;
    num_t alpha;
    num_t beta;
public:
    AlgorithmHCLaplacian():DenoisingAlgorithm(){}
    AlgorithmHCLaplacian(TriMesh & _inputMesh):DenoisingAlgorithm(_inputMesh){}
    void run()
    {
        printf("Running HC Laplacian ... %i %f %f \n",num_iterations,alpha,beta);
        printf("Parameter order: %s \n",hc_parameter_order.c_str());
        outputMesh = HCLaplacian(inputMesh,num_iterations,alpha,beta);
        executed = 1;
    }
    void setParameters(vector<num_t> & _parameters)
    {
        vector<num_t> parameters = hc_parameters_default;
        for(size_t i = 0; i<_parameters.size() && i<parameters.size() ; i++ )
            parameters[i] = _parameters[i];
        num_iterations = (int)parameters[0];
        alpha = parameters[1];
        beta = parameters[2];
    }
    void setDefaultParameters()
    {
        num_iterations = (int)hc_num_iterations_default;
        alpha = hc_alpha_default;
        beta = hc_beta_default;
    }
};


class AlgorithmBilateral : public DenoisingAlgorithm
{
    int num_iterations;
public:
    AlgorithmBilateral():DenoisingAlgorithm(){}
    AlgorithmBilateral(TriMesh & _inputMesh):DenoisingAlgorithm(_inputMesh){}
    void run()
    {
        printf("Running Bilateral mesh denoising ... %i \n",num_iterations);
        printf("Parameter order: %s \n",bilateral_parameter_order.c_str());
        outputMesh = bilateral(inputMesh,num_iterations);
        executed = 1;
    }
    void setParameters(vector<num_t> & _parameters)
    {
        vector<num_t> parameters = bilateral_parameters_default;
        for(size_t i = 0; i<_parameters.size() && i<parameters.size() ; i++ )
            parameters[i] = _parameters[i];
        num_iterations = (int)parameters[0];
    }
    void setDefaultParameters()
    {
        num_iterations = (int)bilateral_num_iterations_default;
    }
};

class AlgorithmFastAndEffective : public DenoisingAlgorithm
{
    int num_normal_iterations;
    int num_vertex_iterations;
    num_t threshold_T;
public:
    AlgorithmFastAndEffective():DenoisingAlgorithm(){}
    AlgorithmFastAndEffective(TriMesh & _inputMesh):DenoisingAlgorithm(_inputMesh){}
    void run()
    {
        printf("Running Fast and effective feature-preserving mesh denoising ... %i %i %f \n",num_normal_iterations,num_vertex_iterations,threshold_T);
        printf("Parameter order: %s \n",fast_parameter_order.c_str());
        outputMesh = FastAndEffectiveFeaturePreserving(inputMesh, num_normal_iterations,num_vertex_iterations,threshold_T,kVertexBased);
        executed = 1;
    }
    void setParameters(vector<num_t> & _parameters)
    {
        vector<num_t> parameters = fast_parameters_default;
        for(size_t i = 0; i<_parameters.size() && i<parameters.size() ; i++ )
            parameters[i] = _parameters[i];
        num_normal_iterations = (int)parameters[0];
        num_vertex_iterations = (int)parameters[1];
        threshold_T = parameters[2];
    }
    void setDefaultParameters()
    {
        num_normal_iterations = (int)fast_num_normal_iterations_default;
        num_vertex_iterations = (int)fast_num_vertex_iterations_default;
        threshold_T = fast_threshold_T_default;
    }
};

class AlgorithmBilateralNormal : public DenoisingAlgorithm
{
    int num_normal_iterations;
    int num_vertex_iterations;
    num_t sigma_c_scalar;
    num_t sigma_s;
public:
    AlgorithmBilateralNormal():DenoisingAlgorithm(){}
    AlgorithmBilateralNormal(TriMesh & _inputMesh):DenoisingAlgorithm(_inputMesh){}
    void run()
    {
        printf("Running Bilateral Normal ... %i %i %f %f \n",num_normal_iterations,num_vertex_iterations,sigma_c_scalar,sigma_s);
        printf("Parameter order: %s \n",bilateralnormal_parameter_order.c_str());
        outputMesh = bilateralNormal(inputMesh,num_normal_iterations, num_vertex_iterations,sigma_c_scalar,sigma_s);
        executed = 1;
    }
    void setParameters(vector<num_t> & _parameters)
    {
        vector<num_t> parameters = bilateralnormal_parameters_default;
        for(size_t i = 0; i<_parameters.size() && i<parameters.size() ; i++ )
            parameters[i] = _parameters[i];
        num_normal_iterations = (int)parameters[0];
        num_vertex_iterations = (int)parameters[1];
        sigma_c_scalar = parameters[2];
        sigma_s = parameters[3];
    }
    void setDefaultParameters()
    {
        num_normal_iterations = (int)bilateralnormal_num_normal_iterations_default;
        num_vertex_iterations = (int)bilateralnormal_num_vertex_iterations_default;
        sigma_c_scalar = bilateralnormal_sigma_c_scalar_default;
        sigma_s = bilateralnormal_sigma_s_default;
    }
};

class AlgorithmGuided : public DenoisingAlgorithm
{
    int num_normal_iterations;
    int num_vertex_iterations;
    num_t sigma_c_scalar;
    num_t sigma_s;
    num_t radius_scalar;
public:
    AlgorithmGuided():DenoisingAlgorithm(){}
    AlgorithmGuided(TriMesh & _inputMesh):DenoisingAlgorithm(_inputMesh){}
    void run()
    {
        printf("Running Guided Mesh Normal Filtering ... %i %i %f %f %f \n",num_normal_iterations,num_vertex_iterations,sigma_c_scalar,sigma_s,radius_scalar);
        printf("Parameter order: %s \n",guided_parameter_order.c_str());
        outputMesh = guided(inputMesh,num_normal_iterations, num_vertex_iterations,sigma_c_scalar,sigma_s,radius_scalar);
        executed = 1;
    }
    void setParameters(vector<num_t> & _parameters)
    {
        vector<num_t> parameters = guided_parameters_default;
        for(size_t i = 0; i<_parameters.size() && i<parameters.size() ; i++ )
            parameters[i] = _parameters[i];
        num_normal_iterations = (int)parameters[0];
        num_vertex_iterations = (int)parameters[1];
        sigma_c_scalar = parameters[2];
        sigma_s = parameters[3];
        radius_scalar = parameters[4];
    }
    void setDefaultParameters()
    {
        num_normal_iterations = (int)guided_num_normal_iterations_default;
        num_vertex_iterations = (int)guided_num_vertex_iterations_default;
        sigma_c_scalar = guided_sigma_c_scalar_default;
        sigma_s = guided_sigma_s_default;
        radius_scalar = guided_radius_scalar_default;
    }
};

#endif // DENOISING_H
