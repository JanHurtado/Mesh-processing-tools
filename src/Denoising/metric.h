#ifndef METRIC_H
#define METRIC_H

#include "util.h"
#include "nanoflann.hpp"
#include "curvature.h"



using namespace nanoflann;


/** @defgroup groupKD-Tree KD-Tree group
 *  This group contains structures, definitions and functions for KD-Tree indexation
 *  @{
 */

/**
 * @brief The PointCloud struct: This struct is used by nanoflann library to represent point clouds indexed in a KD-Tree. See nanoflann documentation for more detail.
 */
struct PointCloud
{
    struct Point
    {
        num_t  x,y,z;
    };

    vector<Point>  pts;

    inline size_t kdtree_get_point_count() const { return pts.size(); }

    inline float kdtree_distance(const float *p1, const size_t idx_p2,size_t size) const
    {
        float d0=p1[0]-pts[idx_p2].x;
        float d1=p1[1]-pts[idx_p2].y;
        float d2=p1[2]-pts[idx_p2].z;
        return d0*d0+d1*d1+d2*d2;
    }

    inline num_t kdtree_distance(const num_t *p1, const size_t idx_p2,size_t size) const
    {
        num_t d0=p1[0]-pts[idx_p2].x;
        num_t d1=p1[1]-pts[idx_p2].y;
        num_t d2=p1[2]-pts[idx_p2].z;
        return d0*d0+d1*d1+d2*d2;
    }

    inline float kdtree_get_pt(const size_t idx, int dim) const
    {
        if (dim==0) return pts[idx].x;
        else if (dim==1) return pts[idx].y;
        else return pts[idx].z;
    }

    template <class BBOX>
    bool kdtree_get_bbox(BBOX &bb) const { return false; }

};

/**
 * @brief generateVertexPointCloud: Point cloud generation using vertex coordinates as points
 * @param mesh: input mesh
 * @param pointCloud: output point cloud
 */
void generateVertexPointCloud(TriMesh & mesh,PointCloud &pointCloud);

/**
 * @brief generateFaceCentroidPointCloud: Point cloud generation using face centroids as points
 * @param mesh: input mesh
 * @param pointCloud: output point cloud
 */
void generateFaceCentroidPointCloud(TriMesh & mesh,PointCloud &pointCloud);

/**
 * @brief kd_tree_t: KD-Tree data structure definition with dimension = 3 and based on L2 distance.
 */
typedef KDTreeSingleIndexAdaptor<
    L2_Simple_Adaptor<num_t, PointCloud > ,
    PointCloud,
    3 /* dim */
    > kd_tree_t;

/** @} */ // end of group1



num_t getMeanSquareAngleError(TriMesh &desired_mesh , TriMesh & mesh);

num_t getAreaError(num_t desiredArea, TriMesh & mesh);

num_t getVolError(num_t desiredVol,TriMesh & mesh);

num_t getMeanVertexDistanceError(TriMesh & desired_mesh, kd_tree_t & dmKDTreeVertices, kd_tree_t &dmKDTreeFaces, TriMesh & mesh);

num_t getL2VertexBasedError(TriMesh & desired_mesh, kd_tree_t &dmKDTreeVertices, kd_tree_t & dmKDTreeFaces, TriMesh & mesh);

num_t getMeanSquareAngleError2(TriMesh & desired_mesh, kd_tree_t &dmKDTreeVertices, kd_tree_t & dmKDTreeFaces, TriMesh & mesh);

num_t getMeanQuadricError(TriMesh & desired_mesh, kd_tree_t &dmKDTreeVertices, TriMesh & mesh);

num_t getL2NormalBasedError(TriMesh & desired_mesh, kd_tree_t &dmKDTreeVertices, kd_tree_t & dmKDTreeFaces, TriMesh & mesh);

num_t getMeanTangentialError(TriMesh & desired_mesh, kd_tree_t &dmKDTreeVertices, TriMesh & mesh);

num_t getMeanDiscreteCurvatureError(TriMesh & desired_mesh, kd_tree_t &dmKDTreeVertices,vector<num_t> & dmCurvature, TriMesh & mesh, vector<num_t> & mCurvature);


inline num_t getMeanSquareAngleError(TriMesh & desired_mesh,kd_tree_t & dmKDTreeVertices, kd_tree_t &dmKDTreeFaces,vector<num_t> & dmCurvature,TriMesh & mesh, vector<num_t> & mCurvature)
{
    return getMeanSquareAngleError(desired_mesh , mesh);
}

inline num_t getAreaError(TriMesh & desired_mesh,kd_tree_t & dmKDTreeVertices, kd_tree_t &dmKDTreeFaces,vector<num_t> & dmCurvature,TriMesh & mesh, vector<num_t> & mCurvature)
{
    return getAreaError(getArea(desired_mesh),mesh);
}

inline num_t getVolError(TriMesh & desired_mesh,kd_tree_t & dmKDTreeVertices, kd_tree_t &dmKDTreeFaces,vector<num_t> & dmCurvature,TriMesh & mesh, vector<num_t> & mCurvature)
{
    return getVolError(getVolume(desired_mesh),mesh);
}

inline num_t getMeanVertexDistanceError(TriMesh & desired_mesh,kd_tree_t & dmKDTreeVertices, kd_tree_t &dmKDTreeFaces,vector<num_t> & dmCurvature,TriMesh & mesh, vector<num_t> & mCurvature)
{
    return getMeanVertexDistanceError(desired_mesh, dmKDTreeVertices,dmKDTreeFaces, mesh);
}

inline num_t getL2VertexBasedError(TriMesh & desired_mesh,kd_tree_t & dmKDTreeVertices, kd_tree_t &dmKDTreeFaces,vector<num_t> & dmCurvature,TriMesh & mesh, vector<num_t> & mCurvature)
{
    return getL2VertexBasedError(desired_mesh,dmKDTreeVertices,dmKDTreeFaces,mesh);
}

inline num_t getMeanSquareAngleError2(TriMesh & desired_mesh,kd_tree_t & dmKDTreeVertices, kd_tree_t &dmKDTreeFaces,vector<num_t> & dmCurvature,TriMesh & mesh, vector<num_t> & mCurvature)
{
    return getMeanSquareAngleError2(desired_mesh,dmKDTreeVertices,dmKDTreeFaces,mesh);
}

inline num_t getMeanQuadricError(TriMesh & desired_mesh,kd_tree_t & dmKDTreeVertices, kd_tree_t &dmKDTreeFaces,vector<num_t> & dmCurvature,TriMesh & mesh, vector<num_t> & mCurvature)
{
    return getMeanQuadricError(desired_mesh,dmKDTreeVertices,mesh);
}

inline num_t getL2NormalBasedError(TriMesh & desired_mesh,kd_tree_t & dmKDTreeVertices, kd_tree_t &dmKDTreeFaces,vector<num_t> & dmCurvature,TriMesh & mesh, vector<num_t> & mCurvature)
{
    return getL2NormalBasedError(desired_mesh,dmKDTreeVertices,dmKDTreeFaces,mesh);
}

inline num_t getMeanTangentialError(TriMesh & desired_mesh,kd_tree_t & dmKDTreeVertices, kd_tree_t &dmKDTreeFaces,vector<num_t> & dmCurvature,TriMesh & mesh, vector<num_t> & mCurvature)
{
    return getMeanTangentialError(desired_mesh,dmKDTreeVertices,mesh);
}

inline num_t getMeanDiscreteCurvatureError(TriMesh & desired_mesh,kd_tree_t & dmKDTreeVertices, kd_tree_t &dmKDTreeFaces,vector<num_t> & dmCurvature,TriMesh & mesh, vector<num_t> & mCurvature)
{
    return getMeanDiscreteCurvatureError(desired_mesh,dmKDTreeVertices,dmCurvature,mesh,mCurvature);
}

#endif // METRIC_H
