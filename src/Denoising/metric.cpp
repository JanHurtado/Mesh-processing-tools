#include "metric.h"

const num_t average_edge_length_factor = 2.0; //KD-Tree nearest neighbor tolerance radius = factor*average_edge_length

void generateVertexPointCloud(TriMesh & mesh, PointCloud &pointCloud)
{
    pointCloud.pts.resize(mesh.n_vertices());
    for(TriMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
    {
        pointCloud.pts[v_it->idx()].x = mesh.point(*v_it)[0];
        pointCloud.pts[v_it->idx()].y = mesh.point(*v_it)[1];
        pointCloud.pts[v_it->idx()].z = mesh.point(*v_it)[2];
    }
}

void generateFaceCentroidPointCloud(TriMesh & mesh,PointCloud &pointCloud)
{
    pointCloud.pts.resize(mesh.n_faces());
    vector<TriMesh::Point> face_centroids;
    getAllFaceCentroids(mesh, face_centroids);

    for (size_t i=0;i<face_centroids.size();i++)
    {
        pointCloud.pts[i].x = face_centroids[i][0];
        pointCloud.pts[i].y = face_centroids[i][1];
        pointCloud.pts[i].z = face_centroids[i][2];
    }
}



num_t getAreaError(num_t desiredArea, TriMesh & mesh)
{
    num_t dif = abs(getArea(mesh) - desiredArea);
    return dif/desiredArea;
}



num_t getVolError(num_t desiredVol,TriMesh & mesh)
{
    num_t dif = abs(getVolume(mesh) - desiredVol);
    return dif/desiredVol;
}

num_t point_triangle_dist(TriMesh & _mesh,TriMesh::FaceHandle fh,TriMesh::Point _point)
{
    TriMesh::Point _p = _point;
    TriMesh::FaceVertexIter fv_iter = _mesh.fv_iter(fh);
    TriMesh::Point _v0 = _mesh.point(*fv_iter);
    fv_iter++;
    TriMesh::Point _v1 = _mesh.point(*fv_iter);
    fv_iter++;
    TriMesh::Point _v2 = _mesh.point(*fv_iter);
    {
      const TriMesh::Point v0v1 = _v1 - _v0;
      const TriMesh::Point v0v2 = _v2 - _v0;
      const TriMesh::Point n = v0v1 % v0v2; // not normalized !
      const num_t d = n.sqrnorm();
      // Check if the triangle is degenerated
      if (d < FLT_MIN && d > -FLT_MIN) {
        return -1.0;
      }
      const num_t invD = 1.0 / d;
      // these are not needed for every point, should still perform
      // better with many points against one triangle
      const TriMesh::Point v1v2 = _v2 - _v1;
      const num_t inv_v0v2_2 = 1.0 / v0v2.sqrnorm();
      const num_t inv_v0v1_2 = 1.0 / v0v1.sqrnorm();
      const num_t inv_v1v2_2 = 1.0 / v1v2.sqrnorm();
      TriMesh::Point v0p = _p - _v0;
      TriMesh::Point t = v0p % n;
      num_t  s01, s02, s12;
      const num_t a = (t | v0v2) * -invD;
      const num_t b = (t | v0v1) * invD;
      if (a < 0)
      {
        // Calculate the distance to an edge or a corner vertex
        s02 = ( v0v2 | v0p ) * inv_v0v2_2;
        if (s02 < 0.0)
        {
          s01 = ( v0v1 | v0p ) * inv_v0v1_2;
          if (s01 <= 0.0) {
            v0p = _v0;
          } else if (s01 >= 1.0) {
            v0p = _v1;
          } else {
            v0p = _v0 + v0v1 * s01;
          }
        } else if (s02 > 1.0) {
          s12 = ( v1v2 | ( _p - _v1 )) * inv_v1v2_2;
          if (s12 >= 1.0) {
            v0p = _v2;
          } else if (s12 <= 0.0) {
            v0p = _v1;
          } else {
            v0p = _v1 + v1v2 * s12;
          }
        } else {
          v0p = _v0 + v0v2 * s02;
        }
      } else if (b < 0.0) {
        // Calculate the distance to an edge or a corner vertex
        s01 = ( v0v1 | v0p ) * inv_v0v1_2;
        if (s01 < 0.0)
        {
      const TriMesh::Point n = v0v1 % v0v2; // not normalized !
          s02 = ( v0v2 |  v0p ) * inv_v0v2_2;
          if (s02 <= 0.0) {
            v0p = _v0;
          } else if (s02 >= 1.0) {
            v0p = _v2;
          } else {
            v0p = _v0 + v0v2 * s02;
          }
        } else if (s01 > 1.0) {
          s12 = ( v1v2 | ( _p - _v1 )) * inv_v1v2_2;
          if (s12 >= 1.0) {
            v0p = _v2;
          } else if (s12 <= 0.0) {
            v0p = _v1;
          } else {
            v0p = _v1 + v1v2 * s12;
          }
        } else {
          v0p = _v0 + v0v1 * s01;
        }
      } else if (a+b > 1.0) {
        // Calculate the distance to an edge or a corner vertex
        s12 = ( v1v2 | ( _p - _v1 )) * inv_v1v2_2;
        if (s12 >= 1.0) {
          s02 = ( v0v2 | v0p ) * inv_v0v2_2;
          if (s02 <= 0.0) {
            v0p = _v0;
          } else if (s02 >= 1.0) {
            v0p = _v2;
          } else {
            v0p = _v0 + v0v2*s02;
          }
        } else if (s12 <= 0.0) {
          s01 = ( v0v1 |  v0p ) * inv_v0v1_2;
          if (s01 <= 0.0) {
            v0p = _v0;
          } else if (s01 >= 1.0) {
            v0p = _v1;
          } else {
            v0p = _v0 + v0v1 * s01;
          }
        } else {
          v0p = _v1 + v1v2 * s12;
        }
      } else {
        // Calculate the distance to an interior point of the triangle
        return ( (_p - n*((n|v0p) * invD)) - _p).sqrnorm();
      }
      return (v0p - _p).sqrnorm();
    }
}

num_t point_triangle_minDist(TriMesh & _mesh,set<size_t> & _targetFaces,TriMesh::Point _point,TriMesh::FaceHandle & min_dist_fh)
{
    num_t min=9e10;
    for(set<size_t>::iterator it=_targetFaces.begin();it!=_targetFaces.end();it++)
    {
        TriMesh::FaceHandle fh((int)(*it));
        num_t d = point_triangle_dist(_mesh,fh,_point);
        if(d<min)
        {
            min=d;
            min_dist_fh=fh;
        }
    }
    return min;
}

num_t getMeanVertexDistanceError(TriMesh & desired_mesh, kd_tree_t &dmKDTreeVertices, kd_tree_t & dmKDTreeFaces, TriMesh & mesh)
{
    desired_mesh.request_face_normals();
    desired_mesh.update_normals();
    vector<num_t> mesh_vertex_areas;
    getAllVertexAreas(mesh,mesh_vertex_areas);

    num_t search_radius = average_edge_length_factor*getAverageEdgeLength(desired_mesh);
    num_t mesh_total_area = 0.0;
    num_t error_sum = 0.0;
    num_t error_max = -9e10;
    for(TriMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
    {
        num_t query_pt[3] = {};

        query_pt[0]=mesh.point(*v_it)[0];
        query_pt[1]=mesh.point(*v_it)[1];
        query_pt[2]=mesh.point(*v_it)[2];
        vector<pair<size_t,num_t> >   ret_matches_v;
        vector<pair<size_t,num_t> >   ret_matches_f;
        nanoflann::SearchParams params;
                //params.sorted = false;

        dmKDTreeVertices.radiusSearch(&query_pt[0],search_radius, ret_matches_v, params);
        dmKDTreeFaces.radiusSearch(&query_pt[0],search_radius, ret_matches_f, params);

        set<size_t> targetFaces;
        for(int i=0;i<ret_matches_f.size();i++)
            targetFaces.insert(ret_matches_f[i].first);
        for(int i=0;i<ret_matches_v.size();i++)
        {
            TriMesh::VertexHandle vh(ret_matches_v[i].first);
            for(TriMesh::VertexFaceIter vf_it = mesh.vf_iter(vh);vf_it.is_valid();vf_it++)
            {
                targetFaces.insert(vf_it->idx());
            }
        }
        if (targetFaces.size()>0)
        {
            TriMesh::FaceHandle min_dist_fh;
            num_t min_dist = point_triangle_minDist(desired_mesh,targetFaces,mesh.point(*v_it),min_dist_fh);
            mesh_total_area+=mesh_vertex_areas[v_it->idx()];
            error_sum+=min_dist*mesh_vertex_areas[v_it->idx()];
            if(min_dist>error_max)
                error_max=min_dist;
        }
    }
    //cout<<endl<<error_sum/mesh_total_area<<" "<<error_max;
    return error_sum/mesh_total_area;
}

num_t getL2VertexBasedError(TriMesh & desired_mesh, kd_tree_t &dmKDTreeVertices, kd_tree_t & dmKDTreeFaces, TriMesh & mesh)
{
    vector<num_t> mesh_vertex_areas;
    getAllVertexAreas(mesh,mesh_vertex_areas);

    num_t search_radius = average_edge_length_factor*getAverageEdgeLength(desired_mesh);
    num_t mesh_total_area = 0.0;
    num_t error_sum = 0.0;
    for(TriMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
    {
        num_t query_pt[3] = {};

        query_pt[0]=mesh.point(*v_it)[0];
        query_pt[1]=mesh.point(*v_it)[1];
        query_pt[2]=mesh.point(*v_it)[2];
        vector<pair<size_t,num_t> >   ret_matches_v;
        vector<pair<size_t,num_t> >   ret_matches_f;
        nanoflann::SearchParams params;

        dmKDTreeVertices.radiusSearch(&query_pt[0],search_radius, ret_matches_v, params);
        dmKDTreeFaces.radiusSearch(&query_pt[0],search_radius, ret_matches_f, params);

        set<size_t> targetFaces;
        for(int i=0;i<ret_matches_f.size();i++)
            targetFaces.insert(ret_matches_f[i].first);
        for(int i=0;i<ret_matches_v.size();i++)
        {
            TriMesh::VertexHandle vh((int)(ret_matches_v[i].first));
            for(TriMesh::VertexFaceIter vf_it = mesh.vf_iter(vh);vf_it.is_valid();vf_it++)
            {
                targetFaces.insert(vf_it->idx());
            }
        }
        if (targetFaces.size()>0)
        {
            TriMesh::FaceHandle min_dist_fh;
            num_t min_dist = point_triangle_minDist(desired_mesh,targetFaces,mesh.point(*v_it),min_dist_fh);
            mesh_total_area+=mesh_vertex_areas[v_it->idx()];
            error_sum+=pow(min_dist,2.0)*mesh_vertex_areas[v_it->idx()];
        }
    }
    return sqrt(error_sum/mesh_total_area);

}

num_t getMeanQuadricError(TriMesh & desired_mesh, kd_tree_t &dmKDTreeVertices, TriMesh & mesh)
{
    desired_mesh.request_face_normals();
    desired_mesh.request_vertex_normals();
    desired_mesh.update_normals();
    vector<num_t> desired_mesh_face_areas;
    vector<num_t> mesh_vertex_areas;
    getAllFaceAreas(desired_mesh,desired_mesh_face_areas);
    getAllVertexAreas(mesh,mesh_vertex_areas);
    num_t error_sum = 0.0;
    num_t total_area = 0.0;
    for(TriMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
    {
        num_t query_pt[3] = {};

        query_pt[0]=mesh.point(*v_it)[0];
        query_pt[1]=mesh.point(*v_it)[1];
        query_pt[2]=mesh.point(*v_it)[2];

        const size_t num_results = 1;
        size_t ret_index;
        num_t out_dist_sqr;
        nanoflann::KNNResultSet<num_t> resultSet(num_results);
        resultSet.init(&ret_index, &out_dist_sqr );
        dmKDTreeVertices.findNeighbors(resultSet, &query_pt[0], nanoflann::SearchParams(10));
        num_t vertex_quadric = 0.0;
        TriMesh::VertexHandle nearest_vertex((int)ret_index);
        for(TriMesh::VertexFaceIter vf_it = desired_mesh.vf_iter(nearest_vertex);vf_it.is_valid();vf_it++)
        {
            TriMesh::Point p = mesh.point(*v_it);
            TriMesh::FaceHandle fh = *vf_it;
            TriMesh::Normal face_normal = desired_mesh.normal(fh);
            face_normal.normalize_cond();

            TriMesh::FaceVertexIter fv_iter = desired_mesh.fv_iter(fh);
            TriMesh::Point p0 = desired_mesh.point(*fv_iter);
            num_t d = -(face_normal|p0);

            num_t face_quadric = pow((face_normal|p)+d,2.0);
            vertex_quadric += face_quadric*desired_mesh_face_areas[vf_it->idx()];
        }
        error_sum += mesh_vertex_areas[v_it->idx()]*vertex_quadric;
        total_area+=mesh_vertex_areas[v_it->idx()];
    }
    return error_sum/total_area;
}

num_t getMeanSquareAngleError(TriMesh &desired_mesh,TriMesh & mesh)
{
    desired_mesh.request_face_normals();
    desired_mesh.update_face_normals();

    mesh.request_face_normals();
    mesh.update_face_normals();
    num_t mean_square_angle_error = 0.0;
    for (TriMesh::FaceIter f_it1 = mesh.faces_begin(), f_it2 = desired_mesh.faces_begin();
        f_it1 != mesh.faces_end(); f_it1++, f_it2++ )
    {
        TriMesh::Normal normal1 = mesh.normal(*f_it1);
        TriMesh::Normal normal2 = desired_mesh.normal(*f_it2);
        num_t dot_product = normal1 | normal2;
        dot_product = min(1.0, max(dot_product, -1.0));
        mean_square_angle_error += acos(dot_product) * 180.0 / M_PI;
    }
    return mean_square_angle_error / (num_t)desired_mesh.n_faces();
}

num_t getMeanSquareAngleError2(TriMesh & desired_mesh, kd_tree_t &dmKDTreeVertices, kd_tree_t & dmKDTreeFaces, TriMesh & mesh)
{
    desired_mesh.request_face_normals();
    desired_mesh.update_normals();
    mesh.request_face_normals();
    mesh.update_normals();
    vector<TriMesh::Point> centroids;
    getAllFaceCentroids(mesh,centroids);

    num_t search_radius = average_edge_length_factor*getAverageEdgeLength(desired_mesh);
    num_t num_vert = 0.0;
    num_t error_sum = 0.0;
    for(size_t i=0;i<centroids.size();i++)
    {
        num_t query_pt[3] = {};

        query_pt[0]=centroids[i][0];
        query_pt[1]=centroids[i][1];
        query_pt[2]=centroids[i][2];
        vector<pair<size_t,num_t> >   ret_matches_v;
        vector<pair<size_t,num_t> >   ret_matches_f;
        nanoflann::SearchParams params;

        dmKDTreeVertices.radiusSearch(&query_pt[0],search_radius, ret_matches_v, params);
        dmKDTreeFaces.radiusSearch(&query_pt[0],search_radius, ret_matches_f, params);

        set<size_t> targetFaces;
        for(int j=0;j<ret_matches_f.size();j++)
            targetFaces.insert(ret_matches_f[j].first);
        for(int j=0;j<ret_matches_v.size();j++)
        {
            TriMesh::VertexHandle vh((int)(ret_matches_v[j].first));
            for(TriMesh::VertexFaceIter vf_it = mesh.vf_iter(vh);vf_it.is_valid();vf_it++)
            {
                targetFaces.insert(vf_it->idx());
            }
        }
        if (targetFaces.size()>0)
        {
            TriMesh::FaceHandle min_dist_fh;
            point_triangle_minDist(desired_mesh,targetFaces,centroids[i],min_dist_fh);
            num_vert+=1.0;
            TriMesh::Normal normal1 = mesh.normal(TriMesh::FaceHandle((int)i));
            TriMesh::Normal normal2 = desired_mesh.normal(min_dist_fh);
            num_t dot_product = normal1 | normal2;
            dot_product = min(1.0, max(dot_product, -1.0));
            error_sum+=acos(dot_product) * 180.0 / M_PI;;
        }
    }
    return error_sum/num_vert;
}

num_t getL2NormalBasedError(TriMesh & desired_mesh, kd_tree_t &dmKDTreeVertices, kd_tree_t & dmKDTreeFaces, TriMesh & mesh)
{
    desired_mesh.request_face_normals();
    desired_mesh.update_normals();
    mesh.request_face_normals();
    mesh.update_face_normals();
    vector<num_t> mesh_face_areas;
    getAllFaceAreas(desired_mesh,mesh_face_areas);
    vector<TriMesh::Point> centroids;
    getAllFaceCentroids(mesh,centroids);

    num_t search_radius = average_edge_length_factor*getAverageEdgeLength(desired_mesh);
    num_t mesh_total_area = 0.0;
    num_t error_sum = 0.0;
    for(size_t i=0;i<centroids.size();i++)
    {
        num_t query_pt[3] = {};

        query_pt[0]=centroids[i][0];
        query_pt[1]=centroids[i][1];
        query_pt[2]=centroids[i][2];
        vector<pair<size_t,num_t> >   ret_matches_v;
        vector<pair<size_t,num_t> >   ret_matches_f;
        nanoflann::SearchParams params;

        dmKDTreeVertices.radiusSearch(&query_pt[0],search_radius, ret_matches_v, params);
        dmKDTreeFaces.radiusSearch(&query_pt[0],search_radius, ret_matches_f, params);

        set<size_t> targetFaces;
        for(int j=0;j<ret_matches_f.size();j++)
            targetFaces.insert(ret_matches_f[j].first);
        for(int j=0;j<ret_matches_v.size();j++)
        {
            TriMesh::VertexHandle vh((int)(ret_matches_v[j].first));
            for(TriMesh::VertexFaceIter vf_it = mesh.vf_iter(vh);vf_it.is_valid();vf_it++)
            {
                targetFaces.insert(vf_it->idx());
            }
        }
        if (targetFaces.size()>0)
        {
            TriMesh::FaceHandle min_dist_fh;
            point_triangle_minDist(desired_mesh,targetFaces,centroids[i],min_dist_fh);
            mesh_total_area+=mesh_face_areas[i];
            TriMesh::Normal normal1 = mesh.normal(TriMesh::FaceHandle((int)i));
            TriMesh::Normal normal2 = desired_mesh.normal(min_dist_fh);
            num_t temp = pow((normal1 - normal2).norm(),2.0);
            error_sum+=temp*mesh_face_areas[i];
        }
    }
    return error_sum/mesh_total_area;
}

num_t getMeanTangentialError(TriMesh & desired_mesh, kd_tree_t &dmKDTreeVertices, TriMesh & mesh)
{
    desired_mesh.request_face_normals();
    desired_mesh.request_vertex_normals();
    desired_mesh.update_normals();
    mesh.request_face_normals();
    mesh.request_vertex_normals();
    mesh.update_normals();
    vector<num_t> mesh_vertex_areas;
    getAllVertexAreas(mesh,mesh_vertex_areas);
    num_t error_sum = 0.0;
    num_t mesh_total_area = 0.0;
    for(TriMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
    {
        num_t query_pt[3] = {};

        query_pt[0]=mesh.point(*v_it)[0];
        query_pt[1]=mesh.point(*v_it)[1];
        query_pt[2]=mesh.point(*v_it)[2];

        const size_t num_results = 1;
        size_t ret_index;
        num_t out_dist_sqr;
        nanoflann::KNNResultSet<num_t> resultSet(num_results);
        resultSet.init(&ret_index, &out_dist_sqr );
        dmKDTreeVertices.findNeighbors(resultSet, &query_pt[0], nanoflann::SearchParams(10));


        TriMesh::Normal normal1 = mesh.normal(*v_it);
        TriMesh::Normal normal2 = desired_mesh.normal(TriMesh::VertexHandle((int)ret_index));
        num_t temp = pow((normal1 - normal2).norm(),2.0);

        mesh_total_area+=mesh_vertex_areas[v_it->idx()];
        error_sum+=temp*mesh_vertex_areas[v_it->idx()];
    }
    return error_sum/mesh_total_area;
}

num_t getMeanDiscreteCurvatureError(TriMesh & desired_mesh, kd_tree_t &dmKDTreeVertices,vector<num_t> & dmCurvature, TriMesh & mesh, vector<num_t> & mCurvature)
{
    desired_mesh.request_face_normals();
    desired_mesh.request_vertex_normals();
    desired_mesh.update_normals();
    mesh.request_face_normals();
    mesh.request_vertex_normals();
    mesh.update_normals();
    vector<num_t> mesh_vertex_areas;
    getAllVertexAreas(mesh,mesh_vertex_areas);
    num_t error_sum = 0.0;
    num_t mesh_total_area = 0.0;
    for(TriMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
    {
        num_t query_pt[3] = {};

        query_pt[0]=mesh.point(*v_it)[0];
        query_pt[1]=mesh.point(*v_it)[1];
        query_pt[2]=mesh.point(*v_it)[2];

        const size_t num_results = 1;
        size_t ret_index;
        num_t out_dist_sqr;
        nanoflann::KNNResultSet<num_t> resultSet(num_results);
        resultSet.init(&ret_index, &out_dist_sqr );
        dmKDTreeVertices.findNeighbors(resultSet, &query_pt[0], nanoflann::SearchParams(10));

        num_t temp = abs(dmCurvature[ret_index]-mCurvature[v_it->idx()]);

        mesh_total_area+=mesh_vertex_areas[v_it->idx()];
        error_sum+=temp*mesh_vertex_areas[v_it->idx()];
    }
    return error_sum/mesh_total_area;
}
