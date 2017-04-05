#include "denoising.h"

bool hasIND(TriMesh::Point & p)
{
    if(p[0]!=p[0]||p[1]!=p[1]||p[2]!=p[2])
        return true;
    else return false;
}

/////////////////////////////////////////////////
/// Auxiliar functions
/////////////////////////////////////////////////

void updateVertexPositions(TriMesh &mesh, vector<TriMesh::Normal> &filtered_normals, int iteration_number, bool fixed_boundary)
{
    vector<TriMesh::Point> new_points(mesh.n_vertices());

    vector<TriMesh::Point> centroids;

    for(int iter = 0; iter < iteration_number; iter++)
    {
        getAllFaceCentroids(mesh, centroids);
        for(TriMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
        {
            TriMesh::Point p = mesh.point(*v_it);
            if(fixed_boundary && mesh.is_boundary(*v_it))
            {
                new_points.at(v_it->idx()) = p;
            }
            else
            {
                num_t face_num = 0.0;
                TriMesh::Point temp_point(0.0, 0.0, 0.0);
                for(TriMesh::VertexFaceIter vf_it = mesh.vf_iter(*v_it); vf_it.is_valid(); vf_it++)
                {
                    TriMesh::Normal temp_normal = filtered_normals[vf_it->idx()];
                    TriMesh::Point temp_centroid = centroids[vf_it->idx()];
                    temp_point += temp_normal * (temp_normal | (temp_centroid - p));
                    face_num++;
                }
                p += temp_point / face_num;
                if (!hasIND(p))
                {
                    new_points.at(v_it->idx()) = p;
                }
                else new_points.at(v_it->idx()) = p;
                    //cout<<"p: "<<p<<" temp_point: "<<temp_point<<" fn: "<<face_num<<endl;

            }
        }

        for(TriMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
            mesh.set_point(*v_it, new_points[v_it->idx()]);
    }
}

num_t getSigmaC(TriMesh &mesh, vector<TriMesh::Point> &face_centroids, num_t sigma_c_scalar)
{
    num_t sigma_c = 0.0;
    num_t num = 0.0;
    for(TriMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
    {
        TriMesh::Point ci = face_centroids[f_it->idx()];
        for(TriMesh::FaceFaceIter ff_it = mesh.ff_iter(*f_it); ff_it.is_valid(); ff_it++)
        {
            TriMesh::Point cj = face_centroids[ff_it->idx()];
            sigma_c += (ci - cj).length();
            num++;
        }
    }
    sigma_c *= sigma_c_scalar / num;
    return sigma_c;
}

num_t getRadius(TriMesh &mesh , num_t scalar )
{
    vector<TriMesh::Point> centroids;
    getAllFaceCentroids(mesh, centroids);

    num_t radius = 0.0;
    num_t num = 0.0;
    for(TriMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
    {
        TriMesh::Point fi = centroids[f_it->idx()];
        for(TriMesh::FaceFaceIter ff_it = mesh.ff_iter(*f_it); ff_it.is_valid(); ff_it++)
        {
            TriMesh::Point fj = centroids[ff_it->idx()];
            radius += (fj - fi).length();
            num++;
        }
    }
    return radius * scalar / num;
}

/////////////////////////////////////////////////
/////////////////////////////////////////////////
/// Isotropic mesh denoising algorithms
/////////////////////////////////////////////////
/////////////////////////////////////////////////

/////////////////////////////////////////////////
/// Uniform Laplacian smoothing
/////////////////////////////////////////////////

TriMesh uniformLaplacian(TriMesh & _mesh,int iteration_number,num_t scale)
{
    TriMesh mesh = _mesh;
    vector<TriMesh::Point> displacement_points;
    displacement_points.resize(mesh.n_vertices());
    for(int iter = 0; iter < iteration_number; iter++)
    {
        for(TriMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
        {
            vector<TriMesh::VertexHandle> vertex_neighbors;
            getVertexNeighbors(_mesh,*v_it,1,vertex_neighbors);
            num_t weight = 1.0/((num_t)vertex_neighbors.size());
            TriMesh::Point sum(0.0,0.0,0.0);
            for(int i=0;i<vertex_neighbors.size();i++)
            {
                sum=sum+((mesh.point(vertex_neighbors[i])-mesh.point(*v_it))*weight);
            }
           displacement_points[v_it->idx()]=sum;
        }

        for(TriMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
        {
            mesh.set_point(*v_it, displacement_points[v_it->idx()]*scale+mesh.point(*v_it));
        }
    }
    return mesh;
}

/////////////////////////////////////////////////
/// Taubin smoothing (G. Taubin)
/////////////////////////////////////////////////
TriMesh taubin(TriMesh & _mesh,int iteration_number,num_t lambda,num_t mu)
{
    /*TriMesh mesh = _mesh;
    for(int iter = 0; iter < iteration_number; iter++)
    {
        mesh=uniformLaplacian(mesh,1,lambda);
        mesh=uniformLaplacian(mesh,1,mu);
    }
    return mesh;*/
    TriMesh mesh = _mesh;

        vector<TriMesh::Point> displacement_points;
        displacement_points.resize(mesh.n_vertices());
        vector<TriMesh::Point> new_points(mesh.n_vertices());

        for(int it = 0; it < iteration_number; it++)
        {
            for (TriMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
            {
                vector<TriMesh::VertexHandle> vertex_neighbors;
                getVertexNeighbors(_mesh, *v_it, 1, vertex_neighbors);
                num_t weight = 1.0 / ((num_t)vertex_neighbors.size());
                TriMesh::Point sum(0.0, 0.0, 0.0);
                for (int i = 0; i<(int)(vertex_neighbors.size()); i++)
                    sum = sum + ((mesh.point(vertex_neighbors[i]) - mesh.point(*v_it))*weight);
                displacement_points[v_it->idx()] = sum;
            }
            for (TriMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
                mesh.set_point(*v_it, displacement_points[v_it->idx()] * lambda + mesh.point(*v_it));
            for (TriMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
            {
                vector<TriMesh::VertexHandle> vertex_neighbors;
                getVertexNeighbors(_mesh, *v_it, 1, vertex_neighbors);
                num_t weight = 1.0 / ((num_t)vertex_neighbors.size());
                TriMesh::Point sum(0.0, 0.0, 0.0);
                for (int i = 0; i<(int)(vertex_neighbors.size()); i++)
                    sum = sum + ((mesh.point(vertex_neighbors[i]) - mesh.point(*v_it))*weight);
                displacement_points[v_it->idx()] = sum;
            }
            for (TriMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
                mesh.set_point(*v_it, displacement_points[v_it->idx()] * mu + mesh.point(*v_it));
        }
        return mesh;
}

/////////////////////////////////////////////////
/// HC Laplacian smoothing (Vollmer et al.)
/////////////////////////////////////////////////
TriMesh HCLaplacian(TriMesh & _mesh,int iteration_number,num_t alpha,num_t beta)
{
    TriMesh mesh = _mesh;
    vector<TriMesh::Point> original_points;
    getAllPoints(mesh,original_points);
    TriMesh temp_mesh_p;
    vector<TriMesh::Point> temp_points_p(original_points);
    for(int iter = 0; iter < iteration_number; iter++)
    {
        vector<TriMesh::Point> temp_points_q(temp_points_p);
        temp_mesh_p=uniformLaplacian(mesh,1,1.0);

        getAllPoints(temp_mesh_p,temp_points_p);
        vector<TriMesh::Point> temp_points_b(temp_points_p.size());
        for(size_t i = 0;i<temp_points_p.size();i++)
            temp_points_b[i]=temp_points_p[i]-(original_points[i]*alpha+(1.0-alpha)*temp_points_q[i]);
        for(TriMesh::VertexIter v_it = mesh.vertices_begin();v_it!=mesh.vertices_end();v_it++)
        {
            TriMesh::Point t_sum(0,0,0);
            num_t num = 0;
            for(TriMesh::VertexVertexIter vv_it = mesh.vv_iter(*v_it);vv_it.is_valid();vv_it++,num+=1.0)
                t_sum+=temp_points_b[vv_it->idx()];
            temp_points_p[v_it->idx()]=temp_points_p[v_it->idx()]-(beta*temp_points_b[v_it->idx()]+((1.0-beta)/num)*t_sum);
        }

    }

    for(TriMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
        mesh.set_point(*v_it, temp_points_p[v_it->idx()]);
    return mesh;
}

/////////////////////////////////////////////////
/// Hildebrandt and Polthier smoohing
/////////////////////////////////////////////////
num_t HildebrandtAndPolthier_weight_function(num_t lambda,num_t transition_width,num_t value)
{
    //return 1.0;
    if(abs(value)<=lambda)
        return 1.0;
    else
        return (pow(lambda,2.0))/(transition_width*pow(lambda-abs(value),2.0)+pow(lambda,2.0));
}

TriMesh HildebrandtAndPolthier(TriMesh & _mesh,int iteration_number,num_t lambda,num_t transition_width,num_t step_width)
{
    TriMesh mesh = _mesh;
    //vector<num_t> edge_mean_curvatures(mesh.n_edges(),0.0);
   // vector<TriMesh::Normal> edge_normals(mesh.n_edges());
    vector<TriMesh::Point> laplacian;
    vector<num_t> areas;
    mesh.request_face_normals();
    mesh.request_halfedge_normals();

    for(int iter = 0;iter<iteration_number;iter++)
    {

        mesh.update_normals();
        laplacian.resize(mesh.n_vertices(),TriMesh::Point(.0,.0,.0));
        getAllVertexAreas(mesh,areas);
        for(TriMesh::EdgeIter e_it = mesh.edges_begin();e_it != mesh.edges_end();e_it++)
        {
            if(!mesh.is_boundary(*e_it))
            {
                //TriMesh::FaceHandle fh1 = mesh.face_handle(mesh.halfedge_handle(*e_it, 0));
                //TriMesh::FaceHandle fh2 = mesh.opposite_face_handle(mesh.halfedge_handle(*e_it, 0));
                //TriMesh::Normal n1 = mesh.normal(fh1);
                //TriMesh::Normal n2 = mesh.normal(fh2);
                TriMesh::Normal Ne = mesh.normal(mesh.halfedge_handle(*e_it, 0));
                //edge_normals[e_it->idx()] = Ne;
                num_t dihedral_angle = mesh.calc_dihedral_angle(*e_it);
                num_t edge_length = mesh.calc_edge_length(*e_it);
                num_t He = 2*edge_length*cos(dihedral_angle/2.0);
                //edge_mean_curvatures[e_it->idx()] = He;
                TriMesh::VertexHandle v1 = mesh.to_vertex_handle(mesh.halfedge_handle(*e_it, 0));
                TriMesh::VertexHandle v2 = mesh.from_vertex_handle(mesh.halfedge_handle(*e_it, 0));
                num_t weight = HildebrandtAndPolthier_weight_function(lambda,transition_width,He);
                TriMesh::Point temp =(weight*He)*Ne;
                //TriMesh::Point temp = Ne;
                //cout<<"dihedral angle: "<<dihedral_angle<<endl;
                //cout<<"Edge length: "<<edge_length<<" He: "<<He<<" Ne: "<<Ne<<" weight: "<<weight<<" Temp: "<<temp<<endl;
                laplacian[v1.idx()]-=temp;
                laplacian[v2.idx()]-=temp;
            }
        }
        for(TriMesh::VertexIter v_it=mesh.vertices_begin();v_it!=mesh.vertices_end();v_it++)
        {
            TriMesh::Point current_point = mesh.point(*v_it);
            cout<<current_point<<" "<<(laplacian[v_it->idx()])<<" ***** "<<(step_width/(2.0*areas[v_it->idx()]))<<endl;
            current_point += (step_width/(2.0*areas[v_it->idx()]))*laplacian[v_it->idx()];
            mesh.set_point(*v_it,current_point);
        }
    }
    return mesh;
}

/////////////////////////////////////////////////
/// Bilateral mesh denoising (Fleishman et al.)
/////////////////////////////////////////////////

TriMesh bilateral(TriMesh & _mesh,int iteration_number)
{
    TriMesh mesh = _mesh;

    vector<TriMesh::Point> points;
    points.resize(mesh.n_vertices());

    for(int iter = 0; iter < iteration_number; iter++)
    {
        mesh.request_face_normals();
        mesh.request_vertex_normals();
        mesh.update_normals();

        num_t sigma_c, sigma_s;

        int index = 0;
        for(TriMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
        {
            // calculate sigma_c
            TriMesh::Point pi = mesh.point(*v_it);
            TriMesh::Normal ni = mesh.normal(*v_it);
            num_t max = 1e10;
            for(TriMesh::VertexVertexIter vv_it = mesh.vv_iter(*v_it); vv_it.is_valid(); vv_it++)
            {
                TriMesh::Point pj = mesh.point(*vv_it);
                num_t length = (pi - pj).length();
                if(length < max) max = length;
            }
            sigma_c = max;

            vector<TriMesh::VertexHandle> vertex_neighbor;
            getAdaptiveVertexNeighbors(mesh, *v_it, sigma_c, vertex_neighbor);

            // calculate sigma_s
            num_t average_off_set = 0;
            vector<num_t> off_set_dis;
            off_set_dis.clear();
            for(int i = 0; i < (int)vertex_neighbor.size(); i++)
            {
                TriMesh::Point pj = mesh.point(vertex_neighbor[i]);

                num_t t = (pj - pi) | ni;
                t = sqrt(t*t);
                average_off_set += t;
                off_set_dis.push_back(t);
            }
            average_off_set = average_off_set / (num_t)vertex_neighbor.size();
            num_t offset = 0;
            for(int j = 0; j < (int)off_set_dis.size(); j++)
                offset += (off_set_dis[j] - average_off_set) * (off_set_dis[j] - average_off_set);
            offset /= (num_t)off_set_dis.size();

            sigma_s = (sqrt(offset) < 1.0e-12) ? (sqrt(offset) + 1.0e-12) : sqrt(offset);

            num_t sum = 0; num_t normalizer = 0;
            for(int iv = 0; iv < (int)vertex_neighbor.size(); iv++)
            {
                TriMesh::Point pj = mesh.point(vertex_neighbor[iv]);

                num_t t = (pi - pj).length();
                num_t h = (pj - pi) | ni;
                num_t wc = GaussianWeight(t,sigma_c);
                num_t ws = GaussianWeight(h,sigma_s);
                sum += wc * ws * h;
                normalizer += wc * ws;
            }
            points[index] = pi + ni * (sum / normalizer);
            index++;
        }
        index = 0;
        for(TriMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
        {
            mesh.set_point(*v_it, points[index]);
            index++;
        }
    }
    return mesh;
}




/////////////////////////////////////////////////
/// Fast and effective feature-presearving mesh denoising (Sun et al.)
/////////////////////////////////////////////////

void updateFilteredNormalsFast(TriMesh &mesh, int normal_iteration_number, num_t threshold_T,FaceNeighborType fnt,vector<TriMesh::Normal> &filtered_normals)
{
    // get parameter for normal update

    vector< vector<TriMesh::FaceHandle> > all_face_neighbors;
    if(fnt==kVertexBased)
        getAllFaceNeighbors_VertexBased(mesh, true, all_face_neighbors);
    else
        getAllFaceNeighbors_EdgeBased(mesh, true, all_face_neighbors);

    filtered_normals.resize(mesh.n_faces());

    vector<TriMesh::Normal> previous_normals;
    getAllFaceNormals(mesh, previous_normals);

    for(int iter = 0; iter < normal_iteration_number; iter++)
    {
        for(TriMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
        {
            TriMesh::Normal ni = previous_normals[f_it->idx()];
            vector<TriMesh::FaceHandle> face_neighbors = all_face_neighbors[f_it->idx()];
            TriMesh::Normal temp_normal(0.0, 0.0, 0.0);
            for(int i = 0; i <(int)face_neighbors.size(); i++)
            {
                TriMesh::Normal nj = previous_normals[face_neighbors[i].idx()];
                num_t value = (ni | nj) - threshold_T;
                num_t weight = (value > 0) ? value * value : 0;
                temp_normal += nj * weight;
            }
            temp_normal.normalize_cond();
            filtered_normals[f_it->idx()] = temp_normal;
        }
        previous_normals = filtered_normals;
    }
}


TriMesh FastAndEffectiveFeaturePreserving(TriMesh & _mesh, int normal_iteration_number, int vertex_iteration_number, num_t threshold_T, FaceNeighborType fnt)
{
    TriMesh mesh = _mesh;
    if(mesh.n_vertices() == 0)
        return mesh;

    // update face normal
    vector<TriMesh::Normal> filtered_normals;
    updateFilteredNormalsFast(mesh,normal_iteration_number,threshold_T,fnt, filtered_normals);

    // update vertex position
    updateVertexPositions(mesh, filtered_normals, vertex_iteration_number, true);

    return mesh;
}

/////////////////////////////////////////////////
/// Bilateral normal filtering for mesh denoising (Zheng et al.)
/////////////////////////////////////////////////

void updateFilteredNormals(TriMesh &mesh, int normal_iteration_number, num_t sigma_c_scalar, num_t sigma_s, vector<TriMesh::Normal> &filtered_normals)
{
    filtered_normals.resize(mesh.n_faces());

    vector< vector<TriMesh::FaceHandle> > all_face_neighbor;
    getAllFaceNeighbors_VertexBased(mesh,0,all_face_neighbor);
    vector<TriMesh::Normal> previous_normals;
    getAllFaceNormals(mesh, previous_normals);
    vector<num_t> face_areas;
    getAllFaceAreas(mesh, face_areas);
    vector<TriMesh::Point> face_centroids;
    getAllFaceCentroids(mesh, face_centroids);

    num_t sigma_c = getSigmaC(mesh, face_centroids, sigma_c_scalar);

    for(int iter = 0; iter < normal_iteration_number; iter++)
    {
        for(TriMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
        {
            int index_i = f_it->idx();
            TriMesh::Normal ni = previous_normals[index_i];
            TriMesh::Point ci = face_centroids[index_i];
            vector<TriMesh::FaceHandle> face_neighbor = all_face_neighbor[index_i];
            int size = (int)face_neighbor.size();
            TriMesh::Normal temp_normal(0.0, 0.0, 0.0);
            num_t weight_sum = 0.0;
            for(int i = 0; i < size; i++)
            {
                int index_j = face_neighbor[i].idx();
                TriMesh::Normal nj = previous_normals[index_j];
                TriMesh::Point cj = face_centroids[index_j];

                num_t spatial_distance = (ci - cj).length();
                num_t spatial_weight = GaussianWeight(spatial_distance,sigma_c);
                num_t range_distance = (ni - nj).length();
                num_t range_weight = GaussianWeight(range_distance,sigma_s);

                num_t weight = face_areas[index_j] * spatial_weight * range_weight;
                weight_sum += weight;
                temp_normal += nj * weight;
            }
            temp_normal /= weight_sum;
            temp_normal.normalize_cond();
            filtered_normals[index_i] = temp_normal;
        }
        previous_normals = filtered_normals;
    }
}

TriMesh bilateralNormal(TriMesh &_mesh, int normal_iteration_number, int vertex_iteration_number, num_t sigma_c_scalar,  num_t sigma_s )
{
    TriMesh mesh = _mesh;

    if(mesh.n_vertices() == 0)
        return mesh;

    // update face normal
    vector<TriMesh::Normal> filtered_normals;
    updateFilteredNormals(mesh, normal_iteration_number, sigma_c_scalar, sigma_s, filtered_normals);

    updateVertexPositions(mesh, filtered_normals, vertex_iteration_number, true);

    return mesh;
}


/////////////////////////////////////////////////
/// Guided mesh normal filtering (Zhang et al.)
/////////////////////////////////////////////////

void getAllFaceNeighborGMNF(TriMesh &mesh, FaceNeighborType face_neighbor_type, num_t radius, bool include_central_face,
                                                       vector<vector<TriMesh::FaceHandle> > &all_face_neighbor)
{
    vector<TriMesh::FaceHandle> face_neighbor;
    for(TriMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
    {
        if(face_neighbor_type == kVertexBased)
            getFaceNeighbors_VertexBased(mesh, *f_it, face_neighbor);
        else if(face_neighbor_type == kRadiusBased)
            getFaceNeighbors_RadiusBased(mesh, *f_it, radius, face_neighbor);

        if(include_central_face)
            face_neighbor.push_back(*f_it);
        all_face_neighbor[f_it->idx()] = face_neighbor;
    }
}

void getAllGuidedNeighborGMNF(TriMesh &mesh, vector<vector<TriMesh::FaceHandle> > &all_guided_neighbor)
{
    vector<TriMesh::FaceHandle> face_neighbor;
    for(TriMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
    {
        getFaceNeighbors_VertexBased(mesh, *f_it, face_neighbor);
        face_neighbor.push_back(*f_it);
        all_guided_neighbor[f_it->idx()] = face_neighbor;
    }
}

void getFaceNeighborInnerEdge(TriMesh &mesh, vector<TriMesh::FaceHandle> &face_neighbor, vector<TriMesh::EdgeHandle> &inner_edge)
{
    inner_edge.clear();
    vector<bool> edge_flag((int)mesh.n_edges(), false);
    vector<bool> face_flag((int)mesh.n_faces(), false);

    for(int i = 0; i < (int)face_neighbor.size(); i++)
        face_flag[face_neighbor[i].idx()] = true;

    for(int i = 0; i < (int)face_neighbor.size(); i++)
    {
        for(TriMesh::FaceEdgeIter fe_it = mesh.fe_iter(face_neighbor[i]); fe_it.is_valid(); fe_it++)
        {
            if((!edge_flag[fe_it->idx()]) && (!mesh.is_boundary(*fe_it)))
            {
                edge_flag[fe_it->idx()] = true;
                TriMesh::HalfedgeHandle heh = mesh.halfedge_handle(*fe_it, 0);
                TriMesh::FaceHandle f = mesh.face_handle(heh);
                TriMesh::HalfedgeHandle heho = mesh.opposite_halfedge_handle(heh);
                TriMesh::FaceHandle fo = mesh.face_handle(heho);
                if(face_flag[f.idx()] && face_flag[fo.idx()])
                    inner_edge.push_back(*fe_it);
            }
        }
    }
}

void getRangeAndMeanNormal(TriMesh &mesh, vector<vector<TriMesh::FaceHandle> > &all_guided_neighbor,
                                                       vector<num_t> &face_areas, vector<TriMesh::Normal> &normals,
                                                       vector<pair<num_t, TriMesh::Normal> > &range_and_mean_normal)
{
    const num_t epsilon = 1.0e-9;

    for(TriMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
    {
        int index = f_it->idx();
        vector<TriMesh::FaceHandle> face_neighbor = all_guided_neighbor[index];
        num_t metric = 0.0;
        TriMesh::Normal average_normal(0.0, 0.0, 0.0);
        num_t maxdiff = -1.0;

        for(int i = 0; i < (int)face_neighbor.size(); i++)
        {
            int index_i = face_neighbor[i].idx();
            num_t area_weight = face_areas[index_i];
            TriMesh::Normal ni = normals[index_i];
            average_normal += ni * area_weight;

            for(int j = i+1; j < (int)face_neighbor.size(); j++)
            {
                int index_j = face_neighbor[j].idx();
                TriMesh::Normal nj = normals[index_j];
                num_t diff = NormalDistance(ni, nj);

                if(diff > maxdiff)
                {
                    maxdiff = diff;
                }
            }
        }

        vector<TriMesh::EdgeHandle> inner_edge_handle;
        getFaceNeighborInnerEdge(mesh, face_neighbor, inner_edge_handle);
        num_t sum_tv = 0.0, max_tv = -1.0;
        for(int i = 0; i < (int)inner_edge_handle.size(); i++)
        {
            TriMesh::HalfedgeHandle heh = mesh.halfedge_handle(inner_edge_handle[i], 0);
            TriMesh::FaceHandle f = mesh.face_handle(heh);
            TriMesh::Normal n1 = normals[f.idx()];
            TriMesh::HalfedgeHandle heho = mesh.opposite_halfedge_handle(heh);
            TriMesh::FaceHandle fo = mesh.face_handle(heho);
            TriMesh::Normal n2 = normals[fo.idx()];
            num_t current_tv = NormalDistance(n1, n2);
            max_tv = (current_tv > max_tv)? current_tv : max_tv;
            sum_tv += current_tv;
        }
        /*TriMesh::Normal temp = average_normal;
        temp.normalize_cond();*/
        average_normal.normalize_cond();
        /*if (!hasIND(temp))
        {
            average_normal = temp;
            //cout<<"average_normal "<<average_normal<<endl;
            //return;
        }
        else
        {
            cout<<" size neighborhood: "<<face_neighbor.size()<<" average_normal "<<average_normal<<endl;
            for(TriMesh::FaceVertexIter fv_it = mesh.fv_begin(*f_it); fv_it.is_valid(); fv_it++)
                cout<<mesh.point(*fv_it)<<" // ";
            cout<<endl;
        }*/
        metric = maxdiff * max_tv / (sum_tv + epsilon);

        range_and_mean_normal[index] = make_pair(metric, average_normal);
    }
}

void getGuidedNormals(TriMesh &mesh, vector<vector<TriMesh::FaceHandle> > &all_guided_neighbor,
                                                 vector<num_t> &face_areas, vector<TriMesh::Normal> &normals,
                                                 vector<pair<num_t, TriMesh::Normal> > range_and_mean_normal,
                                                 vector<TriMesh::Normal> &guided_normals)
{
    getRangeAndMeanNormal(mesh, all_guided_neighbor, face_areas, normals, range_and_mean_normal);

    for(TriMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
    {
        vector<TriMesh::FaceHandle> face_neighbor = all_guided_neighbor[f_it->idx()];
        num_t min_range = 1.0e8;
        int min_idx = 0;
        for(int i = 0; i < (int)face_neighbor.size(); i++)
        {
            num_t current_range = range_and_mean_normal[face_neighbor[i].idx()].first;
            if(min_range > current_range){
                min_range = current_range;
                min_idx = i;
            }
        }
        TriMesh::FaceHandle min_face_handle = face_neighbor[min_idx];
        guided_normals[f_it->idx()] = range_and_mean_normal[min_face_handle.idx()].second;
    }
}







void updateFilteredNormalsGuided(TriMesh &mesh, vector<TriMesh::Normal> &filtered_normals, num_t radius_scalar,
                                 num_t sigma_c_scalar, int normal_iteration_number, num_t sigma_s, int vertex_iteration_number)
{
    filtered_normals.resize((int)mesh.n_faces());
    // get parameter for local scheme normal update
    int face_neighbor_index = kRadiusBased;
    bool include_central_face = 1;

    /*radius_scalar = 2.0;
    sigma_c_scalar = 1.0;
    normal_iteration_number = 6;
    sigma_s = 0.35;
    vertex_iteration_number = 10;*/

    FaceNeighborType face_neighbor_type = face_neighbor_index == 0 ? kRadiusBased : kVertexBased;

    num_t radius;
    if(face_neighbor_type == kRadiusBased)
        radius = getRadius(mesh,radius_scalar);

    vector<vector<TriMesh::FaceHandle> > all_face_neighbor((int)mesh.n_faces());
    getAllFaceNeighborGMNF(mesh, face_neighbor_type, radius, include_central_face, all_face_neighbor);
    vector<vector<TriMesh::FaceHandle> > all_guided_neighbor((int)mesh.n_faces());
    getAllGuidedNeighborGMNF(mesh, all_guided_neighbor);
    getAllFaceNormals(mesh, filtered_normals);

    vector<num_t> face_areas((int)mesh.n_faces());
    vector<TriMesh::Point> face_centroids((int)mesh.n_faces());
    vector<TriMesh::Normal> previous_normals((int)mesh.n_faces());
    vector<TriMesh::Normal> guided_normals((int)mesh.n_faces());
    vector<pair<num_t, TriMesh::Normal> > range_and_mean_normal((int)mesh.n_faces());
    for(int iter = 0; iter < normal_iteration_number; iter++)
    {
        getAllFaceCentroids(mesh, face_centroids);
        num_t sigma_c = getSigmaC(mesh, face_centroids, sigma_c_scalar);
        getAllFaceAreas(mesh, face_areas);
        getAllFaceNormals(mesh, previous_normals);

        getGuidedNormals(mesh, all_guided_neighbor, face_areas, previous_normals, range_and_mean_normal, guided_normals);

        for(TriMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
        {
            int index = f_it->idx();
            const vector<TriMesh::FaceHandle> face_neighbor = all_face_neighbor[index];
            TriMesh::Normal filtered_normal(0.0, 0.0, 0.0);
            for(int j = 0; j < (int)face_neighbor.size(); j++)
            {
                int current_face_index = face_neighbor[j].idx();

                num_t spatial_dis = (face_centroids[index] - face_centroids[current_face_index]).length();
                num_t spatial_weight = GaussianWeight(spatial_dis, sigma_c);
                num_t range_dis = (guided_normals[index] - guided_normals[current_face_index]).length();
                num_t range_weight = GaussianWeight(range_dis, sigma_s);

                filtered_normal += previous_normals[current_face_index] * (face_areas[current_face_index] * spatial_weight * range_weight);
            }
            if(face_neighbor.size())
                filtered_normals[index] = filtered_normal.normalize_cond();
        }


        updateVertexPositions(mesh, filtered_normals, vertex_iteration_number, false);
    }
}

TriMesh guided(TriMesh _mesh,int normal_iteration_number, int vertex_iteration_number, num_t sigma_c_scalar, num_t sigma_s, num_t radius_scalar)
{
    // get data
    TriMesh mesh = _mesh;

    if(mesh.n_vertices() == 0)
        return mesh;

    // update face normal
    vector<TriMesh::Normal> filtered_normals;
    updateFilteredNormalsGuided(mesh, filtered_normals, radius_scalar, sigma_c_scalar, normal_iteration_number, sigma_s, vertex_iteration_number);

    // get parameter for vertex update
    //vertex_iteration_number = 10;

    updateVertexPositions(mesh, filtered_normals, vertex_iteration_number, true);

    return mesh;
}
