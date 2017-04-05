#include "iomesh.h"

bool importMesh(TriMesh &mesh, string filename)
{
    OpenMesh::IO::Options opt;
    //opt += OpenMesh::IO::Options::VertexColor;
    //opt += OpenMesh::IO::Options::ColorFloat;
    if(!OpenMesh::IO::read_mesh(mesh, filename, opt))
        return false;
    return true;
}

bool exportMesh(TriMesh &mesh, string filename)
{
    OpenMesh::IO::Options opt;
    //opt += OpenMesh::IO::Options::VertexColor;
    //opt += OpenMesh::IO::Options::ColorFloat;
    return OpenMesh::IO::write_mesh(mesh, filename,opt);
}
