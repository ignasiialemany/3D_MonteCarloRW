#include "polyhedronSet.h"

polyhedronSet::polyhedronSet(Eigen::MatrixXd &vertices_input, Eigen::MatrixXd &faces_input)
{
    _unstrainedvertices = vertices_input;
    _faces = faces_input;
    //Find the mean of each col of vertces_input
    //Sum the first column of vertices_input
    //Sum the second column of vertices_input
    _unstrainedcentroid(0) = vertices_input.col(0).sum()/vertices_input.rows();
    _unstrainedcentroid(1) = vertices_input.col(1).sum()/vertices_input.rows();
    _unstrainedcentroid(2) = vertices_input.col(2).sum()/vertices_input.rows();
}

polyhedronSet::polyhedronSet(const std::string& filename){

    offMeshReader reader;
    reader.read(filename);
    _unstrainedvertices = reader.getVertices();
    _faces = reader.getFaces();
    _unstrainedcentroid(0) = _unstrainedvertices.col(0).sum()/_unstrainedvertices.rows();
    _unstrainedcentroid(1) = _unstrainedvertices.col(1).sum()/_unstrainedvertices.rows();
    _unstrainedcentroid(2) = _unstrainedvertices.col(2).sum()/_unstrainedvertices.rows();
}

void polyhedronSet::precomputePolygons(Eigen::Vector3d &centroid, Eigen::VectorXd &sequence_dt,std::function<double(double)> strain){
    double total_time=0.;
    _strain = strain;
    _polygonSet.reserve(sequence_dt.size()+1);
    for (int i = 0; i < sequence_dt.size()+1; i++)
    {
        double strain = _strain(total_time);
        Eigen::MatrixXd strained_vertices = compute_vertex(strain, centroid);
        polygon polygon_curr = polygon(strained_vertices, _faces);
        _polygonSet.emplace_back(strained_vertices, _faces);
        if (i<sequence_dt.size()){
            total_time += sequence_dt(i);
        }
    }
}
    
Eigen::MatrixXd polyhedronSet::compute_vertex(double strain, Eigen::Vector3d &centroid_strain)
{
    //Create EigenMatrixXd displacement based on number of vertices and 3 columns
    Eigen::MatrixXd strained_vertices(_unstrainedvertices.rows(), 3);
    Eigen::MatrixXd displacement(_unstrainedvertices.rows(), 3);

    //Iterate over rows of _unstrainedvertices
    for (int i=0; i<_unstrainedvertices.rows(); i++)
    {
        Eigen::Vector3d vertex_eigen = _unstrainedvertices.row(i);
        Eigen::Vector3d relative_position = vertex_eigen - centroid_strain;

        //For displacements row component 0 and 1 do not change
        displacement(i, 0) = (1/std::sqrt(1+strain))*relative_position(0) - relative_position(0);
        displacement(i, 1) = (1/std::sqrt(1+strain))*relative_position(1) - relative_position(1);
        displacement(i, 2) = strain*relative_position(2);

        strained_vertices(i, 0) = vertex_eigen(0) + displacement(i, 0);
        strained_vertices(i, 1) = vertex_eigen(1) + displacement(i, 1);
        strained_vertices(i, 2) = vertex_eigen(2) + displacement(i, 2);
    }

    //Iterate over _faces and check that the 3 vertices are not colinear
    for(int j=0; j<_faces.rows();j++){
    Eigen::Vector3d index_vertices = _faces.row(j);
    Eigen::Vector3d vertex1 = strained_vertices.row(index_vertices(0)-1);
    Eigen::Vector3d vertex2 = strained_vertices.row(index_vertices(1)-1);
    Eigen::Vector3d vertex3 = strained_vertices.row(index_vertices(2)-1);
    Eigen::Vector3d edge1 = vertex2 - vertex1;
    Eigen::Vector3d edge2 = vertex3 - vertex1;
    Eigen::Vector3d cross_product = edge1.cross(edge2);
    if (cross_product.norm() < 1e-9) { // you might want to choose a suitable small value for comparison
        double xs =1234;
    }
    }

    return strained_vertices;
}