#ifndef MESH_H
#define MESH_H

#include "Types.h"
#include "Vertex.h"
#include "Edge.h"
#include "Face.h"
#include "HalfEdge.h"
#include <unordered_map>
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>

class Mesh {
public:
    // default constructor
    Mesh();
        
    // read mesh from file
    bool read(const std::string& fileName);
    
    // write mesh to file
    bool write(const std::string& fileName) const;
    
    // deforms the mesh
    void deform(const int iterations);
        
    // member variables
    std::vector<HalfEdge> halfEdges;
    std::vector<Vertex> vertices;
    std::vector<Eigen::Vector3d> uvs;
    std::vector<Eigen::Vector3d> normals;
    std::vector<Edge> edges;
    std::vector<Face> faces;
    std::vector<HalfEdgeIter> boundaries;

private:
    // center mesh about origin and rescale to unit radius
    void normalize();
    
    // builds Laplace Beltrami operator
    void buildLaplacian(Eigen::SparseMatrix<double>& L);
    
    // prefactors Laplace Beltrami operator
    void setup();
    
    // computes rotations using SVD
    void computeRotations();
    
    // member variable
    std::unordered_map<int, double> weights;
    Eigen::MatrixXd deformedCoords;
    std::vector<Eigen::Matrix3d> rotations;
    Eigen::SparseMatrix<double> LT;
    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> solver;
    bool didSetup;
};

#endif