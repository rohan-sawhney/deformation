#include "Mesh.h"
#include "MeshIO.h"
#include <Eigen/SVD>

Mesh::Mesh()
{
    
}

bool Mesh::read(const std::string& fileName)
{
    std::ifstream in(fileName.c_str());

    if (!in.is_open()) {
        std::cerr << "Error: Could not open file for reading" << std::endl;
        return false;
    }
    
    bool readSuccessful = false;
    if ((readSuccessful = MeshIO::read(in, *this))) {
        normalize();
        setup();
    }
    
    return readSuccessful;
}

bool Mesh::write(const std::string& fileName) const
{
    std::ofstream out(fileName.c_str());
    
    if (!out.is_open()) {
        std::cerr << "Error: Could not open file for writing" << std::endl;
        return false;
    }
    
    MeshIO::write(out, *this);
    
    return false;
}

void Mesh::buildLaplacian(Eigen::SparseMatrix<double>& L)
{
    std::vector<Eigen::Triplet<double>> LTriplet;
    
    for (VertexCIter v = vertices.begin(); v != vertices.end(); v++) {
    
        double sumCoefficients = 0.0;
        if (!v->anchor && !v->handle) {
            
            HalfEdgeCIter he = v->he;
            do {
                // (cotA + cotB) / 2
                double coefficient = 0.5 * (he->cotan() + he->flip->cotan());
                sumCoefficients += coefficient;
                
                LTriplet.push_back(Eigen::Triplet<double>(v->index, he->flip->vertex->index, -coefficient));
                weights(v->index, he->flip->vertex->index) = coefficient;
                
                he = he->flip->next;
            } while (he != v->he);
            
        } else {
            sumCoefficients = 1.0;
        }
        
        LTriplet.push_back(Eigen::Triplet<double>(v->index, v->index, sumCoefficients));
    }
    
    L.setFromTriplets(LTriplet.begin(), LTriplet.end());
}

void Mesh::setup()
{
    int v = (int)vertices.size();
    weights.resize(v, v);
    weights.setZero();
    
    deformedCoords.resize(v, 3);
    deformedCoords.setZero();
    
    rotations.clear();
    rotations.resize(v, Eigen::Matrix3d::Identity());
    
    Eigen::SparseMatrix<double> L(v, v);
    buildLaplacian(L);
    
    LT = L.transpose();
    solver.compute(LT*L);
}

void Mesh::computeRotations()
{
    Eigen::Matrix3d id = Eigen::Matrix3d::Identity();
    for (VertexCIter v = vertices.begin(); v != vertices.end(); v++) {
        
        int column = 0;
        int valence = v->valence();
        Eigen::MatrixXd P(3, valence);
        Eigen::MatrixXd Pdash(3, valence);
        
        HalfEdgeCIter he = v->he;
        do {
            VertexCIter v2 = he->flip->vertex;
            
            P.col(column) = (v->position - v2->position) * weights(v->index, v2->index);
            Pdash.col(column) = deformedCoords.row(v->index) - deformedCoords.row(v2->index);
            
            column ++;
            he = he->flip->next;
        
        } while (he != v->he);
        
        // compute covariance matrix
        Eigen::MatrixXd S = P * P.transpose();
        
        // compute svd
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(S, Eigen::ComputeThinU | Eigen::ComputeThinV);
        
        Eigen::MatrixXd U = svd.matrixU();
        Eigen::MatrixXd V = svd.matrixV();
        
        // compute determinant to determine the sign of the column of U corresponding to the smallest
        // singular value
        id(2, 2) = (V * U.transpose()).determinant();
        rotations[v->index] = V * id * U.transpose();
    }
}

void Mesh::deform(const int iterations)
{
    Eigen::MatrixXd b((int)vertices.size(), 3);
    b.setZero();
    
    for (int i = 0; i < iterations; i++) {
        // build b
        for (VertexCIter v = vertices.begin(); v != vertices.end(); v++) {
            
            Eigen::Vector3d p = v->position;
            
            if (!v->anchor && !v->handle) {
                p.setZero();
                HalfEdgeCIter he = v->he;
                do {
                    VertexCIter v2 = he->flip->vertex;
                    p += 0.5 * weights(v->index, v2->index) *
                        ((rotations[v->index] + rotations[v2->index]) * (v->position - v2->position));
                    
                    he = he->flip->next;
                    
                } while (he != v->he);
            }
            
            b.row(v->index) = p;
        }
    
        // update deformed positions
        for (int j = 0; j < 3; j++) {
            deformedCoords.col(j) = solver.solve(LT*b.col(j));
        }
        
        // compute rotations
        if (iterations > 1) computeRotations();
    }
    
    // update positions
    for (VertexIter v = vertices.begin(); v != vertices.end(); v++) {
        v->position = deformedCoords.row(v->index);
    }
}

void Mesh::normalize()
{
    // compute center of mass
    Eigen::Vector3d cm = Eigen::Vector3d::Zero();
    for (VertexCIter v = vertices.begin(); v != vertices.end(); v++) {
        cm += v->position;
    }
    cm /= (double)vertices.size();
    
    // translate to origin
    for (VertexIter v = vertices.begin(); v != vertices.end(); v++) {
        v->position -= cm;
    }
    
    // determine radius
    double rMax = 0;
    for (VertexCIter v = vertices.begin(); v != vertices.end(); v++) {
        rMax = std::max(rMax, v->position.norm());
    }
    
    // rescale to unit sphere
    for (VertexIter v = vertices.begin(); v != vertices.end(); v++) {
        v->position /= rMax;
    }
}
