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
        
        HalfEdgeCIter he = v->he;
        double sumCoefficients = 0.0;
        if (!v->anchor && !v->handle) {
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
    
    solver.compute(L);
}

void Mesh::computeRotationsSVD()
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
        Eigen::Matrix3d S = P * P.transpose();
        
        // compute svd
        Eigen::JacobiSVD<Eigen::Matrix3d> svd(S, Eigen::ComputeThinU | Eigen::ComputeThinV);
        
        Eigen::Matrix3d U = svd.matrixU();
        Eigen::Matrix3d V = svd.matrixV();
        
        // compute determinant to determine the sign of the column of U corresponding to the smallest
        // singular value
        id(2, 2) = (V * U.transpose()).determinant();
        rotations[v->index] = V * id * U.transpose();
    }
}

void Mesh::deform(const int iterations)
{
    for (int i = 0; i < iterations; i++) {
        
        // build b
        Eigen::VectorXd b((int)vertices.size(), 3);
        for (VertexCIter v = vertices.begin(); v != vertices.end(); v++) {
            
            if (!v->anchor && !v->handle) {
                Eigen::Vector3d p = Eigen::Vector3d::Zero();
                HalfEdgeCIter he = v->he;
                do {
                    VertexCIter v2 = he->flip->vertex;
                    
                    p += 0.5 * weights(v->index, v2->index) *
                        ((rotations[v->index] + rotations[v2->index]) * (v->position - v2->position));
                    
                    he = he->flip->next;
                    
                } while (he != v->he);
                
                b.row(v->index)[0] = p.x();
                b.row(v->index)[1] = p.y();
                b.row(v->index)[2] = p.z();
            }
        }
        
        // update deformed positions
        for (int j = 0; j < 3; j++) {
            deformedCoords.col(j) = solver.solve(b.col(j));
        }
        
        // compute rotations
        if (iterations > 0) computeRotationsSVD();
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
