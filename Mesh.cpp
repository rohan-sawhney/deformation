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
        didSetup = false;
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
                double coefficient = (he->cotan() + he->flip->cotan()) / 2.0;
                coefficient = std::max(0.0, coefficient);
                weights[he->edge->index] = coefficient;
                sumCoefficients += coefficient;
                
                LTriplet.push_back(Eigen::Triplet<double>(v->index, he->flip->vertex->index, -coefficient));
                
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
    
    weights.clear();
    
    Eigen::SparseMatrix<double> L(v, v);
    buildLaplacian(L);
    
    LT = L.transpose();
    solver.compute(LT*L);
    
    deformedCoords.resize(v, 3);
    deformedCoords.setZero();
    
    rotations.clear();
    rotations.resize(v, Eigen::Matrix3d::Identity());
    
    didSetup = true;
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
            
            P.col(column) = (v->reference - v2->reference) * weights[he->edge->index];
            Pdash.col(column) = deformedCoords.row(v->index) - deformedCoords.row(v2->index);
            
            column ++;
            he = he->flip->next;
        
        } while (he != v->he);
        
        // compute covariance matrix
        Eigen::MatrixXd S = P * Pdash.transpose();
        
        // compute svd
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(S, Eigen::ComputeThinU | Eigen::ComputeThinV);
        
        Eigen::MatrixXd UT = svd.matrixU().transpose();
        Eigen::MatrixXd V = svd.matrixV();
        
        // compute determinant to determine the sign of the column of U corresponding to the smallest
        // singular value
        id(2, 2) = (V * UT).determinant();
        rotations[v->index] = V * id * UT;
    }
}

void Mesh::deform(const int iterations)
{
    if (!didSetup) setup();
    
    Eigen::MatrixXd b((int)vertices.size(), 3);    
    for (int i = 0; i < iterations; i++) {
        // build b
        for (VertexCIter v = vertices.begin(); v != vertices.end(); v++) {
            
            if (!v->anchor && !v->handle) {
                
                Eigen::Vector3d p = Eigen::Vector3d::Zero();
                HalfEdgeCIter he = v->he;
                do {
                    VertexCIter v2 = he->flip->vertex;
                    p += ((rotations[v->index] + rotations[v2->index]) * (v->reference - v2->reference) *
                          weights[he->edge->index] / 2.0);
                    
                    he = he->flip->next;
                    
                } while (he != v->he);
            
                b.row(v->index) = p;
                
            } else {
                b.row(v->index) = v->position;
            }
        }
        
        // update deformed positions
        b = LT*b;
        for (int j = 0; j < 3; j++) {
            deformedCoords.col(j) = solver.solve(b.col(j));
        }
        
        // compute rotations after initial naive laplacian guess
        if (i > 0) computeRotations();
    }

    // update positions
    for (VertexIter v = vertices.begin(); v != vertices.end(); v++) {
        v->position = deformedCoords.row(v->index).transpose();
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
