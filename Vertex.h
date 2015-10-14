#ifndef VERTEX_H
#define VERTEX_H

#include "Types.h"

class Vertex {
public:
    // outgoing halfedge
    HalfEdgeIter he;
    
    // location in 3d
    Eigen::Vector3d position;
    
    // id between 0 and |V|-1
    int index;
    
    // flag for anchor vertex
    bool anchor;
    
    // flag for control vertex
    bool handle;
    
    // returns valence
    int valence() const;
    
    // checks if vertex is contained in any edge or face
    bool isIsolated() const;
};

#endif
