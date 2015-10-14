#ifndef FACE_H
#define FACE_H

#include "Types.h"

class Face {
public:
    // one of the halfedges associated with this face
    HalfEdgeIter he;
    
    // id between 0 and |F|-1
    int index;
    
    // checks if this face lies on boundary
    bool isBoundary() const;
    
    // returns face area
    double area() const;
    
    // returns normal to face
    Eigen::Vector3d normal() const;
};

#endif
