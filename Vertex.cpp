#include "Vertex.h"
#include "HalfEdge.h"
#include "Face.h"

std::vector<HalfEdge> isolated;

bool Vertex::isIsolated() const
{
    return he == isolated.begin();
}

int Vertex::valence() const
{
    int v = 0;
    HalfEdgeCIter h = he;
    do {
        v ++;
        h = h->flip->next;
        
    } while (h != he);
    
    return v;
}

double Vertex::dualArea() const
{
    double area = 0.0;
    
    HalfEdgeCIter h = he;
    do {
        area += h->face->area();
        h = h->flip->next;
        
    } while (h != he);
    
    return area / 3.0;
}