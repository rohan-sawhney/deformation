#include "Vertex.h"
#include "HalfEdge.h"

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