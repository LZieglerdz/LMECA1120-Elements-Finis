#include"fem.h"

# ifndef NOCOMPARE

int edgesCompare(const void *edgeOne, const void *edgeTwo)
{
    int *nodeOne = ((femEdge*) edgeOne)->node;
    int *nodeTwo = ((femEdge*) edgeTwo)->node;  
    int  diffMin = fmin(nodeOne[0],nodeOne[1]) - fmin(nodeTwo[0],nodeTwo[1]);
    int  diffMax = fmax(nodeOne[0],nodeOne[1]) - fmax(nodeTwo[0],nodeTwo[1]);
    
    if (diffMin < 0)    return  1;
    if (diffMin > 0)    return -1;
    if (diffMax < 0)    return  1;
    if (diffMax > 0)    return -1; 
    return  0;
}


# ifndef NOEXPAND
# endif

void edgesExpand(femEdges *theEdges)
{
    femMesh *theMesh = theEdges->mesh;
    femEdge *edges = theEdges->edges;
    int i,j,nLoc = theMesh->nLocalNode;
 
    for (i = 0; i < theMesh->nElem; i++) {
        int *elem = &(theMesh->elem[i*nLoc]);
        for (j = 0; j < nLoc; j++) {
            int id = i * nLoc + j;
            edges[id].elem[0] = i;
            edges[id].elem[1] = -1;
            edges[id].node[0] = elem[j];
            edges[id].node[1] = elem[(j + 1) % nLoc]; }}
}

# endif
# ifndef NOSORT

void edgesSort(femEdges *theEdges)
{
    qsort(theEdges->edges, theEdges->nEdge, sizeof(femEdge), edgesCompare);
}

# endif
# ifndef NOSHRINK
    
void edgesShrink(femEdges *theEdges)
{
    int i,index = 0;          
    int nBoundary = 0;
    
    femEdge *edges = theEdges->edges;
    for (i=0; i < theEdges->nEdge; i++) {
      if (i == theEdges->nEdge - 1 || edgesCompare(&edges[i],&edges[i+1]) != 0) {
              edges[index] = edges[i];
              nBoundary++; }
      else {  edges[index] = edges[i];
              edges[index].elem[1] = edges[i+1].elem[0];
              i = i+1;}
      index++; }
      
    theEdges->edges = realloc(edges, index * sizeof(femEdge));
    theEdges->nEdge = index;
    theEdges->nBoundary = nBoundary;
}


# endif
# ifndef NOBOUNDARYLENGTH

double edgesBoundaryLength(femEdges *theEdges)
{
    double L = 0;
    double x[2],y[2],dx,dy;

    int i,j;
    
    for (i=0; i < theEdges->nEdge; ++i) {
      if (theEdges->edges[i].elem[1] == -1) {
        for (j=0; j < 2; ++j) {
            x[j] = theEdges->mesh->X[theEdges->edges[i].node[j]];
            y[j] = theEdges->mesh->Y[theEdges->edges[i].node[j]];
              }
        dx = x[1]-x[0]; dy = y[1]-y[0];
        L += sqrt(dx*dx+dy*dy); }}
    return L;
}

# endif
