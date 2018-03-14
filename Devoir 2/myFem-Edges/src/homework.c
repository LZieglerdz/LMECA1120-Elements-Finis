#include "fem.h"


# ifndef NOEXPAND

void edgesExpand(femEdges *theEdges){	

/* //debug print
	printf("elem %d\n", theEdges->mesh->elem[0]);			//
	printf("coord node in X %f\n", theEdges->mesh->X[0]);		//coord X
	printf("coord node in Y %f\n", theEdges->mesh->Y[0]);		//coord Y
	printf("number of nodes %d\n", theEdges->mesh->nNode);		//# of nodes
	printf("number of triangles %d\n", theEdges->mesh->nElem); 	//# of triangles
	
	*/
	
	int i;
	for (i = 0; i < ((theEdges->nEdge)/3); ++i) {
//	set first node
		theEdges->edges[3*i].node[0] = theEdges->mesh->elem[3*i];
		theEdges->edges[3*i+1].node[0] = theEdges->mesh->elem[3*i+1] ;
		theEdges->edges[3*i+2].node[0] = theEdges->mesh->elem[3*i+2];
//	set second node
	       	theEdges->edges[3*i].node[1] = theEdges->mesh->elem[3*i+1];		
		theEdges->edges[3*i+1].node[1] = theEdges->mesh->elem[3*i+2] ;			     
		theEdges->edges[3*i+2].node[1] = theEdges->mesh->elem[3*i];
//	set triangle number
		theEdges->edges[3*i].elem[0] = i;
		theEdges->edges[3*i+1].elem[0] = i;
		theEdges->edges[3*i+2].elem[0] = i;
//	fill last column
		theEdges->edges[3*i].elem[1] = -1;
		theEdges->edges[3*i+1].elem[1] = -1;
		theEdges->edges[3*i+2].elem[1] = -1;

	}

}
# endif
# ifndef NOSORT

void edgesSort(femEdges *theEdges)
{
	qsort(theEdges->edges, theEdges->nEdge, sizeof(femEdge), edgesCompare);
}

# endif
# ifndef NOCOMPARE

int edgesCompare(const void* e0, const void *e1)
{   
 	int diagnostic = ((femEdge*) e0)->node[0] - ((femEdge*) e1)->node[0];
    	int diagnosticBis = ((femEdge*) e0)->node[1] - ((femEdge*) e1)->node[1];
	if (diagnostic > 0) return -1;
	if (diagnostic < 0) return 1;
	if (diagnostic == 0) {
		if (diagnosticBis > 0) return -1;
		if (diagnosticBis < 0) return 1;
	}
	else return 0;
}

# endif
# ifndef NOSHRINK



void edgesShrink(femEdges *theEdges)
{
	int l = theEdges->nEdge;
//	printf("%d\n",l);
	int i,j,k;
	int doublon = 0;	//ne pas oublier de compter ces PUTAINS DE DOUBLONS
	int n = 0;          // Nouveau nombre total de segments :(sizeof(tab)/sizeoftab[0] )
	int nBoundary = 0;  // Nombre de segments frontieres : A MODIFIER


//	int test=0;
//	int t;

	//boucle pour fixer indices doubles
	for (i = 0; i < l; i++){
		for (j = 0; j < l; j++){
			if (theEdges->edges[i].node[0] == theEdges->edges[j].node[1]	&&	theEdges->edges[j].node[0] == theEdges->edges[i].node[1]) {
				theEdges->edges[i].elem[1] = theEdges->edges[j].elem[0];

				for (k = j; k < l -1 ; k++) {
					theEdges->edges[k] = theEdges->edges[k+1];
				}
				doublon++;
	
			}
	
		}
	
	}

	n = l - doublon;
//	printf("%d\n", n);
	for (i = 0; i < n; i++) {
		if (theEdges->edges[i].elem[1] == -1) {
			nBoundary++;
		}
	}
    
    // Reallocation du tableau des edges
    
    theEdges->edges = realloc(theEdges->edges, n * sizeof(femEdge));
    theEdges->nEdge = n;
    theEdges->nBoundary = nBoundary;


}

# endif
# ifndef NOBOUNDARYLENGTH

double edgesBoundaryLength(femEdges *theEdges)
{
	double L = 0;
	int i;
	for (i = 0; i< theEdges->nEdge; i++){
		if (theEdges->edges[i].elem[1] == -1){
			double nodeX1, nodeX2, nodeY1, nodeY2;	
			nodeX1 = theEdges->mesh->X[theEdges->edges[i].node[0]];
			nodeX2 = theEdges->mesh->X[theEdges->edges[i].node[1]];
			nodeY1 = theEdges->mesh->Y[theEdges->edges[i].node[0]];
			nodeY2 = theEdges->mesh->Y[theEdges->edges[i].node[1]];
			
			L += sqrt( pow((nodeX1-nodeX2),2) + pow((nodeY1-nodeY2),2) ); 
		}
	}
	return L;
}

# endif
