/*
 * AUTHORS
 * Ziegler Laurent
 * Lacroix Arthur
 */

#include"fem.h"

# ifndef NOPOISSONCREATE

femPoissonProblem *femPoissonCreate(const char *filename)
{
    femPoissonProblem *theProblem = malloc(sizeof(femPoissonProblem));
    theProblem->mesh  = femMeshRead(filename);           
    theProblem->edges = femEdgesCreate(theProblem->mesh);  
    if (theProblem->mesh->nLocalNode == 4) {
        theProblem->space = femDiscreteCreate(4,FEM_QUAD);
        theProblem->rule = femIntegrationCreate(4,FEM_QUAD); }
    else if (theProblem->mesh->nLocalNode == 3) {
        theProblem->space = femDiscreteCreate(3,FEM_TRIANGLE);
        theProblem->rule = femIntegrationCreate(3,FEM_TRIANGLE); }
    theProblem->system = femFullSystemCreate(theProblem->mesh->nNode);
    return theProblem;
}

# endif
# ifndef NOPOISSONFREE

void femPoissonFree(femPoissonProblem *theProblem)
{
    femFullSystemFree(theProblem->system);
    femIntegrationFree(theProblem->rule);
    femDiscreteFree(theProblem->space);
    femEdgesFree(theProblem->edges);
    femMeshFree(theProblem->mesh);
    free(theProblem);
}
    
# endif
# ifndef NOMESHLOCAL

void femMeshLocal(const femMesh *theMesh, const int i, int *map, double *x, double *y)
{
	int j;
	for (j = 0; j < theMesh->nLocalNode; j++) {
	 	map[j] = theMesh->elem[i*theMesh->nLocalNode + j];
		x[j] = theMesh->X[map[j]];
		y[j] = theMesh->Y[map[j]];
	}
}

# endif
# ifndef NOPOISSONSOLVE


void femPoissonSolve(femPoissonProblem *theProblem)
{
	int elem, locNode, edge, i,j, map[4];
	double x[4], y[4], phi[4], dphidx[4], dphidy[4], dphidxi[4], dphideta[4];
	for (elem = 0; elem < theProblem->mesh->nElem; elem++){
		femMeshLocal(theProblem->mesh, elem, map, x, y);
		for (locNode = 0; locNode < theProblem->rule->n; locNode++) {
			double xi, eta, weight, xLocal, yLocal, dxdxi, dydxi, dxdeta, dydeta;
			xi	= theProblem->rule->xsi[locNode];
			eta	= theProblem->rule->eta[locNode];
			weight	= theProblem->rule->weight[locNode];
			xLocal	= 0;
			yLocal	= 0;
			dxdxi	= 0;
			dydxi	= 0;
			dxdeta	= 0;
			dydeta	= 0;
			femDiscretePhi2(theProblem->space, xi, eta, phi);
			femDiscreteDphi2(theProblem->space, xi, eta, dphidxi, dphideta); 
			for (i = 0; i < theProblem->space->n; i++) {
				xLocal	+= x[i]*phi[i];
				yLocal	+= y[i]*phi[i];
				dxdxi 	+= x[i]*dphidxi[i];
				dydxi	+= y[i]*dphidxi[i];
				dxdeta 	+= x[i]*dphideta[i];
				dydeta 	+= y[i]*dphideta[i];
			}	
			double jacobian = fabs(dxdxi * dydeta - dydxi * dxdeta);
			for (i = 0; i < theProblem->space->n; i++) {							//calcul des dphidx et dphidy
				dphidx[i] = (1/jacobian) * (dphidxi[i]*dydeta - dphideta[i]*dydxi);
				dphidy[i] = (1/jacobian) * (dphidxi[i]*dxdeta - dphideta[i]*dxdxi);
			}
			for (i = 0; i < theProblem->space->n; i++) {							//remplissage de la matrice A par integration 
				for (j = 0; j < theProblem->space->n; j++) {
					double f = (dphidx[i]*dphidx[j] + dphidy[i]*dphidy[j]) * jacobian;
					theProblem->system->A[map[i]][map[j]] += weight * f;
				}
				theProblem->system->B[map[i]] += weight * phi[i] * jacobian;
			}
		}
	}
		for (edge = 0; edge < theProblem->edges->nEdge; edge++) {
			if (theProblem->edges->edges[edge].elem[1] < 0) {
				for (i = 0; i < 2; i++) {
					int node = theProblem->edges->edges[edge].node[i];
					femFullSystemConstrain(theProblem->system,node,0);
				}
			}
		}
 	femFullSystemEliminate(theProblem->system);

}

# endif
