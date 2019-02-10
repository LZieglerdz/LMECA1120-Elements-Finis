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
    theProblem->systemX = femFullSystemCreate(theProblem->mesh->nNode);
    theProblem->systemY = femFullSystemCreate(theProblem->mesh->nNode);
    return theProblem;
}

# endif
# ifndef NOPOISSONFREE

void femPoissonFree(femPoissonProblem *theProblem)
{
    femFullSystemFree(theProblem->systemX);
    femFullSystemFree(theProblem->systemY);
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

double radius(femPoissonProblem *theProblem){
  int i;
  femMesh *theMesh = theProblem->mesh;
  double rad = 0;
  for (i=0; i < theMesh->nNode; i++){
    double x = theMesh->X[i];
    double y = theMesh->Y[i];
    printf("%d %f %f\n", i, x, y );
    double dist = sqrt(pow(x,2)+pow(y,2));
    if (dist > rad){
      rad = dist;
    }
  }
  return 0.9*rad;
}




void femPoissonSolve(femPoissonProblem *theProblem)
{

  double omega = 1;   //En rad/s

  int elem, locNode, edge, i,j, map[4];
	double x[4], y[4], phi[4], dphidx[4], dphidy[4], dphidxi[4], dphideta[4];

  femMesh *theMesh = theProblem->mesh;

  double critRad = radius(theProblem);
  printf( "%f \n", critRad);


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
					theProblem->systemX->A[map[i]][map[j]] += weight * f;
          theProblem->systemY->A[map[i]][map[j]] += weight * f;
				}
				theProblem->systemX->B[map[i]] += weight * phi[i] * jacobian;
        theProblem->systemY->B[map[i]] += weight * phi[i] * jacobian;
			}
		}
	}
		for (edge = 0; edge < theProblem->edges->nEdge; edge++) {
      double x = theMesh->X[theProblem->edges->edges[edge].node[0]];
      double y = theMesh->Y[theProblem->edges->edges[edge].node[1]];
      double dist = sqrt(pow(x,2)+pow(y,2));

      if (theProblem->edges->edges[edge].elem[1] < 0 && dist < critRad ) {
				for (i = 0; i < 2; i++) {
					int node = theProblem->edges->edges[edge].node[i];
					femFullSystemConstrain(theProblem->systemX,node,0);
          femFullSystemConstrain(theProblem->systemY,node,0);
				}
			}
      if (theProblem->edges->edges[edge].elem[1] < 0 && dist > critRad && x > 0 && y > 0) {
				for (i = 0; i < 2; i++) {
					int node = theProblem->edges->edges[edge].node[i];
					femFullSystemConstrain(theProblem->systemX,node,-omega);
          femFullSystemConstrain(theProblem->systemY,node,omega);
				}
			}
      if (theProblem->edges->edges[edge].elem[1] < 0 && dist > critRad && x < 0 && y > 0) {
				for (i = 0; i < 2; i++) {
					int node = theProblem->edges->edges[edge].node[i];
					femFullSystemConstrain(theProblem->systemX,node,-omega);
          femFullSystemConstrain(theProblem->systemY,node,-omega);
				}
			}
      if (theProblem->edges->edges[edge].elem[1] < 0 && dist > critRad && x < 0 && y < 0) {
				for (i = 0; i < 2; i++) {
					int node = theProblem->edges->edges[edge].node[i];
					femFullSystemConstrain(theProblem->systemX,node,omega);
          femFullSystemConstrain(theProblem->systemY,node,-omega);
				}
			}
      if (theProblem->edges->edges[edge].elem[1] < 0 && dist > critRad && x > 0 && y < 0) {
				for (i = 0; i < 2; i++) {
					int node = theProblem->edges->edges[edge].node[i];
					femFullSystemConstrain(theProblem->systemX,node,omega);
          femFullSystemConstrain(theProblem->systemY,node,omega);
				}
			}
		}
 	femFullSystemEliminate(theProblem->systemX);
  femFullSystemEliminate(theProblem->systemY);

}

# endif
