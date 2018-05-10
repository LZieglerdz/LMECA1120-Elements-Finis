#include "fem.h"
#include <stdbool.h>

# ifndef NOPOISSONCREATE

femPoissonProblem *femPoissonCreate(const char *filename)
{
    femPoissonProblem *theProblem = malloc(sizeof(femPoissonProblem));
    theProblem->mesh  = femMeshRead(filename);
    femMeshClean(theProblem->mesh);
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


double innerRadius(femPoissonProblem *theProblem){
  int i;
  double iRad,x,y,dist;

  x = theProblem->mesh->X[0];
  y = theProblem->mesh->Y[0];
  iRad = sqrt(pow(x,2)+pow(y,2));
  for(i = 1; i < theProblem->mesh->nNode; i++){
    x = theProblem->mesh->X[i];
    y = theProblem->mesh->Y[i];
    dist = sqrt(pow(x,2)+pow(y,2));
    if (dist < iRad)
      iRad = dist;
  }
  return iRad;
}
double outerRadius(femPoissonProblem *theProblem){
  int i;
  double iRad,x,y,dist;

  x = theProblem->mesh->X[0];
  y = theProblem->mesh->Y[0];
  iRad = sqrt(pow(x,2)+pow(y,2));
  for(i = 1; i < theProblem->mesh->nNode; i++){
    x = theProblem->mesh->X[i];
    y = theProblem->mesh->Y[i];
    dist = sqrt(pow(x,2)+pow(y,2));
    if (dist > iRad)
      iRad = dist;
  }
  return iRad;
}

bool whereIsBrian(femMesh *theMesh, double xGrain, double yGrain){
  int i;
  for(i = 0; i < theMesh->nElem; i++){
    double x0 = theMesh->X[theMesh->elem[i*theMesh->nLocalNode]];
    double x1 = theMesh->X[theMesh->elem[i*theMesh->nLocalNode+1]];
    double x2 = theMesh->X[theMesh->elem[i*theMesh->nLocalNode+2]];
    double y0 = theMesh->Y[theMesh->elem[i*theMesh->nLocalNode+0]];
    double y1 = theMesh->Y[theMesh->elem[i*theMesh->nLocalNode+1]];
    double y2 = theMesh->Y[theMesh->elem[i*theMesh->nLocalNode+2]];

    bool s_ab = (x1-x0)*(yGrain-y0)-(y1-y0)*(xGrain-x0) > 0;

    if((x2-x0)*(yGrain-y0)-(y2-y0)*(xGrain-x0) > 0 == s_ab)
      return false;
    if((x2-x1)*(yGrain-y1)-(y2-y1)*(xGrain-x1) > 0 != s_ab)
      return false;
    return true;
  }
}

double tau(femGrains *myGrains, femMesh *theMesh, int iElem, double U, double *phi){
  int i;
  double xGrain = myGrains->x[iElem];
  double yGrain = myGrains->y[iElem];
  double elem, tau;
  for (i = 0; i < theMesh->nElem; i++){
    if (whereIsBrian(theMesh, xGrain, yGrain) == true){
      elem = i;
      break;
    }
  }
  for (i = 0; i < theMesh->nLocalNode; i++){
    tau += phi[i] * U;
  }
}

void femPoissonSolve(femPoissonProblem *theProblem, femGrains *myGrains)
{
	int elem, locNode, edge, i,j, k, map[4];
	double x[4], y[4], phi[4], dphidx[4], dphidy[4], dphidxi[4], dphideta[4], mu, gamma;
  double *vx = myGrains->vx;
  double *vy = myGrains->vy;

  mu = 1;
  gamma = 0.5; // checker wiki pour d'autres coeffs

  for (elem = 0; elem < theProblem->mesh->nElem; elem++){
		femMeshLocal(theProblem->mesh, elem, map, x, y);
		for (locNode = 0; locNode < theProblem->rule->n; locNode++) {
			double xi, eta, weight, xLocal, yLocal, dxdxi, dydxi, dxdeta, dydeta;
			xi       = theProblem->rule->xsi[locNode];
			eta      = theProblem->rule->eta[locNode];
			weight   = theProblem->rule->weight[locNode];
			xLocal   = 0;
			yLocal   = 0;
			dxdxi    = 0;
			dydxi	   = 0;
			dxdeta   = 0;
			dydeta   = 0;
			femDiscretePhi2(theProblem->space, xi, eta, phi);
			femDiscreteDphi2(theProblem->space, xi, eta, dphidxi, dphideta);
			for (i = 0; i < theProblem->space->n; i++) {
				xLocal	+= x[i]*phi[i];
				yLocal	+= y[i]*phi[i];
				dxdxi 	+= x[i]*dphidxi[i];
				dydxi   += y[i]*dphidxi[i];
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
					theProblem->systemX->A[map[i]][map[j]] += weight * jacobian * mu * (dphidx[i]*dphidx[j] + dphidy[i]*dphidy[j]);
          theProblem->systemY->A[map[i]][map[j]] += weight * jacobian * mu * (dphidx[i]*dphidx[j] + dphidy[i]*dphidy[j]);

          for (k = 0; k < myGrains->n; i++){
            theProblem->systemX->A[map[i]][map[j]] += gamma * tau(myGrains, theProblem->mesh, k, xLocal, phi) * tau(myGrains, theProblem->mesh, k, xLocal, phi);
            theProblem->systemY->A[map[i]][map[j]] += gamma * tau(myGrains, theProblem->mesh, k, xLocal, phi) * tau(myGrains, theProblem->mesh, k, xLocal, phi);
          }
        }
        for (k = 0; k< myGrains->n; i++){
          theProblem->systemX->B[map[i]] += gamma * tau(myGrains, theProblem->mesh, k, xLocal, phi) * vx[k];
          theProblem->systemY->B[map[i]] += gamma * tau(myGrains, theProblem->mesh, k, xLocal, phi) * vy[k];
       }
      }
		}
	}
  // determiner si edge est au Rmax

  printf("%f \n", innerRadius(theProblem));
  printf("%f \n", outerRadius(theProblem));
  double mRad = outerRadius(theProblem) - innerRadius(theProblem);

/*
  for (i = 0; i<theProblem->system->size; i++){
    double node1X = theProblem->mesh->X[theProblem->edges->edges[i].node[0]];
    double node2X = theProblem->mesh->X[theProblem->edges->edges[i].node[1]];
    double node1Y = theProblem->mesh->Y[theProblem->edges->edges[i].node[0]];
    double node2Y = theProblem->mesh->Y[theProblem->edges->edges[i].node[1]];
    double medX, medY, medDist;

    if (node1X < node2X){
      medX = ((node1X - node2X)/2) + node1X;
    }
    if (node1X > node2X){
      medX = ((node1X - node2X)/2) + node2X;
    }
    if (node1Y < node2Y){
      medY = ((node1Y - node2Y)/2) + node1Y;
    }
    if (node1Y < node2Y){
      medY = ((node1Y - node2Y)/2) + node2Y;
    }

    medDist = sqrt( pow(medX,2) + pow(medY,2) );

    if ( theProblem->edges->edges[edge].elem[1] == -1 && medDist < mRad) //accès noeud, pas encore coord noeuds ><
      theProblem->edges->edges[i].elem[1] = -2;
  }*/

  for (edge = 0; edge < theProblem->edges->nEdge; edge++) {
    if (theProblem->edges->edges[edge].elem[1] < 0  ) { //&& sqrt(pow(theProblem->edges->edges[edge].node[0],2) + pow(theProblem->edges->edges[edge].node[1],2)) > mRad
      printf("%d %d \n", edge, theProblem->edges->edges[edge].elem[1]);
			for (i = 0; i < 2; i++) {
				int node = theProblem->edges->edges[edge].node[i];
				femFullSystemConstrain(theProblem->systemX,node,0);
        femFullSystemConstrain(theProblem->systemY,node,0);
			}
		}
    /*if (theProblem->edges->edges[edge].elem[1] == -2  ) { //&& sqrt(pow(theProblem->edges->edges[edge].node[0],2) + pow(theProblem->edges->edges[edge].node[1],2)) > mRad
      printf("%d %d \n", edge, theProblem->edges->edges[edge].elem[1]);
			for (i = 0; i < 2; i++) {
				int node = theProblem->edges->edges[edge].node[i];
				femFullSystemConstrain(theProblem->system,node,-1);
			}
		}*/
	}
 femFullSystemEliminate(theProblem->systemX);
 femFullSystemEliminate(theProblem->systemY);

}

# endif


# ifndef NOCONTACTITERATE

double femGrainsContactIterate(femGrains *myGrains, double dt, int iter)
{
    int i,j,iContact;
    int n = myGrains->n;


    double *x          = myGrains->x;
    double *y          = myGrains->y;
    double *m          = myGrains->m;
    double *r          = myGrains->r;
    double *vy         = myGrains->vy;
    double *vx         = myGrains->vx;
    double *dvBoundary = myGrains->dvBoundary;
    double *dvContacts = myGrains->dvContacts;
    double rIn         = myGrains->radiusIn;
    double rOut        = myGrains->radiusOut;

    double error = 0.0;
    double gap,rr, rx, ry, nx, ny, vn, dv, dvIn, dvOut;

    error = 0;
    iContact = 0;
    for(i = 0; i < n; i++) {
        for(j = i+1; j < n; j++) {
            rx = (x[j]-x[i]);
            ry = (y[j]-y[i]);
            rr = sqrt(rx*rx+ry*ry);
            nx = rx/rr;
            ny = ry/rr;
            if (iter == 0) {
                  dv = dvContacts[iContact]; }
            else {
                  vn = (vx[i]-vx[j])*nx + (vy[i]-vy[j])*ny ;
                  gap = rr - (r[i]+r[j]);
                  dv = fmax(0.0, vn + dvContacts[iContact] - gap/dt);
                  dv = dv - dvContacts[iContact];
                  dvContacts[iContact] += dv;
                  error = fmax(fabs(dv),error); }
            vx[i] -= dv * nx * m[j] / ( m[i] + m[j] );
            vy[i] -= dv * ny * m[j] / ( m[i] + m[j] );
            vx[j] += dv * nx * m[i] / ( m[i] + m[j] );
            vy[j] += dv * ny * m[i] / ( m[i] + m[j] );
            iContact++; }}

    for(i = 0; i < n; i++) {
        rr = sqrt(x[i]*x[i]+y[i]*y[i]);
        nx = x[i]/rr;
        ny = y[i]/rr;
        if (iter == 0) {
            dv = dvBoundary[i]; }
        else {
            vn = vx[i]*nx + vy[i]*ny ;
            gap = rOut - rr - r[i];
            dvOut = fmax(0.0, vn + dvBoundary[i] - gap/dt);
            gap = rr - rIn - r[i];
            dvIn  = fmax(0.0,-vn - dvBoundary[i] - gap/dt);
            dv = dvOut - dvIn - dvBoundary[i];
            dvBoundary[i] += dv;
            error = fmax(fabs(dv),error); }
        vx[i] -= dv * nx;
        vy[i] -= dv * ny; }
    return error;
}

# endif
# ifndef NOUPDATE


void femGrainsUpdate(femGrains *myGrains, double dt, double tol, double iterMax)
{
    int i;
    int n = myGrains->n;

    double *x          = myGrains->x;
    double *y          = myGrains->y;
    double *m          = myGrains->m;
    double *vy         = myGrains->vy;
    double *vx         = myGrains->vx;
    double gamma       = myGrains->gamma;
    double gx          = myGrains->gravity[0];
    double gy          = myGrains->gravity[1];

//
// -1- Calcul des nouvelles vitesses des grains sur base de la gravit� et de la trainee
//

    for(i = 0; i < n; i++) {
      double fx = m[i] * gx - gamma * (vx[i] );
      double fy = m[i] * gy - gamma * (vy[i] );
      vx[i] += fx * dt / m[i];
      vy[i] += fy * dt / m[i];
    }

//
// -2- Correction des vitesses pour tenir compte des contacts
//

    int iter = 0;
    double error;

    do {
        error = femGrainsContactIterate(myGrains,dt,iter);
        iter++; }
    while ((error > tol/dt && iter < iterMax) || iter == 1);
    printf("iterations = %4d : error = %14.7e \n",iter-1,error);

//
// -3- Calcul des nouvelles positions sans penetrations de points entre eux
//

    for (i = 0; i < n; ++i) {
        x[i] += vx[i] * dt;
        y[i] += vy[i] * dt; }
}

# endif
