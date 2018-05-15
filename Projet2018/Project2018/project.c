#include "fem.h"
#include <stdbool.h>


/*
 * AUTHORS
 * Ziegler Laurent
 * Lacroix Arthur
 */



bool checkTriangle(femMesh *theMesh, int elem, double xGrain, double yGrain){
 double x0 = theMesh->X[theMesh->elem[elem*theMesh->nLocalNode]];
 double x1 = theMesh->X[theMesh->elem[elem*theMesh->nLocalNode+1]];
 double x2 = theMesh->X[theMesh->elem[elem*theMesh->nLocalNode+2]];
 double y0 = theMesh->Y[theMesh->elem[elem*theMesh->nLocalNode]];
 double y1 = theMesh->Y[theMesh->elem[elem*theMesh->nLocalNode+1]];
 double y2 = theMesh->Y[theMesh->elem[elem*theMesh->nLocalNode+2]];

 bool s_ab = (x1-x0)*(yGrain-y0)-(y1-y0)*(xGrain-x0) > 0;

 if((x2-x0)*(yGrain-y0)-(y2-y0)*(xGrain-x0) > 0 == s_ab)
   return false;
 if((x2-x1)*(yGrain-y1)-(y2-y1)*(xGrain-x1) > 0 != s_ab)
   return false;
 return true;
}

int whereIsBrian(femMesh *theMesh, double xGrain, double yGrain){
  int i;
  for (i = 0; i < theMesh->nElem; i++){
    if(checkTriangle(theMesh, i, xGrain, yGrain)){
      return i;
    }
  }
}

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

double radOut(femPoissonProblem *theProblem){
  int i;
  femMesh *theMesh = theProblem->mesh;
  double rad = 0;
  for (i=0; i < theMesh->nNode; i++){
    double x = theMesh->X[i];
    double y = theMesh->Y[i];
    double dist = sqrt(pow(x,2)+pow(y,2));
    if (dist > rad){
      rad = dist;
    }
  }
  return rad;
}

double radIn(femPoissonProblem *theProblem){
  int i;
  femMesh *theMesh = theProblem->mesh;
  double x = theMesh->X[0];
  double y = theMesh->Y[0];
  double rad = sqrt(pow(x,2)+pow(y,2));
  for (i=1; i < theMesh->nNode; i++){
    x = theMesh->X[i];
    y = theMesh->Y[i];
    double dist = sqrt(pow(x,2)+pow(y,2));
    if (dist < rad){
      rad = dist;
    }
  }
  return rad;
}

double tau(femPoissonProblem *theProblem, int i, int iElem){
  femMesh *theMesh = theProblem->mesh;
  int nLocalNode = theMesh->nLocalNode;
  double phi[3];

  double x    = theMesh->X[i];
  double y    = theMesh->Y[i];
  double x1   = theMesh->X[theMesh->elem[iElem*nLocalNode+0]];
  double y1   = theMesh->Y[theMesh->elem[iElem*nLocalNode+0]];
  double x2   = theMesh->X[theMesh->elem[iElem*nLocalNode+1]];
  double y2   = theMesh->Y[theMesh->elem[iElem*nLocalNode+1]];
  double x3   = theMesh->X[theMesh->elem[iElem*nLocalNode+2]];
  double y3   = theMesh->Y[theMesh->elem[iElem*nLocalNode+2]];

  double xi   = (x1*(y3-y)-x3*(y1-y)-x*(y3-y1))/(x1*(y3-y2)+x2*(y1-y3)+x3*(y2-y1));
  double eta  = (x1*(y2-y)-x3*(y1-y)-x*(y2-y1))/(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2));
  femDiscretePhi2(theProblem->space, xi, eta, phi);
  double tau = x1*phi[0] + x2*phi[1] + x3*phi[2];
  return tau;
}


void femPoissonSolve(femPoissonProblem *theProblem, double omega, double mu, femGrains *myGrains)
{
  int elem, locNode, edge, i, j, k, map[4];
	double x[4], y[4], phi[4], dphidx[4], dphidy[4], dphidxi[4], dphideta[4];

  femMesh *theMesh = theProblem->mesh;
  int nLocalNode = theProblem->mesh->nLocalNode;
  int n = myGrains->n;

  double outerRad = radOut(theProblem);
  double critRad = .9*outerRad;
  double gamma       = myGrains->gamma;
  double *vy         = myGrains->vy;
  double *vx         = myGrains->vx;

  int grainTriangle[n];
  for (i = 0; i < n; i++){
    grainTriangle[i] = whereIsBrian(theProblem->mesh, myGrains->x[i], myGrains->y[i]);
  }

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
				dydxi	  += y[i]*dphidxi[i];
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
					theProblem->systemX->A[map[i]][map[j]] += weight * f * mu;
          //printf("first %f \n", theProblem->systemX->A[map[i]][map[j]]);
          theProblem->systemY->A[map[i]][map[j]] += weight * f * mu;
          for (k = 0; k < myGrains->n; k++){
            int elem = grainTriangle[k];
            //if ( elem*3+0 == i || elem*3+1 == i || elem*3+2 == i  ){
            //  if ( elem*3+0 == j || elem*3+1 == j || elem*3+2 == j){
                theProblem->systemX->A[map[i]][map[j]] += gamma * tau(theProblem, i, elem ) * tau(theProblem, j, elem );
                //printf("second %f \n", theProblem->systemX->A[map[i]][map[j]]);
                theProblem->systemY->A[map[i]][map[j]] += gamma * tau(theProblem, i, elem ) * tau(theProblem, j, elem );
            //  }
          //  }
          }
				}
        for (k = 0; k < myGrains->n; k++){
				  theProblem->systemX->B[map[i]] += gamma * vx[k] * tau(theProblem, i, elem )/*tauX(theMesh, i, nLocalNode,theProblem->rule->xsi, theProblem->rule->eta)*/;
          theProblem->systemY->B[map[i]] += gamma * vy[k] * tau(theProblem, i, elem )/*tauY(theMesh, i, nLocalNode,theProblem->rule->xsi, theProblem->rule->eta)*/;
        }
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
      if (theProblem->edges->edges[edge].elem[1] < 0 && dist > critRad && x >= 0 && y >= 0) {
			 	for (i = 0; i < 2; i++) {
			 		int node = theProblem->edges->edges[edge].node[i];
			 		femFullSystemConstrain(theProblem->systemX,node,+fabs(y/outerRad)*omega);
           femFullSystemConstrain(theProblem->systemY,node,-fabs(x/outerRad)*omega);
			 	}
			 }
       if (theProblem->edges->edges[edge].elem[1] < 0 && dist > critRad && x <= 0 && y >= 0) {
			 	for (i = 0; i < 2; i++) {
			 		int node = theProblem->edges->edges[edge].node[i];
			 		femFullSystemConstrain(theProblem->systemX,node,+fabs(y/outerRad)*omega);
           femFullSystemConstrain(theProblem->systemY,node,+fabs(x/outerRad)*omega);
			 	}
			 }
       if (theProblem->edges->edges[edge].elem[1] < 0 && dist > critRad && x <= 0 && y <= 0) {
			 	for (i = 0; i < 2; i++) {
			 		int node = theProblem->edges->edges[edge].node[i];
			 		femFullSystemConstrain(theProblem->systemX,node,-fabs(y/outerRad)*omega);
           femFullSystemConstrain(theProblem->systemY,node,+fabs(x/outerRad)*omega);
			 	}
			 }
       if (theProblem->edges->edges[edge].elem[1] < 0 && dist > critRad && x >= 0 && y <= 0) {
			 	for (i = 0; i < 2; i++) {
			 		int node = theProblem->edges->edges[edge].node[i];
			 		femFullSystemConstrain(theProblem->systemX,node,-fabs(y/outerRad)*omega);
           femFullSystemConstrain(theProblem->systemY,node,-fabs(x/outerRad)*omega);
			 	}
			 }
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


void femGrainsUpdate(femPoissonProblem *theProblem, femGrains *myGrains, double dt, double tol, double iterMax, double omega, double mu)
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

    int map;

    femFullSystem *theSystemX = theProblem->systemX;
    femFullSystem *theSystemY = theProblem->systemY;

    int grainTriangle[n];
    for (i = 0; i < n; i++){
      grainTriangle[i] = whereIsBrian(theProblem->mesh, myGrains->x[i], myGrains->y[i]);
    }

//
// -1- Calcul des nouvelles vitesses des grains sur base de la gravit� et de la trainee
//

    for(i = 0; i < n; i++) {
      int elem = grainTriangle[i];
      int locNode = theProblem->mesh->nLocalNode;
      int x0 = theProblem->mesh->elem[elem*locNode];
      int x1 = theProblem->mesh->elem[elem*locNode+1];
      int x2 = theProblem->mesh->elem[elem*locNode+2];
      int y0 = theProblem->mesh->elem[elem*locNode];
      int y1 = theProblem->mesh->elem[elem*locNode+1];
      int y2 = theProblem->mesh->elem[elem*locNode+2];
      double fx = m[i] * gx - gamma * (vx[i] - (theSystemX->B[x0] + theSystemX->B[x1] + theSystemX->B[x2])/3 );
      double fy = m[i] * gy - gamma * (vy[i] - (theSystemY->B[y0] + theSystemY->B[y1] + theSystemY->B[y2])/3 );
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
        y[i] += vy[i] * dt;
      }

}

# endif
