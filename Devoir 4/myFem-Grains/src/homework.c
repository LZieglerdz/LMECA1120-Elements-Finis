/*
 * Auteurs: ZIEGLER Laurent & LACROIX Arthur
 *
 */


#include"fem.h"
#include"math.h"


#ifndef NOCONTACTITERATE

double femGrainsContactIterate(femGrains *myGrains, double dt, int iter)  
{
	int i, j, k = 0; 
	int n 			= myGrains->n; 			//# particules
	double *x		= myGrains->x;			//abscisse particule
	double *y		= myGrains->y;			//ordonnée particule
	double *m		= myGrains->m;			//masse particule
	double *r		= myGrains->r;			//rayon particule
	double *vy		= myGrains->vy;			//vitesse en ordonnée
	double *vx		= myGrains->vx;			//vitesse en abscisse
	double *dvBoundary	= myGrains->dvBoundary;		//en français dans le texte
	double *dvContacts	= myGrains->dvContacts;		//idem
	double rIn		= myGrains->radiusIn;		//rayon de la frontière interne (?)
	double rOut		= myGrains->radiusOut;		//rayon de la frontière externe (logiquement...)
    	double *g		= myGrains->gravity;		//vecteur des accélérations subies par les particule ( g[0]->ax; g[1]->ay )

	double zeta = 0.0;	
	
	double dOut, dIn, norm;
	double gamma = 0;
	double vn = 0;
	double dv = 0;
	double nrml[2] = {0,0};
	double nContacts = n*(n-1)/2;


	//remise à zero des arrays dvBoundary et dvContact après la dernière itération
//	if (iter == 0) {
		for (i = 0; i < nContacts; i++) {
			if (i < n) {
				dvBoundary[i] = 0;
			}
			dvContacts[i] = 0;
		}
//	}
	
	for (i = 0; i < n-1; i++) {
		for (j = i + 1; j < n; j++) {
			norm = sqrt( pow(x[i]-x[j],2) + pow(y[i]-y[j],2) );
			gamma = norm - (r[i] + r[j]);
			nrml[0] = (x[j] - x[i])/norm;
			nrml[1] = (y[j] - y[i])/norm;
			vn = (vx[i] - vx[j])*nrml[0] + (vy[i] - vy[j])*nrml[1];
			dv = fmax(0, vn + dvContacts[k] - gamma/dt) - dvContacts[k];
			vx[i] -= dv*nrml[0]*m[j]/(m[i]+m[j]);
			vy[i] -= dv*nrml[1]*m[j]/(m[i]+m[j]);
			vx[j] += dv*nrml[0]*m[i]/(m[i]+m[j]);
			vy[j] += dv*nrml[1]*m[i]/(m[i]+m[j]);
			dvContacts[k] += dv;

			zeta = fmax(zeta, fabs(dv));
			k++;	
		}
	}	
	
	for (i = 0; i < n; i++) {
		norm 		 = sqrt( pow(x[i],2) + pow(y[i],2) );
		nrml[0] 	 = x[i]/norm;
		nrml[1] 	 = y[i]/norm;
		vn 		 = vx[i]*nrml[0] + vy[i]*nrml[1];
		dOut 		 = fmax(0, vn + dvBoundary[i] - (rOut - norm - r[i])/dt);
		dIn 		 = fmax(0, -vn - dvBoundary[i] - (norm  -rIn - r[i])/dt);
		dv		 = dOut - dIn - dvBoundary[i];
		vx[i]		-= dv*nrml[0];
		vy[i]		-= dv*nrml[1];
		dvBoundary[i]	+= dv;
		zeta 		 = fmax(zeta, fabs(dv));

	}

    return zeta;

}

#endif
#ifndef NOUPDATE

void femGrainsUpdate(femGrains *myGrains, double dt, double tol, double iterMax)
{
    int n = myGrains->n;
    int i,iter = 0;
    double zeta;
    double *x          = myGrains->x;
    double *y          = myGrains->y;
    double *m          = myGrains->m;
    double *vy         = myGrains->vy;
    double *vx         = myGrains->vx;
    double gamma       = myGrains->gamma;
    double gx          = myGrains->gravity[0];
    double gy          = myGrains->gravity[1];

 
// 
// -1- Calcul des nouvelles vitesses des grains sur base de la gravité et de la trainee
//
//
	for (i=0; i<n; i++) {
	    	vy[i] += gy*dt - gamma*vy[i]*dt/m[i];
		vx[i] += gx*dt - gamma*vx[i]*dt/m[i];
	}

//
// -2- Correction des vitesses pour tenir compte des contacts        
//       
    do {
        zeta = femGrainsContactIterate(myGrains,dt,iter);
        iter++; }
    while ((zeta > tol/dt && iter < iterMax) || iter == 1);
    printf("iterations = %4d : error = %14.7e \n",iter-1,zeta);
 
//  
// -3- Calcul des nouvelles positions sans penetrations de points entre eux
//
    for (i = 0; i < n; ++i) {
        x[i] += vx[i] * dt;
        y[i] += vy[i] * dt; }
}


#endif
