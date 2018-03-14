nclude"fem.h"
#include"math.h"

# ifndef NOCONTACTITERATE

double femGrainsContactIterate(femGrains *myGrains, double dt, int iter)  {
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

	if (iter == 0)	{
		int i;
		for (i = 0;i < n*(n - 1)/2;i++)
		{
			if (i < n)
			{
				dvBoundary[i] = 0;
			}
			dvContacts[i] = 0;
		}
	}
	double dist,dv,dvout,dvin,norme;
	double vn = 0;
	double normal[2] = {0,0};
	int i;
	int k = 0;
	for (i = 0; i < n - 1; i++)
	{
		int j;
		for (j = i + 1; j < n; j++)
		{
			norme = sqrt(pow(x[i] - x[j], 2) + pow(y[i] - y[j], 2));
			normal[0] = (x[j] - x[i])/norme;
			normal[1] = (y[j] - y[i])/norme;
			vn = (vx[i] - vx[j])*normal[0]+ (vy[i] - vy[j])*normal[1];
			dist = norme - (r[i] + r[j]);
			dv = fmax(0, vn + dvContacts[k]-dist/dt);
			dv = dv - dvContacts[k];
			vx[i] = vx[i] - (m[j] / (m[i] + m[j]))*dv*normal[0];
			vy[i] = vy[i] - (m[j] / (m[i] + m[j]))*dv*normal[1];
			vx[j] = vx[j] + (m[i] / (m[i] + m[j]))*dv*normal[0];
			vy[j] = vy[j] + (m[i] / (m[i] + m[j]))*dv*normal[1];
			error = fmax(fabs(dv), error);
			dvContacts[k] = dvContacts[k] + dv;
			k++;
		}

	}


	
	for(i=0; i < n ; i++)

	{

		norme = sqrt(pow(x[i], 2) + pow(y[i], 2));

		normal[0] =  x[i] / norme;

		normal[1] =  y[i] / norme;

		vn = vx[i] * normal[0] + vy[i] * normal[1];

		dvout = fmax(0, vn + dvBoundary[i] - (rOut - norme - r[i]) / dt);

		dvin = fmax(0,  -vn - dvBoundary[i] - (norme - rIn - r[i]) / dt);

		dv=dvout-dvin;

		dv = dv - dvBoundary[i];



		vx[i] = vx[i] - 1 * dv*normal[0];

		vy[i] = vy[i] - 1 * dv*normal[1];



		error = fmax(fabs(dv), error);

		dvBoundary[i] = dvBoundary[i] + dv;

	}

	return error;

}

# endif
# ifndef NOUPDATE

void femGrainsUpdate(femGrains *myGrains, double dt, double tol, double iterMax)
{
	int n = myGrains->n;
	int i,iter = 0;
	double error;
	double *x          = myGrains->x;
	double *y          = myGrains->y;
	double *m          = myGrains->m;
	double *vy         = myGrains->vy;
	double *vx         = myGrains->vx;
	double gamma       = myGrains->gamma;
	double gx          = myGrains->gravity[0];
	double gy          = myGrains->gravity[1];

	// 
	// // -1- Calcul des nouvelles vitesses des grains sur base de la gravité et de la trainee
	// //
	
	double fact;
	int j;
	for (j = 0;j < n;j++)
	{
		fact = 1 - gamma*dt / (m[j]);
		vx[j] = fact * vx[j] + gx * dt;
		vy[j] = fact * vy[j] + gy * dt;
	}

	//
	// -2- Correction des vitesses pour tenir compte des contacts        
	//       

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
