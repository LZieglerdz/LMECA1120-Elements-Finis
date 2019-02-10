//
//Authors: Ziegler Laurent & Lacroix Arthur
//Base sur le cours LMECA1120:Introduction aux éléments finis donné par le Pr. Vincent Legat
//


#include <stdio.h>
#include <math.h>


#ifdef graphic
#include "glfem.h"
#endif

double integrate(double x[3], double y[3], double (*f) (double, double)){
	int i;
	double I = 0;
	double xLoc[3];
	double yLoc[3];

	double xik[3] = {.166666667,.166666667,.666666667};
	double etak[3] = {.166666667,.666666667,.166666667};
	double wk[3] = {.166666667,.166666667,.166666667};

	double Jack = fabs((x[1]-x[0])*(y[2]-y[0]) - (x[2]-x[0])*(y[1]-y[0]));


	for (i=0;i<3;i++) {
		double xiLoc = xik[i];
		double etaLoc = etak[i];
		double phi[3] = {1-xiLoc-etaLoc, xiLoc, etaLoc};
		xLoc[i] = x[0] * phi[0] + x[1] * phi[1] + x[2] * phi[2];
		yLoc[i] = y[0] * phi[0] + y[1] * phi[1] + y[2] * phi[2];
		I += f(xLoc[i],yLoc[i]) * wk[i];
	}

#ifdef graphic

    glColor3f (1.0,1.0,1.0); glfemDrawElement(x,y,3);
    glColor3f (1.0,0.0,0.0); glfemDrawNodes(x,y,3);
    glColor3f (0.0,0.0,1.0); glfemDrawNodes(xLoc,yLoc,3);
    
#endif

    return Jack * I;
}

double integrateRecursive(double x[3], double y[3], double (*f)(double,double), int n)
{
	double I=0;

	if (n == 0) {
		I += integrate(x,y,f);
		return I;
	}
	else {	
		// new nodes
		double N0[2] = { x[0], y[0] };
		double N1[2] = { x[1], y[1] };
		double N2[2] = { x[2], y[2] };
		double N3[2] = { (x[0] + x[1])/2, (y[0]+y[1])/2 };
		double N4[2] = { (x[0] + x[2])/2, (y[0]+y[2])/2 };
		double N5[2] = { (x[1] + x[2])/2, (y[1]+y[2])/2 };

		//new triangles
		double x1[3] = {N0[0], N3[0], N4[0]};
		double y1[3] = {N0[1], N3[1], N4[1]};
		double x2[3] = {N3[0], N1[0], N5[0]};
		double y2[3] = {N3[1], N1[1], N5[1]};
		double x3[3] = {N3[0], N5[0], N4[0]};
		double y3[3] = {N3[1], N5[1], N4[1]};
		double x4[3] = {N4[0], N5[0], N2[0]};
		double y4[3] = {N4[1], N5[1], N2[1]};

		I += integrateRecursive(x1,y1,f,n-1) + integrateRecursive(x2,y2,f,n-1) + integrateRecursive(x3,y3,f,n-1) + integrateRecursive(x4,y4,f,n-1);

		return I;
	}
}
