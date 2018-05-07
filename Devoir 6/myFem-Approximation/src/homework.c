#include "math.h"
#include "fem.h"


  

# ifndef NOAPPROXPHI

void femApproxPhi(double xsi, double eta, double *phi)
{
    phi[0] = (1.0 - xsi - eta) * (1.0/3.0 - xsi - eta) * (2.0/3.0 - xsi - eta) *  9.0/2.0; 
    phi[1] = xsi * (xsi - 1.0/3.0) * (xsi - 2.0/3.0) * 9.0/2.0;
    phi[2] = eta * (eta - 1.0/3.0) * (eta - 2.0/3.0) * 9.0/2.0;
    phi[3] = xsi * (1.0 - xsi - eta) * (2.0/3.0 - xsi - eta) * 27.0/2.0;
    phi[4] = xsi * (xsi - 1.0/3.0) * (1.0 - xsi - eta) * 27.0/2.0;
    phi[5] = xsi * eta * (xsi - 1.0/3.0) * 27.0/2.0; 
    phi[6] = xsi * eta * (eta - 1.0/3.0) * 27.0/2.0;		//ATTENTION NUMEROTATION NOEUDS SUR S
    phi[7] = eta * (1.0 - xsi - eta) * (eta - 1.0/3.0) * 27.0/2.0;		//SEGMENTS
    phi[8] = eta * (1.0 - xsi - eta) * (2.0/3.0 - xsi - eta) * 27.0/2.0;
    phi[9] = xsi * eta * (1.0 - xsi - eta) * 27.0;
    
    // La premiere fonction vous est donnee gracieusement....
    // Il vous reste les 9 autres a obtenir !
    // A modifier :-)

}

# endif
# ifndef NOAPPROXDPHI

void femApproxDphi(double xsi, double eta, double *dphidxsi, double *dphideta)
{
	dphidxsi[0] = (-27*pow(xsi,2) + xsi*(36 - 54*eta) - 27*pow(eta,2) + 36*eta - 11)/2;
	dphideta[0] = (-27*pow(xsi,2) + xsi*(36 - 54*eta) - 27*pow(eta,2) + 36*eta - 11)/2;

	dphidxsi[1] = 1.0 - 9.0*xsi + 13.5*pow(xsi,2);
	dphideta[1] = 0;

	dphidxsi[2] = 0;
	dphideta[2] = 1.0 - 9.0*eta + 13.5*pow(eta,2);

	dphidxsi[3] =9.0 + 40.5*pow(xsi,2) - 22.5*eta + 13.5*pow(eta,2) + xsi*(-45.0 + 54.0*eta) ;
	dphideta[3] = xsi*(-22.5 + 27.0*xsi + 27.0*eta);
	
	dphidxsi[4] = -4.5 - 40.5*pow(xsi,2) + xsi*(36.0 - 27.0*eta) + 4.5*eta;
	dphideta[4] = - 13.5*(-(1.0/3.0) + xsi)*xsi;
	
	dphidxsi[5] = (-4.5 + 27.0*xsi)*eta;
	dphideta[5] = 13.5*(-(1.0/3.0) + xsi)*xsi;
	
	dphidxsi[6] = 13.5*(-(1.0/3.0) + eta)*eta;
	dphideta[6] = xsi*(-4.5 + 27.0*eta);
	
	dphidxsi[7] = -13.5* (eta-1/3)*eta;
	dphideta[7] = 4.5*xsi-40.5*pow(eta,2)+(36-27*xsi)*eta - 4.5;		
	
	dphideta[8] = eta*(-22.5 + 27.0*xsi + 27.0*eta);
	dphidxsi[8] = 13.5*pow(xsi,2) - 22.5*xsi+40.5*pow(eta,2)+(54*xsi-45)*eta+9;
		
	dphidxsi[9] = 2*(13.5 - 27.0*xsi - 13.5*eta)*eta;
	dphideta[9] = 2*xsi*(13.5 - 13.5*xsi - 27.0*eta);





    // A faire a titre de bonus, car on n'a pas vraiment besoin de ces derivees...
    // Il n'y a qu'un unique point a gagner avec ceci
}

# endif
# ifndef NOAPPROXLOCAL

void femApproxLocal(const femApproxProblem *theProblem, const int iElem, int *map)
{
    
    int i, j;
    // Pour une approximation P1-CO, remplacer ce qui suit par :
    //
    //  femMesh *theMesh = theProblem->mesh;
    //  for (int j=0; j < 3; j++) {  
    //      map[j] = theMesh->elem[iElem*3 + j]; }
    // 
    if (iElem > 1)  Error("It only works with meshes one.txt and two.txt ");
    int mapElem[2][10] = {{0,1,2,3,4,5,6,7,8,9},{1,10,2,11,12,13,14,6,5,15}};
    for (j=0; j < 10; j++) { 
        map[j] = mapElem[iElem][j]; }  


    // A generaliser pour un maillage quelconque avec une interpolation P3-C0
    // Avec la version actuelle, le code ne fonctionne qu'avec les maillages one.txt et two.text
    //
    // Reflechir aussi pour obtenir une version efficace de cette fonction qui
    // est appelee TRES souvent lors de l'execution du code :-)
  
       
}

# endif
# ifndef NOAPPROXSOLVE

void femApproxSolve(femApproxProblem *theProblem)
{
	
    femMesh *theMesh = theProblem->mesh;
    femFullSystem *theSystem = theProblem->system;
    femIntegration *theRule = theProblem->rule;
    femDiscrete *theSpace = theProblem->space;
 
    double x[3],y[3],phi[10];
    int iElem,iInteg,i,j,map[10];

    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        femApproxLocal(theProblem,iElem,map); 
        for (j=0; j < 3; ++j) {
            int *elem = &theMesh->elem[iElem*3];
            x[j]   = theMesh->X[elem[j]];
            y[j]   = theMesh->Y[elem[j]];
	}           
        for (iInteg=0; iInteg < theRule->n; iInteg++) {    
		double xsi	= theRule->xsi[iInteg];
		double eta	= theRule->eta[iInteg];
		double weight	= theRule->weight[iInteg];

		double dxdxsi   = 0;
		double dydxsi   = 0;
		double dxdeta  = 0;
		double dydeta  = 0;

		double phiLoc[3]  = {0,0,0};
		double dphi[2][3] = {0,0,0,0,0,0} ;
		_p1c0_phi(theRule->xsi[iInteg], theRule->eta[iInteg], phiLoc);
		_p1c0_dphidx(theRule->xsi[iInteg], theRule->eta[iInteg], dphi[0], dphi[1]);


		double coord[2]	= {0,0};
		for (i = 0; i < 3; i++) {
			dxdxsi   += x[i]*dphi[0][i];
			dydxsi   += y[i]*dphi[0][i];
			dxdeta	 += x[i]*dphi[1][i];
			dydeta   += y[i]*dphi[1][i];
			coord[0] += x[i]*phiLoc[i];
			coord[1] += y[i]*phiLoc[i];
		}		
		double jacobian = fabs(dxdxsi * dydeta - dydxsi * dxdeta);
//		printf("%.2f ", jacobian);		

		femApproxPhi(theRule->xsi[iInteg], theRule->eta[iInteg], phi);
		for (i = 0; i < theSpace->n; i++) {
	            for (j = 0; j < theSpace->n; j++) { 
        	        theSystem->A[map[i]][map[j]] += phi[i] * phi[j] * weight * jacobian;
/*			if (j < theSpace->n-1){
				printf("%.2f ", theSystem->A[map[i]][map[j]]);
			}
			else{
				printf("%.2f \n", theSystem->A[map[i]][map[j]]);
			}
*/		    }
                  theSystem->B[map[i]] += phi[i] * femApproxStommel(coord[0], coord[1]) * weight * jacobian;
		}
	}
    }      
    femFullSystemEliminate(theSystem);
}


# endif

