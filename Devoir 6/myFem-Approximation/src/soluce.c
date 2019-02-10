
#include "fem.h"


# ifndef NOAPPROXPHI

void femApproxPhi(double xsi, double eta, double *phi)
{
    phi[0] = (1.0 - xsi - eta) * (1.0/3.0 - xsi - eta) * (2.0/3.0 - xsi - eta) *  9.0/2.0; 
    phi[1] =               xsi * (xsi - 1.0/3.0)       * (xsi - 2.0/3.0)       *  9.0/2.0;  
    phi[2] =               eta * (eta - 1.0/3.0)       * (eta - 2.0/3.0)       *  9.0/2.0;  
    phi[3] = (1.0 - xsi - eta) * (2.0/3.0 - xsi - eta) * xsi                   * 27.0/2.0;
    phi[4] = (1.0 - xsi - eta) * (xsi - 1.0/3.0)       * xsi                   * 27.0/2.0;
    phi[5] =               xsi * eta                   * (xsi - 1.0/3.0)       * 27.0/2.0;;
    phi[6] =               xsi * eta                   * (eta - 1.0/3.0)       * 27.0/2.0;  
    phi[7] = (1.0 - xsi - eta) * (eta - 1.0/3.0)       * eta                   * 27.0/2.0;
    phi[8] = (1.0 - xsi - eta) * (2.0/3.0 - xsi - eta) * eta                   * 27.0/2.0;
    phi[9] = (1.0 - xsi - eta) * xsi                   * eta                   * 27.0;   
}

# endif
# ifndef NOAPPROXDPHI

void femApproxDphi(double xsi, double eta, double *dphidxsi, double *dphideta)
{
    int i;
    double phi[10];
    double phi_xsi[10];
    double phi_eta[10];
    double epsilon = 1e-10;
    femApproxPhi(xsi,eta,phi);
    femApproxPhi(xsi + epsilon,eta,phi_xsi);
    femApproxPhi(xsi,eta + epsilon,phi_eta);
    for (i=0;i < 10; i++) {
      dphidxsi[i] = (phi_xsi[i] - phi[i])/epsilon;
      dphideta[i] = (phi_eta[i] - phi[i])/epsilon;  }
      
}

# endif
# ifndef NOAPPROXLOCAL

int *mapGlobal = NULL;

void femApproxLocal(const femApproxProblem *theProblem, const int iElem, int *map)
{
          
  // Calcul de la numerotation globale de tous les elements lors du premier appel
  // de la fonction....  Ce calcul ne sera donc effectue qu'UNE UNIQUE fois !
  // Evidemment, ce gain necssite l'allocation d'un tableau en memoire.
  // CPU versus MEMORY : l'eternel compromis impossible :-)
  //
  
    if (mapGlobal == NULL) {
        printf("Allocating mapGlobal for all elements...\n");  
        const int segment[6][4] = {{0,1,3,4},{1,0,4,3},{1,2,5,6},{2,1,6,5},{2,0,7,8},{0,2,8,7}};
        femMesh *theMesh = theProblem->mesh;
        int nNode = theMesh->nNode;
        int nElem = theMesh->nElem;
        femEdges *theEdges = theProblem->edges;
        int nEdge = theEdges->nEdge;        
        mapGlobal = malloc(sizeof(int)*nElem*10);        
        for (int jElem=0; jElem < theMesh->nElem; jElem++) {
            int *theMapGlobal = &mapGlobal[jElem*10];            
            int *node = &theMesh->elem[jElem*3];
            for (int j=0; j < 3; j++) 
                theMapGlobal[j] = node[j];        
            for (int i=0; i < nEdge; i++) {
                femEdge theEdge = theEdges->edges[i];
                for (int j=0; j < 6; j++) {
                    if (theEdge.node[0] == node[segment[j][0]] 
                                  && theEdge.node[1] == node[segment[j][1]]) {
                        theMapGlobal[segment[j][2]] = nNode + i*2    ;
                        theMapGlobal[segment[j][3]] = nNode + i*2 + 1;  }}}    
            theMapGlobal[9] = nNode + 2*nEdge + jElem;  }}
        
  // Assigne la bonne valeur au vecteur :-)
  // Il serait plus efficace de simplement envoyer le pointeur vers la case
  // memoire du tableau, mais cela necessiterait de changer la signature
  // de la fonction.
   
    int *theMapGlobal = &mapGlobal[iElem*10];
    for (int j=0; j < 10; j++) 
        map[j] = theMapGlobal[j]; 
        
  // Pour une approximation P1-CO, remplacer ce qui precede par :
  //
  //  femMesh *theMesh = theProblem->mesh;
  //  for (int j=0; j < 3; j++) {  
  //      map[j] = theMesh->elem[iElem*3 + j]; }
     

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
            y[j]   = theMesh->Y[elem[j]]; } 
          
        for (iInteg=0; iInteg < theRule->n; iInteg++) {    
            double xsi    = theRule->xsi[iInteg];
            double eta    = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];
            
            femDiscretePhi2(theSpace,xsi,eta,phi);
            double xloc = x[2]*eta + x[1]*xsi + x[0]*(1.0-xsi-eta);
            double yloc = y[2]*eta + y[1]*xsi + y[0]*(1.0-xsi-eta);            
            double uloc = femApproxStommel(xloc,yloc);
            double jac = (x[1]-x[0])*(y[2]-y[0])-(y[1]-y[0])*(x[2]-x[0]);
             
            for (i = 0; i < theSpace->n; i++) { 
                for(j = 0; j < theSpace->n; j++) {
                    theSystem->A[map[i]][map[j]] += phi[i] * phi[j] * jac * weight; }}                                                                                            
            for (i = 0; i < theSpace->n; i++) {
                theSystem->B[map[i]] += phi[i] * uloc * jac *weight; }}} 
                
    femFullSystemEliminate(theSystem);
}


# endif
