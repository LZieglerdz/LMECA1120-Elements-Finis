#include"stdlib.h"
#include"fem.h"
#include"math.h"

femMesh *globalMesh;

int xComp(const void* e0, const void *e1) {
	int *left = ((int *) e0);
	int *right = ((int *) e1);
	double diagnostic = globalMesh->X[*left] - globalMesh->X[*right];
	if (diagnostic < 0) {
		return 1;
	}
	if (diagnostic > 0) {
		return -1;
	}
	else return 0;
} 

int yComp(const void* e0, const void *e1) {
	int *left = ((int *) e0);
	int *right = ((int *) e1);
	double diagnostic = globalMesh->Y[*left] - globalMesh->Y[*right];
	if (diagnostic < 0) {
		return 1;
	}
	if (diagnostic > 0) {
		return -1;
	}	
	else return 0;
}

void femDiffusionRenumber(femDiffusionProblem *theProblem, femRenumType renumType)
{
    int i;
    int nNode = theProblem->mesh->nNode;
    int buffer[nNode];
    
    switch (renumType) {
        case FEM_NO :
            for (i = 0; i < nNode; i++) {
                theProblem->number[i] = i;
	    }    
	break;
        case FEM_XNUM :
		for(i = 0; i < nNode; i++) {
			buffer[i] = i;
		}
		globalMesh = theProblem->mesh;
		qsort(buffer, nNode, sizeof(int), xComp);
		for (i = 0; i < nNode; i++) {
			theProblem->number[buffer[i]] = i;
		}
	break;
        case FEM_YNUM : 
		for(i = 0; i < nNode; i++) {
			        buffer[i] = i;
		}
		globalMesh = theProblem->mesh;
		qsort(buffer, nNode, sizeof(int), yComp);
		for (i = 0; i < nNode; i++) {
			        theProblem->number[buffer[i]] = i;
		}
	break;            
        default : Error("Unexpected renumbering option"); }
}

int femDiffusionComputeBand(femDiffusionProblem *theProblem) {

	femMesh *theMesh = theProblem->mesh;
	int i,j, mapMax, mapMin;
	int nLocNode 	= theMesh->nLocalNode;
	int nElem 	= theMesh->nElem;
	int myBand 	= 0;
	int map[(1 + nElem) * nLocNode];
	for (i = 0; i < nElem; i++) {
		for (j = 0; j < nLocNode; j++) {
			map[j] = theProblem->number[theMesh->elem[i*nLocNode + j]];
		}
		mapMax = map[0];
	    	mapMin = map[0];
//		printf("myBand mapMax mapMin diff  %.2d %.2d %.2d %.2d \n", myBand, mapMax, mapMin, abs(mapMax-mapMin));

		for (j = 1; j < nLocNode; j++) {
			mapMax = fmax(mapMax, map[j]);
			mapMin = fmin(mapMin, map[j]);
		}
    		myBand = fmax(myBand, fabs(mapMax - mapMin));	
//		printf("%.2d\n", myBand);
	}
	return(myBand + 1);
}

void femIterativeSolverAssemble(femIterativeSolver* mySolver, double *Aloc, double *Bloc, double *Uloc, int *map, int nLoc) {
   
       	
	int i,j;
	for (i = 0; i < nLoc; i++) { 
		int myRow = map[i];
		for(j = 0; j < nLoc; j++) {       
			mySolver->R[myRow] -= Aloc[i*nLoc+j]*Uloc[j];
		}
		mySolver->R[myRow] += Bloc[i];		//https://en.wikipedia.org/wiki/Band_matrix#Band_storage
	}
	for (i = 0; i < nLoc; i++){
		int row =  map[i];
		for (j = 0; j < nLoc; j++) {
			int col = map[j];
			mySolver->S[row] -= Aloc[i*nLoc+j] * mySolver->D[col];
	}
   } 



}

void femIterativeSolverConstrain(femIterativeSolver* mySolver, int myNode, double myValue) {   
	mySolver->X[myNode] = 0;
	mySolver->D[myNode] = 0;
	mySolver->R[myNode] = 0; 
	mySolver->S[myNode] = 0;
}

double scalar(double array1[], double array2[], int size) {
	int i;
	double scalr = 0;
	for (i = 0; i < size; i++) {
		        scalr += array1[i] * array2[i];
	}
	return scalr;
}


double *femIterativeSolverEliminate(femIterativeSolver *mySolver) {

    mySolver->iter++;
    int i, size = mySolver->size;

    
    double error = 0.0;
    for (i = 0; i < size; i++) {
	mySolver->S[i] = 0;
        error = (mySolver->R[i])*(mySolver->R[i]);
    }
    
    if (mySolver->iter==1) {
	for (i = 0; i < size; i++) {
		mySolver->X[i] = 0;
		mySolver->D[i] = mySolver->R[i];
	}
    }
    else {
	double alpha = -error / scalar(mySolver->R, mySolver->S, size);
	for (i = 0; i < size; i++) {
		mySolver->R[i] += alpha * mySolver->S[i];
	}
	double beta = scalar(mySolver->R, mySolver->R, size) / error;
	for (i = 0; i < size; i++) {
		mySolver->X[i] += alpha * mySolver->D[i];
		mySolver->D[i] = mySolver->R[i] + beta * mySolver->D[i];
//		mySolver->S[i] = 0;
	}

    }
     
    mySolver->error = sqrt(error);
/*
	double alpha = - ( scalar(mySolver->R, mySolver->R, size) / scalar(mySolver->S, mySolver->R, size) );
	double rk1[size];
	for (i = 0; i < size; i++) {
		rk1[i] = mySolver->R[i] + alpha * mySolver->S[i];		
	}
	double beta = scalar(rk1, rk1, size) / scalar(mySolver->R, mySolver->R, size);

	for (i = 0; i < size; i++) {
		mySolver->X[i] += alpha*mySolver->D[i];
	}
	for (i = 0; i < size; i++) {
		mySolver->D[i] = rk1[i] + beta*mySolver->D[i];
	}
	for (i = 0; i < size; i++) {
		mySolver->R[i] = rk1[i];
	}

*/

	return(mySolver->R);

}
